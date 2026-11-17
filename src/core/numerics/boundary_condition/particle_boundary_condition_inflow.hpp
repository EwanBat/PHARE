#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/utilities/random.hpp"
#include "core/def.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <vector>

namespace PHARE::core
{
/**
 * @brief Particle boundary condition for an inflow boundary.
 *
 * Particles of a prescribed plasma are drawn in the ghost cells beyond the boundary (the
 * "stock"), and those that cross it during a push enter the domain.
 *
 * It holds the inflow plasma given by the user (density, bulk velocity, thermal velocity) and,
 * from the population, the number of particles per cell and the random seed.
 *
 *     Init              draw the stock in the box
 *     PredictedMoments  draw a temp stock at the new time, deposit it, drop it
 *     Refill            empty the box, draw it again (the push left a gap at the far edge)
 *     Moments           deposit the stock
 *     Leavers           erase the strays: stock particles the real push took into a neighbour's
 *                       cells. The neighbour drew the same particle (same key), and took it into
 *                       its domain or keeps it in its own front; this copy would be counted twice
 *                       if it came back (see strayBoxes)
 *
 * Why a temp stock at PredictedMoments: after the push of copies (domain_only), the real stock is
 * still at the start of the step, so depositing it would mix two times and count twice the
 * particles that came in (the updater already deposited their pushed copies). The temp draw uses
 * the same key as the Refill that follows at the same time, so both stages deposit the same
 * particles. Only the part in front of the patch's own cells is deposited (frontOfOwnCells).
 */
template<typename ParticleArrayT, typename FieldT, typename GridLayoutT>
class ParticleBoundaryConditionInflow
    : public IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>
{
public:
    using Super     = IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>;
    using data_type = Super::data_type;

    /**
     * @param density inflow density.
     * @param velocity inflow bulk velocity (x, y, z).
     * @param thermalVelocity inflow thermal velocity (x, y, z).
     * @param population the population: name, particles per cell, seed.
     */
    ParticleBoundaryConditionInflow(double const density, std::array<double, 3> const& velocity,
                                    std::array<double, 3> const& thermalVelocity,
                                    ParticleBoundaryInfo const& population)
        : density_{density}
        , velocity_{velocity}
        , thermalVelocity_{thermalVelocity}
        , nbrParticlesPerCell_{population.nbrParticlesPerCell}
        , seed_{MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>::baseSeed(
              population.seed)}
        , populationKey_{fnv1a64(population.name)}
    {
        if (nbrParticlesPerCell_ == 0)
            throw std::runtime_error("Inflow boundary: population '" + population.name
                                     + "' has no particles per cell.");
    }

    NO_DISCARD double density() const { return density_; }
    NO_DISCARD std::array<double, 3> const& velocity() const { return velocity_; }
    NO_DISCARD std::array<double, 3> const& thermalVelocity() const { return thermalVelocity_; }
    NO_DISCARD std::uint32_t nbrParticlesPerCell() const { return nbrParticlesPerCell_; }

    /** @brief Inflow acts at every stage. */
    bool actsAt(ParticleBoundaryStage const /*stage*/) const override { return true; }

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Inflow;
    }

    void apply(ParticleBoundaryStage const stage, data_type& data,
               BoundaryLocation const /*location*/, Box<int, Super::dimension> const& ghostBox,
               GridLayoutT const& gridLayout, double const time) override
    {
        auto& stock = data.levelGhostParticles;

        switch (stage)
        {
            case ParticleBoundaryStage::Init: drawStock_(stock, ghostBox, gridLayout, time); break;
            case ParticleBoundaryStage::Refill:
                stock.erase_in(ghostBox); // visits only the cells of ghostBox
                drawStock_(stock, ghostBox, gridLayout, time);
                break;
            case ParticleBoundaryStage::PredictedMoments:
            {
                // the real stock is still at the start of the step: draw the one of the new time
                auto const ownBox = frontOfOwnCells(ghostBox, gridLayout.AMRBox());
                toDeposit_.clear();
                drawStock_(toDeposit_, ownBox, gridLayout, time);
                deposit_(data, gridLayout);
                break;
            }
            case ParticleBoundaryStage::Moments:
            {
                auto const ownBox = frontOfOwnCells(ghostBox, gridLayout.AMRBox());
                toDeposit_.clear();
                stock.export_particles(ownBox, toDeposit_, [](auto const& p) { return p; });
                deposit_(data, gridLayout);
                break;
            }
            case ParticleBoundaryStage::Leavers:
                for (auto const& strays : strayBoxes(ghostBox, gridLayout.AMRBox()))
                    stock.erase_in(strays);
                break;
        }
    }

private:
    using particle_type = ParticleArrayT::value_type;
    using interpolator_type
        = Interpolator<Super::dimension, GridLayoutT::options.interp_order>;

    /** @brief Add the particles of toDeposit_ to the density and flux of the population. */
    void deposit_(data_type& data, GridLayoutT const& gridLayout)
    {
        interpolator_(toDeposit_, data.density, data.flux, gridLayout);
    }

    /** @brief Draw the inflow Maxwellian in every cell of @p box, into @p stock (a ParticleArray
     * or a std::vector of particles).
     *
     * The time is mixed into the key: without it, each Refill would draw exactly the same
     * particles in a cell, and the same pattern would enter at every step. */
    template<typename Particles>
    void drawStock_(Particles& stock, Box<int, Super::dimension> const& box,
                    GridLayoutT const& gridLayout, double const time) const
    {
        auto const seed = hash_combine(seed_, std::bit_cast<std::uint64_t>(time));
        loadMaxwellianInBox(stock, box, density_, velocity_, thermalVelocity_,
                            nbrParticlesPerCell_, seed, populationKey_, gridLayout.levelNumber());
    }

    double density_;
    std::array<double, 3> velocity_;
    std::array<double, 3> thermalVelocity_;
    std::uint32_t nbrParticlesPerCell_;
    std::uint64_t seed_;          //!< key base, drawn once if the population has no seed
    std::uint64_t populationKey_; //!< two populations must not get the same particles

    //! the particles a Moments stage deposits, kept between calls to reuse its memory
    std::vector<particle_type> toDeposit_;
    interpolator_type interpolator_;
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
