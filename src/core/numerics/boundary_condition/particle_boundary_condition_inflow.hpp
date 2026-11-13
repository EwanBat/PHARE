#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/utilities/random.hpp"
#include "core/def.hpp"

#include <array>
#include <bit>
#include <cstdint>

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
 *     Init     draw the stock in the box
 *     Refill   empty the box, draw it again (the push left a gap at the far edge)
 *     Moments  deposit the stock                              (not written yet)
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

    /** @brief Inflow acts at Init, Refill and Moments. */
    bool actsAt(ParticleBoundaryStage const stage) const override
    {
        return stage == ParticleBoundaryStage::Init or stage == ParticleBoundaryStage::Refill
               or stage == ParticleBoundaryStage::Moments;
    }

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
            case ParticleBoundaryStage::Moments:
                // TODO: deposit the stock, only the part in front of this patch's own cells.
                break;
            case ParticleBoundaryStage::Leavers: break; // inflow does not act then
        }
    }

private:
    /** @brief Draw the inflow Maxwellian in every cell of @p box.
     *
     * The time is mixed into the key: without it, each Refill would draw exactly the same
     * particles in a cell, and the same pattern would enter at every step. */
    void drawStock_(ParticleArrayT& stock, Box<int, Super::dimension> const& box,
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
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
