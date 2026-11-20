#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/utilities/random.hpp"
#include "core/def.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <ranges>
#include <stdexcept>
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
 *     PredictedRefill   draw the stock of the new time in the box, in a temp stock (the real one
 *                       must stay at the start of the step for the real push), kept for Refill
 *     PredictedMoments  deposit the part of the temp stock in front of the patch's own cells
 *     Refill            empty the box, fill it again (the push left a gap at the far edge): with
 *                       the temp stock if it is for this box and this time, else with a new draw
 *     Moments           deposit the part of the stock in front of the patch's own cells
 *
 * The two parts of the step match:
 *
 *                       refill                      deposit
 *     push of copies    PredictedRefill → temp      PredictedMoments ← temp
 *     real push         Refill ← temp → stock       Moments ← stock
 *     Leavers           erase the strays: stock particles the real push took into a neighbour's
 *                       cells. The neighbour drew the same particle (same key), and took it into
 *                       its domain or keeps it in its own front; this copy would be counted twice
 *                       if it came back (see strayBoxes)
 *
 * Why a temp stock after the push of copies (domain_only): the real stock is still at the start of
 * the step, so depositing it would mix two times and count twice the particles that came in (the
 * updater already deposited their pushed copies). Only the part in front of the patch's own cells
 * is deposited (frontOfOwnCells).
 *
 * Why Refill reuses it: draws are keyed by (seed, time, population, cell), and both stages run at
 * the same new time, so a Refill draw would give exactly the particles PredictedRefill drew, in
 * the same order. Drawing is most of the cost of this condition, so it is done once per step.
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
            case ParticleBoundaryStage::Init:
                kept_.clear(); // the level is new: a kept stock may be for a box gone elsewhere
                drawStock_(stock, ghostBox, gridLayout, time);
                break;
            case ParticleBoundaryStage::PredictedRefill:
            {
                // the real stock is still at the start of the step: draw the one of the new time,
                // in the whole box, in the temp stock, and keep it for the Refill at the same time
                auto kept = findKept_(ghostBox);
                if (kept == kept_.end())
                    kept = kept_.insert(kept_.end(), Kept{ghostBox, time, {}});
                kept->time = time;
                kept->particles.clear();
                drawStock_(kept->particles, ghostBox, gridLayout, time);
                break;
            }
            case ParticleBoundaryStage::PredictedMoments:
            {
                auto const kept = findKept_(ghostBox);
                if (kept == kept_.end() or kept->time != time)
                    throw std::runtime_error(
                        "Inflow boundary: PredictedMoments without a PredictedRefill before it.");

                // deposited in place, without a copy: the temp stock has no cell index, so it is
                // filtered particle by particle
                auto const ownBox = frontOfOwnCells(ghostBox, gridLayout.AMRBox());
                auto inFront      = kept->particles | std::views::filter([&](auto const& p) {
                                   return isIn(Point<int, Super::dimension>{p.iCell}, ownBox);
                               });
                interpolator_(inFront, data.density, data.flux, gridLayout);
                break;
            }
            case ParticleBoundaryStage::Refill:
            {
                stock.erase_in(ghostBox); // visits only the cells of ghostBox
                auto const kept = findKept_(ghostBox);
                if (kept != kept_.end() and kept->time == time)
                {
                    stock.reserve(stock.size() + kept->particles.size());
                    for (auto const& particle : kept->particles)
                        stock.push_back(particle);
                    kept->time = std::numeric_limits<double>::quiet_NaN(); // used: never again
                }
                else
                    drawStock_(stock, ghostBox, gridLayout, time);
                break;
            }
            case ParticleBoundaryStage::Moments:
            {
                // the stock has a cell index: only the cells of ownBox are visited, copied out
                auto const ownBox = frontOfOwnCells(ghostBox, gridLayout.AMRBox());
                toDeposit_.clear();
                stock.export_particles(ownBox, toDeposit_, [](auto const& p) { return p; });
                interpolator_(toDeposit_, data.density, data.flux, gridLayout);
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

    //! the particles the Moments stage deposits, kept between calls to reuse its memory
    std::vector<particle_type> toDeposit_;

    /** A temp stock drawn by PredictedRefill, deposited by PredictedMoments, used by the Refill
     * of the same box at the same time. One condition serves every patch of its face on this
     * rank (and a corner box it rules), so it keeps one per ghost box. Refill uses it once (its time becomes NaN, its memory stays for
     * the next step); Init clears them all. */
    struct Kept
    {
        Box<int, Super::dimension> box;
        double time;
        std::vector<particle_type> particles;
    };
    std::vector<Kept> kept_; //!< a few entries: one per ghost box of this face on this rank

    auto findKept_(Box<int, Super::dimension> const& box)
    {
        return std::find_if(kept_.begin(), kept_.end(), [&](Kept const& k) { return k.box == box; });
    }

    interpolator_type interpolator_;
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
