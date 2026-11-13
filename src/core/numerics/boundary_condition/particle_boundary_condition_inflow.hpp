#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/def.hpp"

#include <array>

namespace PHARE::core
{
/**
 * @brief Particle boundary condition for an inflow boundary.
 *
 * Particles of a prescribed plasma are drawn in the ghost cells beyond the boundary, and those
 * that cross it during a push enter the domain.
 *
 * It holds the inflow plasma given by the user: density, bulk velocity and thermal velocity.
 *
 * Physics not implemented yet: the moments are chosen, the work is not written.
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
     */
    ParticleBoundaryConditionInflow(double const density, std::array<double, 3> const& velocity,
                                    std::array<double, 3> const& thermalVelocity)
        : density_{density}
        , velocity_{velocity}
        , thermalVelocity_{thermalVelocity}
    {
    }

    NO_DISCARD double density() const { return density_; }
    NO_DISCARD std::array<double, 3> const& velocity() const { return velocity_; }
    NO_DISCARD std::array<double, 3> const& thermalVelocity() const { return thermalVelocity_; }

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

    void apply(ParticleBoundaryStage const stage, data_type& /*data*/,
               BoundaryLocation const /*location*/,
               Box<int, Super::dimension> const& /*ghostBox*/, GridLayoutT const& /*gridLayout*/,
               double const /*time*/) override
    {
        switch (stage)
        {
            case ParticleBoundaryStage::Init:
                // TODO: fill the box with particles drawn from the inflow plasma.
                break;
            case ParticleBoundaryStage::Refill:
                // TODO: the push left an empty strip at the far edge of the stock.
                // Remove the stock in the box and draw it again, full.
                break;
            case ParticleBoundaryStage::Moments:
                // TODO: deposit the stock, only the part in front of this patch's own cells.
                break;
            case ParticleBoundaryStage::Leavers: break; // inflow does not act then
        }
    }

private:
    double density_;
    std::array<double, 3> velocity_;
    std::array<double, 3> thermalVelocity_;
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
