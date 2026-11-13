#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

namespace PHARE::core
{
/**
 * @brief Particle boundary condition for an open (outflow) boundary.
 *
 * Particles that leave the domain through the boundary are kept in the ghost cells beyond it, and
 * deposited there, until they leave the ghost layer.
 *
 * Physics not implemented yet: the moments are chosen, the work is not written.
 */
template<typename ParticleArrayT, typename FieldT, typename GridLayoutT>
class ParticleBoundaryConditionOpen : public IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>
{
public:
    using Super     = IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>;
    using data_type = Super::data_type;

    /** @brief Open acts at Moments and Leavers. */
    bool actsAt(ParticleBoundaryStage const stage) const override
    {
        return stage == ParticleBoundaryStage::Moments or stage == ParticleBoundaryStage::Leavers;
    }

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Open;
    }

    void apply(ParticleBoundaryStage const stage, data_type& /*data*/,
               BoundaryLocation const /*location*/,
               Box<int, Super::dimension> const& /*ghostBox*/, GridLayoutT const& /*gridLayout*/,
               double const /*time*/) override
    {
        switch (stage)
        {
            case ParticleBoundaryStage::Moments:
                // TODO: deposit the stock (the particles that left in earlier steps).
                break;
            case ParticleBoundaryStage::Leavers:
                // TODO: move the particles that just left through the face from
                // patchGhostParticles to levelGhostParticles (the stock).
                break;
            default: break; // open does not act at the other moments
        }
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP
