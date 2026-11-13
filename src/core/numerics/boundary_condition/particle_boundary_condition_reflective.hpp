#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_REFLECTIVE_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_REFLECTIVE_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

namespace PHARE::core
{
/**
 * @brief Particle boundary condition for a reflective (perfectly conducting) boundary.
 *
 * Particles that cross the boundary are sent back into the domain, their normal velocity
 * reversed.
 *
 * Physics not implemented yet: the moments are chosen, the work is not written.
 */
template<typename ParticleArrayT, typename FieldT, typename GridLayoutT>
class ParticleBoundaryConditionReflective
    : public IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>
{
public:
    using Super     = IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>;
    using data_type = Super::data_type;

    /** @brief Reflective acts at Moments and Leavers. */
    bool actsAt(ParticleBoundaryStage const stage) const override
    {
        return stage == ParticleBoundaryStage::Moments or stage == ParticleBoundaryStage::Leavers;
    }

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Reflective;
    }

    void apply(ParticleBoundaryStage const stage, data_type& /*data*/,
               BoundaryLocation const /*location*/,
               Box<int, Super::dimension> const& /*ghostBox*/, GridLayoutT const& /*gridLayout*/,
               double const /*time*/) override
    {
        switch (stage)
        {
            case ParticleBoundaryStage::Moments:
                // TODO: make density and flux symmetric about the face: add each ghost node
                // to its mirror node inside, then copy the result back to the ghost node.
                break;
            case ParticleBoundaryStage::Leavers:
                // TODO: mirror the particles that just left through the face, and move those
                // that land in the patch back to domainParticles.
                break;
            default: break; // reflective does not act at the other moments
        }
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_REFLECTIVE_HPP
