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
 * Physics not implemented yet: does nothing.
 */
template<typename ParticleArrayT, typename GridLayoutT>
class ParticleBoundaryConditionReflective
    : public IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>
{
public:
    using Super = IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>;

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Reflective;
    }

    void apply(ParticleArrayT& /*particleArray*/, BoundaryLocation const /*boundaryLocation*/,
               Box<std::uint32_t, Super::dimension> const& /*localGhostBox*/,
               GridLayoutT const& /*gridLayout*/, [[maybe_unused]] double const /*time*/) override
    {
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_REFLECTIVE_HPP
