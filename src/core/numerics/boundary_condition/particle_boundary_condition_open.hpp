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
 * Physics not implemented yet: does nothing.
 */
template<typename ParticleArrayT, typename GridLayoutT>
class ParticleBoundaryConditionOpen : public IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>
{
public:
    using Super = IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>;

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Open;
    }

    void apply(ParticleArrayT& /*particleArray*/, BoundaryLocation const /*boundaryLocation*/,
               Box<std::uint32_t, Super::dimension> const& /*localGhostBox*/,
               GridLayoutT const& /*gridLayout*/, [[maybe_unused]] double const /*time*/) override
    {
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP
