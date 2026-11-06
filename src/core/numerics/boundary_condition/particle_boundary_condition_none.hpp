#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP

#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

namespace PHARE::core
{
/**
 * @brief No-op particle boundary condition: leaves particles untouched.
 */
template<typename ParticleArrayT, typename GridLayoutT>
class ParticleBoundaryConditionNone : public IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>
{
public:
    using Super = IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>;

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::None;
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP
