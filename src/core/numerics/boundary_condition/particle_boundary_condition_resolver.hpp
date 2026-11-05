#ifndef PHARE_CORE_NUMERICS_PARTICLE_BOUNDARY_CONDITION_RESOLVER
#define PHARE_CORE_NUMERICS_PARTICLE_BOUNDARY_CONDITION_RESOLVER

#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_none.hpp"

namespace PHARE::core
{

template<ParticleBoundaryConditionType type>
struct ParticleBoundaryConditionSelector;

template<>
struct ParticleBoundaryConditionSelector<ParticleBoundaryConditionType::None>
{
    template<typename ParticleArrayT, typename GridLayoutT>
    using type = ParticleBoundaryConditionNone<ParticleArrayT, GridLayoutT>;
};

template<ParticleBoundaryConditionType type, typename ParticleArrayT, typename GridLayoutT>
using ParticleBoundaryCondition
    = ParticleBoundaryConditionSelector<type>::template type<ParticleArrayT, GridLayoutT>;

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_PARTICLE_BOUNDARY_CONDITION_RESOLVER
