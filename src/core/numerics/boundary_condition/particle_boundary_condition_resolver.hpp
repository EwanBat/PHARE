#ifndef PHARE_CORE_NUMERICS_PARTICLE_BOUNDARY_CONDITION_RESOLVER
#define PHARE_CORE_NUMERICS_PARTICLE_BOUNDARY_CONDITION_RESOLVER

#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_inflow.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_none.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_open.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_reflective.hpp"

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

template<>
struct ParticleBoundaryConditionSelector<ParticleBoundaryConditionType::Reflective>
{
    template<typename ParticleArrayT, typename GridLayoutT>
    using type = ParticleBoundaryConditionReflective<ParticleArrayT, GridLayoutT>;
};

template<>
struct ParticleBoundaryConditionSelector<ParticleBoundaryConditionType::Open>
{
    template<typename ParticleArrayT, typename GridLayoutT>
    using type = ParticleBoundaryConditionOpen<ParticleArrayT, GridLayoutT>;
};

template<>
struct ParticleBoundaryConditionSelector<ParticleBoundaryConditionType::Inflow>
{
    template<typename ParticleArrayT, typename GridLayoutT>
    using type = ParticleBoundaryConditionInflow<ParticleArrayT, GridLayoutT>;
};

template<ParticleBoundaryConditionType type, typename ParticleArrayT, typename GridLayoutT>
using ParticleBoundaryCondition
    = ParticleBoundaryConditionSelector<type>::template type<ParticleArrayT, GridLayoutT>;

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_PARTICLE_BOUNDARY_CONDITION_RESOLVER
