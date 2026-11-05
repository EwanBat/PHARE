#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP

#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

namespace PHARE::core
{
// TODO: implement (forward declaration only, so the resolver compiles)
template<typename ParticleArrayT, typename GridLayoutT>
class ParticleBoundaryConditionNone : public IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>
{
public:
    using Super = IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>;
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP
