#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

namespace PHARE::core
{
/**
 * @brief No-op particle boundary condition: leaves particles untouched.
 */
template<typename ParticleArrayT, typename FieldT, typename GridLayoutT>
class ParticleBoundaryConditionNone
    : public IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>
{
public:
    using Super     = IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>;
    using data_type = Super::data_type;

    /** @brief Never acts. */
    bool actsAt(ParticleBoundaryStage const /*stage*/) const override { return false; }

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::None;
    }

    void apply(ParticleBoundaryStage const /*stage*/, data_type& /*data*/,
               BoundaryLocation const /*location*/,
               Box<int, Super::dimension> const& /*ghostBox*/, GridLayoutT const& /*gridLayout*/,
               double const /*time*/) override
    {
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_NONE_HPP
