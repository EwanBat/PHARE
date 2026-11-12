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
 * Physics not implemented yet: does nothing.
 */
template<typename ParticleArrayT, typename GridLayoutT>
class ParticleBoundaryConditionInflow
    : public IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>
{
public:
    using Super = IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>;

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

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Inflow;
    }

    void apply(ParticleArrayT& /*particleArray*/, BoundaryLocation const /*boundaryLocation*/,
               Box<std::uint32_t, Super::dimension> const& /*localGhostBox*/,
               GridLayoutT const& /*gridLayout*/, [[maybe_unused]] double const /*time*/) override
    {
    }

private:
    double density_;
    std::array<double, 3> velocity_;
    std::array<double, 3> thermalVelocity_;
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_INFLOW_HPP
