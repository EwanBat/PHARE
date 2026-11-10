#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"


#include <cstddef>

namespace PHARE::core
{
/**
 * @brief Supported types of particle boundary conditions.
 */
enum class ParticleBoundaryConditionType : int {
    None,
    Reflective,
    Open,
    Inflow,
};

/**
 * @brief Interface for applying boundary conditions to particle arrays.
 *
 * Concrete particle boundary conditions are provided by implementating this interface.
 *
 * @tparam ParticleArrayT The type of the particle array.
 * @tparam GridLayoutT The grid layout type .
 *
 */
template<typename ParticleArrayT, typename GridLayoutT>
class IParticleBoundaryCondition
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using This = IParticleBoundaryCondition<ParticleArrayT, GridLayoutT>;
    using particle_array_type = ParticleArrayT;
    using grid_layout_type    = GridLayoutT;

    /** @brief Return the type of the boundary condition. */
    virtual ParticleBoundaryConditionType getType() const = 0;

    virtual ~IParticleBoundaryCondition() = default;

    /** @brief Apply the boundary condition to the given particle array. */
    virtual void apply(ParticleArrayT& particleArray, BoundaryLocation const boundaryLocation,
                      Box<std::uint32_t, dimension> const& localGhostBox,
                      GridLayoutT const& gridLayout, double const time) = 0;
};

} // namespace PHARE::core  

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP
