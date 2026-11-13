#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/utilities/box/box.hpp"

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
 * @brief The moments of a time step where a particle boundary condition may act.
 *
 * Each condition acts at some of these moments only, and says which ones in actsAt(). At each
 * moment, only the conditions that act then are called. Order in a time step:
 *
 *     Init      when the root level is created, or regridded
 *     Moments   then, once the domain particles are deposited
 *     ----      push (part 1 or part 2 of the step)
 *     Refill    after the real push only (part 2)
 *     Moments   after each push, before the border sums
 *     ----      border sums
 *     Leavers   after the real push only, before the particle exchange
 *     ----      particle exchange
 *
 * What each condition does:
 *
 *                 Inflow              Open                    Reflective
 *     Init        draw the stock      -                       -
 *     Refill      redraw the stock    -                       -
 *     Moments     deposit the stock   deposit the stock       make rho and F symmetric
 *     Leavers     -                   move leavers to stock   mirror leavers, give them back
 *
 * The "stock" is the set of particles kept in the ghost cells beyond the face.
 */
enum class ParticleBoundaryStage : int {
    Init,
    Refill,
    Moments,
    Leavers,
};


/**
 * @brief What a particle boundary condition may read or change: the particle arrays and the
 * moments of one population on one patch.
 *
 * It only holds references. It is built just before a call and thrown away after.
 */
template<typename ParticleArrayT, typename FieldT, typename VecFieldT>
struct ParticleBoundaryData
{
    ParticleArrayT& domainParticles;     //!< particles inside the patch
    ParticleArrayT& patchGhostParticles; //!< copies of the particles that just left the patch
    ParticleArrayT& levelGhostParticles; //!< on the root level: the stock beyond the faces
    FieldT& density;                     //!< particle density of the population
    VecFieldT& flux;                     //!< flux of the population
};


/**
 * @brief Interface of the particle boundary conditions.
 *
 * A concrete condition (inflow, open, ...) says at which moments of the step it acts (actsAt)
 * and what it does then (apply). apply is only called at the moments where actsAt is true.
 *
 * @tparam ParticleArrayT The type of the particle array.
 * @tparam FieldT The type of a scalar field (the density).
 * @tparam GridLayoutT The grid layout type.
 */
template<typename ParticleArrayT, typename FieldT, typename GridLayoutT>
class IParticleBoundaryCondition
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using particle_array_type = ParticleArrayT;
    using field_type          = FieldT;
    using grid_layout_type    = GridLayoutT;
    using physical_quantity_type = decltype(GridLayoutT::options.field_options)::Quantity;
    using vector_field_type      = VecField<FieldT, physical_quantity_type>;
    using data_type = ParticleBoundaryData<ParticleArrayT, FieldT, vector_field_type>;

    virtual ~IParticleBoundaryCondition() = default;

    /** @brief Return the type of the boundary condition. */
    virtual ParticleBoundaryConditionType getType() const = 0;

    /** @brief True if the condition has something to do at this moment of the step. */
    virtual bool actsAt(ParticleBoundaryStage const stage) const = 0;

    /**
     * @brief Do the work of this moment of the step, in one box of ghost cells beyond a face.
     *
     * @param stage the moment of the step. actsAt(stage) is true.
     * @param data the particles and moments of the population on the patch.
     * @param location the face that rules the box. It gives the normal direction and the side.
     * @param ghostBox the box of ghost cells, in AMR indices (like a particle's iCell).
     * @param gridLayout the layout of the patch.
     * @param time the time of the step.
     */
    virtual void apply(ParticleBoundaryStage const stage, data_type& data,
                       BoundaryLocation const location, Box<int, dimension> const& ghostBox,
                       GridLayoutT const& gridLayout, double const time)
        = 0;
};

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP
