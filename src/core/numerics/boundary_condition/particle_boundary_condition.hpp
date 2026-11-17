#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/def.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

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
 *     Init              when the root level is created, or regridded
 *     Moments           then, once the domain particles are deposited
 *     ----              push of copies (domain_only, part 1 of the step)
 *     PredictedMoments  after it, before the border sums
 *     ----              real push (all, part 2 of the step)
 *     Refill            after it
 *     Moments           then, before the border sums
 *     ----              border sums
 *     Leavers           before the particle exchange (part 2 only)
 *     ----              particle exchange
 *
 * PredictedMoments and Moments both complete the moments at the time of the push. They differ in
 * what the particle arrays hold: after the push of copies, the arrays have not moved (they are
 * still at the start of the step); after the real push, they have.
 *
 * What each condition does:
 *
 *                     Inflow                  Open                    Reflective
 *     Init            draw the stock          -                       -
 *     PredictedMoments draw a temp stock at   (to decide)             make rho and F symmetric
 *                     the new time, deposit
 *     Refill          redraw the stock        -                       -
 *     Moments         deposit the stock       deposit the stock       make rho and F symmetric
 *     Leavers         erase the strays        move leavers to stock   mirror leavers, give them back
 *                                             strays → patch ghosts
 *
 * The "stock" is the set of particles kept in the ghost cells beyond the face. The "strays" are
 * stock particles that the real push took into a neighbour's cells (see strayBoxes): a copy for
 * inflow (the neighbour drew the same particle), the only one for open (a leaver).
 */
enum class ParticleBoundaryStage : int {
    Init,
    PredictedMoments,
    Refill,
    Moments,
    Leavers,
};


/**
 * @brief The part of @p ghostBox this patch deposits: the cells in front of its own cells.
 *
 * A box of ghost cells beyond a face also covers, along the face, the cells in front of a
 * same-level neighbour (up to the ghost width). The neighbour holds the same stock there (draws
 * are keyed by cell), so if both deposited it, the border sums would count it twice. Each patch
 * therefore deposits only the cells in front of its own cells, as it does for its domain
 * particles, and the border sums give the rest.
 *
 *              neighbour  │   this patch    │  neighbour
 *     domain   ...........│.................│...........
 *     ──────────────────── face ──────────────────────────
 *     ghostBox        [ g │ o o o o o o o o │ g ]
 *     result              [ o o o o o o o o ]
 *
 * In a direction where the whole box is beyond the patch (the normal of a face, or both normals
 * of a corner box), nothing is cut. In the other directions, the box is cut to the patch.
 *
 * Same-level patches that overlap (non-owned cells) are not handled: on the root level they do
 * not happen in practice.
 */
template<std::size_t dim>
NO_DISCARD Box<int, dim> frontOfOwnCells(Box<int, dim> ghostBox, Box<int, dim> const& domainBox)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        bool const overlaps = ghostBox.upper[d] >= domainBox.lower[d]
                              and ghostBox.lower[d] <= domainBox.upper[d];
        if (overlaps)
        {
            ghostBox.lower[d] = std::max(ghostBox.lower[d], domainBox.lower[d]);
            ghostBox.upper[d] = std::min(ghostBox.upper[d], domainBox.upper[d]);
        }
    }
    return ghostBox;
}


/**
 * @brief The cells of @p box that are not in @p hole, as boxes that do not overlap.
 *
 * Direction by direction, the slabs of @p box below and above @p hole are cut off, and what is
 * left is narrowed to @p hole in that direction before the next one. In 2D:
 *
 *           ┌───┬───────┬───┐
 *           │   │   4   │   │    d = 0 (x): 1 and 2, left and right of the hole,
 *           │   ├───────┤   │               over the whole height of the box
 *         y │ 1 │ hole  │ 2 │
 *         ↑ │   ├───────┤   │    d = 1 (y): 3 and 4, below and above the hole,
 *         │ │   │   3   │   │               only between its x ends (1 and 2 are taken)
 *         │ └───┴───────┴───┘
 *         └──→ x
 *
 * In 3D the z slabs come last, only over the x and y extent of the hole. A slab of zero width is
 * not made, so a hole touching a side of the box gives fewer boxes. The hole is first cut to the
 * box: a hole outside the box gives {box}, a hole covering it gives {}.
 */
template<std::size_t dim>
NO_DISCARD std::vector<Box<int, dim>> boxMinus(Box<int, dim> box, Box<int, dim> const& hole)
{
    auto const inter = box * hole;
    if (!inter)
        return {box};

    std::vector<Box<int, dim>> parts;
    for (std::size_t d = 0; d < dim; ++d)
    {
        if (box.lower[d] < inter->lower[d])
        {
            auto below     = box;
            below.upper[d] = inter->lower[d] - 1;
            parts.push_back(below);
        }
        if (box.upper[d] > inter->upper[d])
        {
            auto above     = box;
            above.lower[d] = inter->upper[d] + 1;
            parts.push_back(above);
        }
        box.lower[d] = inter->lower[d];
        box.upper[d] = inter->upper[d];
    }
    return parts;
}


/**
 * @brief The cells where a stock particle of @p ghostBox may land, after one real push, inside the
 * simulation domain but in front of a same-level neighbour: the "strays".
 *
 * The ghost box of a face also covers the cells in front of a neighbour (see frontOfOwnCells). A
 * stock particle there can move inward, across the face, into a cell of that neighbour. The
 * updater keeps it in levelGhostParticles (it is in the ghost layer and not owned), where nobody
 * would ever deposit or send it. Each condition decides what to do with it at the Leavers moment.
 *
 *              x=-1 │ x=0
 *     neighbour [o] │ [s]        o = ghostBox
 *     ──────────────┼────        s = result
 *     patch     [o] │ patch
 *     ──────────────┼────
 *     neighbour [o] │ [s]
 *
 * The box is shifted inward by its own depth (a particle moves less than the ghost width per
 * step), then the patch is removed. No special case per codimension:
 *  - a corner box (beyond the patch in every direction) lands inside the patch: no stray, and
 *    boxMinus returns an empty vector at once (no allocation);
 *  - an edge box in 3D (beyond in two directions) can still reach a neighbour along the third,
 *    so it may have strays: skipping every box of codim > 1 would be wrong in 3D.
 */
template<std::size_t dim>
NO_DISCARD std::vector<Box<int, dim>> strayBoxes(Box<int, dim> const& ghostBox,
                                                 Box<int, dim> const& domainBox)
{
    std::array<int, dim> inward{}; // zero in the directions where the box faces the patch
    for (std::size_t d = 0; d < dim; ++d)
    {
        int const depth = ghostBox.upper[d] - ghostBox.lower[d] + 1;
        if (ghostBox.upper[d] < domainBox.lower[d]) // beyond the lower side
            inward[d] = depth;
        else if (ghostBox.lower[d] > domainBox.upper[d]) // beyond the upper side
            inward[d] = -depth;
    }
    return boxMinus(shift(ghostBox, Point<int, dim>{inward}), domainBox);
}


/**
 * @brief What a particle boundary condition needs to know about one ion population.
 *
 * Taken from the population's particle initializer, so that the particles drawn at a boundary
 * look like those drawn in the domain (same number per cell, so the same weight for the same
 * density).
 */
struct ParticleBoundaryInfo
{
    std::string name;                      //!< name of the population
    std::uint32_t nbrParticlesPerCell = 0; //!< particles drawn per cell
    std::optional<std::size_t> seed{};     //!< random seed; none means a random one
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
