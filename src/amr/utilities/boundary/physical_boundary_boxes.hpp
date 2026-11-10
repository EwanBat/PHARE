#ifndef PHARE_AMR_UTILITIES_BOUNDARY_PHYSICAL_BOUNDARY_BOXES_HPP
#define PHARE_AMR_UTILITIES_BOUNDARY_PHYSICAL_BOUNDARY_BOXES_HPP

#include "amr/utilities/box/amr_box.hpp"

#include "core/boundary/boundary_defs.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"

#include <cassert>
#include <cstddef>
#include <stdexcept>

namespace PHARE::amr
{

/**
 * @brief A box of ghost cells of a patch that lies outside the physical domain, together with the
 * boundary whose conditions apply there.
 *
 * @tparam dim Dimension of the box.
 * @tparam BoundaryT Type of the boundary held by the boundary manager.
 */
template<std::size_t dim, typename BoundaryT>
struct PhysicalBoundaryBox
{
    /** Ghost cells to treat, in AMR indices. Along the boundary it also covers the ghost cells
     * that sit in front of a same-level neighbour, up to the ghost width. */
    core::Box<int, dim> box;

    /** The 1-codimensional boundary that rules this box. It gives the normal direction and the
     * side. Where several boundaries meet (codim > 1) it is the master one. */
    core::BoundaryLocation location;

    BoundaryT& boundary; //!< the boundary found at @c location

    std::size_t codim; //!< 1 for a face, 2 or 3 where several non-periodic faces meet
};


/**
 * @brief Call @p fn once for each box of ghost cells of @p patch that lies outside the physical
 * domain.
 *
 * This is the "which box, on which boundary" loop shared by everything that applies a physical
 * boundary condition on a patch, for fields and for particles. Boxes of codimension 2 and 3 only
 * exist where two or three non-periodic boundaries meet; the boundary manager chooses which one
 * rules them.
 *
 * @tparam dim Dimension of the patch.
 * @param patch The patch. Nothing is called if it does not touch a non-periodic boundary.
 * @param ghostWidth Width of the ghost layer to cover, in cells.
 * @param boundaryManager Manager that knows the boundaries and which one prevails at corners.
 * @param fn Callable taking a @c PhysicalBoundaryBox<dim, BoundaryManagerT::boundary_type>.
 */
template<std::size_t dim, typename BoundaryManagerT, typename Fn>
void forEachPhysicalBoundaryBox(SAMRAI::hier::Patch const& patch,
                                SAMRAI::hier::IntVector const& ghostWidth,
                                BoundaryManagerT const& boundaryManager, Fn&& fn)
{
    using boundary_type = BoundaryManagerT::boundary_type;

    auto const patchGeom = patch.getPatchGeometry();
    assert(patchGeom && "patch has no geometry.");

    // must be retrieved to pass as argument to patchGeom->getBoundaryFillBox later
    SAMRAI::hier::Box const& patchBox = patch.getBox();

    // iterations on potential boundary codimensions in [[1, dim]]
    core::for_N<dim>([&](auto tag) {
        constexpr auto codim = tag.value + 1;

        for (SAMRAI::hier::BoundaryBox const& bBox :
             patchGeom->getCodimensionBoundaries(static_cast<int>(codim)))
        {
            SAMRAI::hier::Box const fillBox
                = patchGeom->getBoundaryFillBox(bBox, patchBox, ghostWidth);

            auto const location
                = static_cast<core::CodimNBoundaryLocation<codim>>(bBox.getLocationIndex());

            // A corner in 2D belongs to two 1-codimensional boundaries, so two boundary conditions
            // compete there. The boundary manager chooses which one prevails. If the box is itself
            // 1-codimensional, the master location is its own location.
            core::BoundaryLocation const masterLocation
                = boundaryManager.getMasterBoundaryLocation(location);

            auto* const boundary = boundaryManager.getBoundary(masterLocation);
            if (!boundary)
                throw std::runtime_error("Boundary not found.");

            fn(PhysicalBoundaryBox<dim, boundary_type>{phare_box_from<dim>(fillBox),
                                                       masterLocation, *boundary, codim});
        }
    });
}

} // namespace PHARE::amr

#endif // PHARE_AMR_UTILITIES_BOUNDARY_PHYSICAL_BOUNDARY_BOXES_HPP
