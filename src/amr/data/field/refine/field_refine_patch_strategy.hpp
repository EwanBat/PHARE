#ifndef PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP

#include "amr/data/field/field_data_traits.hpp"
#include "amr/data/tensorfield/tensor_field_data_traits.hpp"
#include "amr/utilities/boundary/physical_boundary_boxes.hpp"

#include "core/boundary/boundary_defs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

namespace PHARE::amr
{

/**
 * @brief Strategy for filling physical boundary conditions and customizing patch refinement.
 *
 * This class implements the SAMRAI::xfer::RefinePatchStrategy interface to
 * specify how physical boundary conditions must be enforced for patches that touch
 * the domain boundaries. Refinement customization via preprocessRefine and postProcessRefine is
 * deferred to child classes.
 *
 * Each quantity a refiner fills is expected to have its own instance of this class (combined by a
 * CompositeRefinePatchStrategy when a refiner fills several), but all these instances will point
 * to the same common boundary manager.
 *
 * @tparam ScalarOrTensorFieldDataT The data type for fields or tensor fields.
 * @tparam BoundaryManagerT Manager responsible for providing boundary condition objects.
 */
template<typename ScalarOrTensorFieldDataT, typename BoundaryManagerT>
    requires(IsFieldData<ScalarOrTensorFieldDataT> || IsTensorFieldData<ScalarOrTensorFieldDataT>)
class FieldRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    static constexpr bool is_scalar        = IsFieldData<ScalarOrTensorFieldDataT>;
    static constexpr bool is_tensor        = !is_scalar;
    static constexpr std::size_t dimension = ScalarOrTensorFieldDataT::dimension;

    using field_geometry_type    = FieldGeometrySelector<ScalarOrTensorFieldDataT, is_scalar>::type;
    using gridlayout_type        = ScalarOrTensorFieldDataT::gridlayout_type;
    using grid_type              = ScalarOrTensorFieldDataT::grid_type;
    using field_type             = grid_type::field_type;
    using physical_quantity_type = BoundaryManagerT::physical_quantity_type;
    using vectorfield_type       = core::VecField<field_type, physical_quantity_type>;
    using scalar_or_tensor_field_type
        = ScalarOrTensorFieldSelector<ScalarOrTensorFieldDataT, is_scalar>::type;
    using scalar_quantity_type = physical_quantity_type::Scalar;
    using vector_quantity_type = physical_quantity_type::Vector;

    using patch_geometry_type = SAMRAI::hier::PatchGeometry;

    using boundary_type = BoundaryManagerT::boundary_type;
    using boundary_condition_type
        = core::IFieldBoundaryCondition<scalar_or_tensor_field_type, gridlayout_type>;

    /**
     * @brief Constructor.
     * @param boundary_manager Manager handling boundary conditions.
     */
    FieldRefinePatchStrategy(BoundaryManagerT& boundaryManager)
        : boundaryManager_{boundaryManager}
        , data_id_{-1}
    {
    }

    /**
     * @brief Check that the patch data identifier is registered.
     */
    void assertIDsSet() const
    {
        assert(data_id_ >= 0 && "FieldRefinePatchStrategy: IDs must be registered before use");
    }

    /**
     * @brief Register the SAMRAI patch data identifier.
     * @param field_id Integer ID from the SAMRAI variable database.
     */
    void registerIDs(int const field_id) { data_id_ = field_id; }

    /**
     * @brief Apply physical boundary conditions via SAMRAI callback.
     *
     * Iterate over patch boundaries that touch a physical domain boundary and apply the appropriate
     * PHARE boundary condition to ghost regions.
     *
     * @param patch The fine patch being refined.
     * @param fill_time Simulation time for BC application.
     * @param ghost_width_to_fill Width of ghost cell layer to be filled.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                       SAMRAI::hier::IntVector const& ghost_width_to_fill) override
    {
        gridlayout_type const& gridLayout = ScalarOrTensorFieldDataT::getLayout(patch, data_id_);

        assert(ghost_width_to_fill <= SAMRAI::hier::IntVector(
                   static_cast<SAMRAI::tbox::Dimension>(static_cast<int>(dimension)),
                   static_cast<int>(gridLayout.options.field_ghost_width)));

        auto scalarOrTensorField = [&]() {
            if constexpr (is_scalar)
            {
                return *(&(ScalarOrTensorFieldDataT::getField(patch, data_id_)));
            }
            else
            {
                return ScalarOrTensorFieldDataT::getTensorField(patch, data_id_);
            };
        }();

        forEachPhysicalBoundaryBox<dimension>(
            patch, ghost_width_to_fill, boundaryManager_, [&](auto const& boundaryBox) {
                // get the boundary condition for the current physical quantity.
                std::shared_ptr<boundary_condition_type> bc
                    = boundaryBox.boundary.getFieldCondition(
                        scalarOrTensorField.physicalQuantity());
                if (!bc)
                    throw std::runtime_error("Field boundary condition not found.");

                // apply the retained boundary condition, as if the current box was belonging to
                // the 1-codimensional master boundary; this essentially defines which Cartesian
                // direction is considered to be the normal one.
                bc->apply(scalarOrTensorField, boundaryBox.location,
                          gridLayout.AMRToLocal(boundaryBox.box), gridLayout, fill_time);
            });
    }


    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector{dim, 1};
    }


    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box,
                          SAMRAI::hier::IntVector const& ratio) override
    {
    }


    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
    }


protected:
    BoundaryManagerT& boundaryManager_; //!< a reference to the boundary manager
    int data_id_; //!< the id of the resource to which this refine patch strategy is attached
};

} // namespace PHARE::amr

#endif // PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
