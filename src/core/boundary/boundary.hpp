#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_HPP
#define PHARE_CORE_BOUNDARY_BOUNDARY_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_resolver.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_resolver.hpp"

#include <concepts>
#include <memory>
#include <unordered_map>
#include <utility>

namespace PHARE::core
{
/**
 * @brief A Boundary is associated to one of the physical boundary (XLower, XUpper) and manages the
 * collection of boundary conditions associated with each physical quantities that requires one.
 * This class is not polymorphic itself, but the different type of boundaries are obtained thanks to
 * the polymorphism of the boundary conditions. Which runtime type to choose for a boundary
 * condition applied to a physical quantity is controled by the @cBoundaryFactory following the
 * desired @c BoundaryType.
 *
 * @tparam GridLayoutT The type for the grid layout.
 * @tparam FieldT The type for the field.
 * @tparam ParticleArrayT The type for the particle array.
 */
template<typename GridLayoutT, IsField FieldT, typename ParticleArrayT = void>
class Boundary
{
public:
    using physical_quantity_type = decltype(GridLayoutT::options.field_options)::Quantity;
    using field_type             = FieldT;
    using This                   = Boundary<GridLayoutT, FieldT, ParticleArrayT>;
    using scalar_quantity_type   = field_type::physical_quantity_type;
    static_assert(std::same_as<scalar_quantity_type, typename physical_quantity_type::Scalar>);
    using vector_quantity_type        = physical_quantity_type::Vector;
    using vector_field_type           = VecField<field_type, physical_quantity_type>;
    using scalar_field_condition_type = IFieldBoundaryCondition<field_type, GridLayoutT>;
    using vector_field_condition_type = IFieldBoundaryCondition<vector_field_type, GridLayoutT>;
    using particle_array_type = ParticleArrayT;
    using particle_condition_type
        = IParticleBoundaryCondition<particle_array_type, field_type, GridLayoutT>;

    Boundary() = delete;
    Boundary(BoundaryType type, BoundaryLocation location)
        : type_{type}
        , location_{location} {};
    // rule of zero: no user-declared destructor, so moves stay implicitly available.

    inline BoundaryType getType() const { return type_; };
    inline BoundaryLocation getLocation() const { return location_; };

    /**
     * @brief Retrieve the registered field boundary condition corresponding to a scalar quantity.
     */
    std::shared_ptr<scalar_field_condition_type>
    getFieldCondition(physical_quantity_type::Scalar quantity)
    {
        auto it = scalar_field_conditions_.find(quantity);
        return (it != scalar_field_conditions_.end()) ? it->second : nullptr;
    }

    /**
     * @brief Retrieve the registered field boundary condition corresponding to a vector quantity.
     */
    std::shared_ptr<vector_field_condition_type>
    getFieldCondition(physical_quantity_type::Vector quantity)
    {
        auto it = vector_field_conditions_.find(quantity);
        return (it != vector_field_conditions_.end()) ? it->second : nullptr;
    }

    /**
     * @brief Retrieve the registered particle boundary condition corresponding to a population.
     */
    std::shared_ptr<particle_condition_type>
    getParticleCondition(std::string const& populationName) const
    {
        auto it = particle_conditions_.find(populationName);
        return (it != particle_conditions_.end()) ? it->second : nullptr;
    }


    /**
     * @brief Register a field boundary condition for a quantity to the boundary.
     *
     * @tparam type The corresponding value of the @c BoundaryType enum corresponding to the desired
     * boundary type.
     * @tparam TensorPhysicalQuantityT Type of the physical quantity, expected to be either @c
     * PhysicalQuantity::Scalar, or @c physical_quantity_type::Vector.
     * @tparam Args Types of the arguments for the FieldConditionT constructor.
     *
     * @param quantity The physical quantity (scalar or vector) to which the field condition should
     * apply.
     * @param args The arguments for the field BC constructor, passed by perfect forwarding to the
     * field BC factory.
     */
    template<FieldBoundaryConditionType type, typename TensorPhysicalQuantityT, typename... Args>
        requires std::same_as<TensorPhysicalQuantityT, scalar_quantity_type>
                 || std::same_as<TensorPhysicalQuantityT, vector_quantity_type>
    void registerFieldCondition(TensorPhysicalQuantityT quantity, Args&&... args)
    {
        if constexpr (std::same_as<TensorPhysicalQuantityT, scalar_quantity_type>)
        {
            scalar_field_conditions_[quantity]
                = std::make_shared<FieldBoundaryCondition<type, field_type, GridLayoutT>>(
                    std::forward<Args>(args)...);
        }
        else
        {
            vector_field_conditions_[quantity]
                = std::make_shared<FieldBoundaryCondition<type, vector_field_type, GridLayoutT>>(
                    std::forward<Args>(args)...);
        }
    }

    /**
     * @brief Register a particle boundary condition for a population to the boundary.
     *
     * @tparam type The corresponding value of the @c ParticleBoundaryConditionType enum
     * corresponding to the desired boundary type.
     * @tparam Args Types of the arguments for the ParticleConditionT constructor.
     *
     * @param populationName The name of the population the particle condition applies to.
     * @param args The arguments for the particle BC constructor, passed by perfect forwarding to
     * the particle BC factory.
     */
    template<ParticleBoundaryConditionType type, typename... Args>
    void registerParticleCondition(std::string const& populationName, Args&&... args)
    {
        using condition_type
            = ParticleBoundaryCondition<type, particle_array_type, field_type, GridLayoutT>;
        particle_conditions_[populationName]
            = std::make_shared<condition_type>(std::forward<Args>(args)...);
    }

private:
    /** The type of the boundary (open, inflow, reflective ...) */
    BoundaryType type_;
    /** The location of the boundary (XLower, XUpper, ...) */
    BoundaryLocation location_;
    /** The list of registered scalar field boundary conditions on the boundary */
    std::unordered_map<scalar_quantity_type, std::shared_ptr<scalar_field_condition_type>>
        scalar_field_conditions_;
    /** The list of registered vector field boundary conditions on the boundary */
    std::unordered_map<vector_quantity_type, std::shared_ptr<vector_field_condition_type>>
        vector_field_conditions_;
    /** The list of registered particle boundary conditions on the boundary */
    std::unordered_map<std::string, std::shared_ptr<particle_condition_type>>
        particle_conditions_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_BOUNDARY_HPP
