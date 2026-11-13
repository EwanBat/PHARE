#ifndef PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP
#define PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP

#include "amr/data/field/field_data.hpp"
#include "amr/data/particles/particles_data.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/utilities/boundary/physical_boundary_boxes.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include "core/def.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace PHARE::amr
{

/**
 * @brief Applies the particle boundary conditions of a population on one patch.
 *
 * It plays for particles the role FieldRefinePatchStrategy plays for fields: it finds the boxes of
 * ghost cells of the patch that lie outside the physical domain, finds the boundary that rules each
 * of them, and hands the work to the particle boundary condition registered there for the
 * population.
 *
 * Unlike FieldRefinePatchStrategy it is not a SAMRAI::xfer::RefinePatchStrategy and no schedule
 * calls it. Particle boundary conditions act at different moments of a time step, so the
 * messenger calls apply() by hand, on patches whose data is already set, and says which moment
 * it is (core::ParticleBoundaryStage).
 *
 * One strategy serves one population. It keeps the name of the population, to find its particle
 * boundary condition, and the patch data ids of its particles, density and flux, to reach them on
 * a patch. So apply() only needs the patch.
 *
 * @tparam GridLayoutT The grid layout type.
 * @tparam GridT The grid type, which holds the data of a field.
 * @tparam BoundaryManagerT Manager responsible for providing boundary condition objects.
 */
template<typename GridLayoutT, typename GridT, typename BoundaryManagerT>
class ParticleBoundaryPatchStrategy
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using boundary_type           = BoundaryManagerT::boundary_type;
    using boundary_box_type       = PhysicalBoundaryBox<dimension, boundary_type>;
    using boundary_condition_type = boundary_type::particle_condition_type;
    using data_type               = boundary_condition_type::data_type;

    using particle_array_type = boundary_type::particle_array_type;
    using particles_data_type = ParticlesData<particle_array_type>;
    using field_data_type     = FieldData<GridLayoutT, GridT>;
    using vecfield_data_type
        = TensorFieldData<1, GridLayoutT, GridT, typename boundary_type::physical_quantity_type>;

    ParticleBoundaryPatchStrategy(BoundaryManagerT& boundaryManager, std::string populationName)
        : boundaryManager_{boundaryManager}
        , populationName_{std::move(populationName)}
    {
    }

    /**
     * @brief Register the SAMRAI patch data ids of the population.
     * @param particlesId id of the population's particles (ParticlesData).
     * @param densityId id of the population's particle density.
     * @param fluxId id of the population's flux.
     */
    void registerIDs(int const particlesId, int const densityId, int const fluxId)
    {
        particlesId_ = particlesId;
        densityId_   = densityId;
        fluxId_      = fluxId;
    }

    NO_DISCARD std::string const& populationName() const { return populationName_; }

    /**
     * @brief Apply, on @p patch, the boundary conditions that act at this moment of the step.
     *
     * For each box of ghost cells beyond a face, the condition of that face is asked if it acts
     * at @p stage. If not, the box is skipped. So, for example, at the Refill moment only the
     * inflow faces do something.
     *
     * @param patch the patch. Nothing happens if it does not touch a physical face.
     * @param stage the moment of the step.
     * @param time the time of the step.
     */
    void apply(SAMRAI::hier::Patch const& patch, core::ParticleBoundaryStage const stage,
               double const time)
    {
        assertIDsSet_();

        forEachPhysicalBoundaryBox<dimension>(
            patch, ghostWidth_, boundaryManager_, [&](boundary_box_type const& boundaryBox) {
                auto const condition = getCondition_(boundaryBox);
                if (!condition->actsAt(stage))
                    return;

                // the data is reached only here, so a face with nothing to do costs nothing
                auto& particles = getParticlesData_(patch);
                auto flux       = vecfield_data_type::getTensorField(patch, fluxId_);
                data_type data{particles.domainParticles, particles.patchGhostParticles,
                               particles.levelGhostParticles,
                               *field_data_type::getField(patch, densityId_), flux};

                condition->apply(stage, data, boundaryBox.location, boundaryBox.box,
                                 field_data_type::getLayout(patch, densityId_), time);
            });
    }

private:
    void assertIDsSet_() const
    {
        assert(particlesId_ >= 0 && densityId_ >= 0 && fluxId_ >= 0
               && "ParticleBoundaryPatchStrategy: IDs must be registered before use");
    }

    /** @brief The particle boundary condition of the population on the boundary of @p boundaryBox. */
    std::shared_ptr<boundary_condition_type> getCondition_(boundary_box_type const& boundaryBox) const
    {
        auto condition = boundaryBox.boundary.getParticleCondition(populationName_);
        if (!condition)
            throw std::runtime_error("Particle boundary condition not found for population "
                                     + populationName_ + ".");
        return condition;
    }

    /** @brief The particle arrays of the population on @p patch. */
    particles_data_type& getParticlesData_(SAMRAI::hier::Patch const& patch) const
    {
        auto particles
            = std::dynamic_pointer_cast<particles_data_type>(patch.getPatchData(particlesId_));
        if (!particles)
            throw std::runtime_error("cannot cast to ParticlesData");
        return *particles;
    }

    BoundaryManagerT& boundaryManager_; //!< a reference to the boundary manager
    std::string populationName_;        //!< the population this strategy serves

    //! width of the particle ghost layer, the layer the boundary conditions act on
    SAMRAI::hier::IntVector const ghostWidth_{
        SAMRAI::tbox::Dimension{static_cast<unsigned short>(dimension)},
        static_cast<int>(GridLayoutT::options.particle_ghost_width)};

    int particlesId_ = -1; //!< patch data id of the population's particles
    int densityId_   = -1; //!< patch data id of the population's particle density
    int fluxId_      = -1; //!< patch data id of the population's flux
};

} // namespace PHARE::amr

#endif // PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP
