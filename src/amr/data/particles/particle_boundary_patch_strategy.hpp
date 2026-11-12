#ifndef PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP
#define PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP

#include "amr/utilities/boundary/physical_boundary_boxes.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include "core/def.hpp"

#include <cassert>
#include <cstddef>
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
 * calls it. Particle boundary conditions act at three different moments of a time step, so the
 * messenger calls the three methods below by hand, on patches whose data is already set.
 *
 * "Boundary particles" are the particles kept in the ghost cells beyond a physical boundary. They
 * live in the levelGhostParticles array of the population, which is otherwise unused on the root
 * level.
 *
 * One strategy serves one population. It keeps the name of the population, to find its particle
 * boundary condition, and the patch data ids of its particles, density and flux, to reach them on
 * a patch. So its methods only need the patch.
 *
 * @tparam GridLayoutT The grid layout type.
 * @tparam BoundaryManagerT Manager responsible for providing boundary condition objects.
 */
template<typename GridLayoutT, typename BoundaryManagerT>
class ParticleBoundaryPatchStrategy
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using boundary_type     = BoundaryManagerT::boundary_type;
    using boundary_box_type = PhysicalBoundaryBox<dimension, boundary_type>;

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
     * @brief Fill the boundary particles of the population on @p patch.
     *
     * An inflow boundary draws its particles in the whole box, the part in front of a same-level
     * neighbour included, so that a particle can enter the patch diagonally.
     */
    void fillBoundaryParticles(SAMRAI::hier::Patch const& patch, double const /*time*/)
    {
        forEachCondition_(patch, [](auto& /*condition*/, boundary_box_type const& /*boundaryBox*/) {
            // nothing yet
        });
    }

    /**
     * @brief Deposit the boundary particles of the population on its density and flux.
     *
     * To be called after the push and before the border sums. Particles drawn by an inflow boundary
     * exist in every patch that sees the box, so each patch must deposit only those in front of its
     * own cells.
     */
    void depositBoundaryParticles(SAMRAI::hier::Patch const& patch, double const /*time*/)
    {
        forEachCondition_(patch, [](auto& /*condition*/, boundary_box_type const& /*boundaryBox*/) {
            // nothing yet
        });
    }

    /**
     * @brief Treat the particles of the population that left @p patch through a physical boundary
     * during the last push.
     *
     * To be called after the real push and before the particle exchange, which clears
     * patchGhostParticles. An outflow boundary moves them to the boundary particles.
     */
    void collectLeavingParticles(SAMRAI::hier::Patch const& patch, double const /*time*/)
    {
        forEachCondition_(patch, [](auto& /*condition*/, boundary_box_type const& /*boundaryBox*/) {
            // nothing yet
        });
    }

private:
    void assertIDsSet_() const
    {
        assert(particlesId_ >= 0 && densityId_ >= 0 && fluxId_ >= 0
               && "ParticleBoundaryPatchStrategy: IDs must be registered before use");
    }

    /**
     * @brief Call @p fn with the particle boundary condition of the population and the box it
     * rules, for each box of particle ghost cells of @p patch that lies outside the physical domain.
     */
    template<typename Fn>
    void forEachCondition_(SAMRAI::hier::Patch const& patch, Fn&& fn)
    {
        assertIDsSet_();

        SAMRAI::hier::IntVector const ghostWidth{
            patch.getDim(), static_cast<int>(GridLayoutT::options.particle_ghost_width)};

        forEachPhysicalBoundaryBox<dimension>(
            patch, ghostWidth, boundaryManager_, [&](boundary_box_type const& boundaryBox) {
                auto const condition = boundaryBox.boundary.getParticleCondition(populationName_);
                if (!condition)
                    throw std::runtime_error(
                        "Particle boundary condition not found for population " + populationName_
                        + ".");

                fn(*condition, boundaryBox);
            });
    }

    BoundaryManagerT& boundaryManager_; //!< a reference to the boundary manager
    std::string populationName_;        //!< the population this strategy serves
    int particlesId_ = -1;              //!< patch data id of the population's particles
    int densityId_   = -1;              //!< patch data id of the population's particle density
    int fluxId_      = -1;              //!< patch data id of the population's flux
};

} // namespace PHARE::amr

#endif // PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP
