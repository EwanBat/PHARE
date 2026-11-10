#ifndef PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP
#define PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP

#include "amr/utilities/boundary/physical_boundary_boxes.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <cstddef>
#include <stdexcept>

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

    explicit ParticleBoundaryPatchStrategy(BoundaryManagerT& boundaryManager)
        : boundaryManager_{boundaryManager}
    {
    }

    /**
     * @brief Fill the boundary particles of @p population on @p patch.
     *
     * An inflow boundary draws its particles in the whole box, the part in front of a same-level
     * neighbour included, so that a particle can enter the patch diagonally.
     */
    template<typename PopulationT>
    void fillBoundaryParticles(SAMRAI::hier::Patch const& patch, PopulationT& population,
                               GridLayoutT const& /*layout*/, double const /*time*/)
    {
        forEachCondition_(patch, population,
                          [](auto& /*condition*/, boundary_box_type const& /*boundaryBox*/) {
                              // nothing yet
                          });
    }

    /**
     * @brief Deposit the boundary particles of @p population on its density and flux.
     *
     * To be called after the push and before the border sums. Particles drawn by an inflow boundary
     * exist in every patch that sees the box, so each patch must deposit only those in front of its
     * own cells.
     */
    template<typename PopulationT>
    void depositBoundaryParticles(SAMRAI::hier::Patch const& patch, PopulationT& population,
                                  GridLayoutT const& /*layout*/, double const /*time*/)
    {
        forEachCondition_(patch, population,
                          [](auto& /*condition*/, boundary_box_type const& /*boundaryBox*/) {
                              // nothing yet
                          });
    }

    /**
     * @brief Treat the particles of @p population that left the patch through a physical boundary
     * during the last push.
     *
     * To be called after the real push and before the particle exchange, which clears
     * patchGhostParticles. An outflow boundary moves them to the boundary particles.
     */
    template<typename PopulationT>
    void collectLeavingParticles(SAMRAI::hier::Patch const& patch, PopulationT& population,
                                 GridLayoutT const& /*layout*/, double const /*time*/)
    {
        forEachCondition_(patch, population,
                          [](auto& /*condition*/, boundary_box_type const& /*boundaryBox*/) {
                              // nothing yet
                          });
    }

private:
    /**
     * @brief Call @p fn with the particle boundary condition of @p population and the box it rules,
     * for each box of particle ghost cells of @p patch that lies outside the physical domain.
     */
    template<typename PopulationT, typename Fn>
    void forEachCondition_(SAMRAI::hier::Patch const& patch, PopulationT& population, Fn&& fn)
    {
        SAMRAI::hier::IntVector const ghostWidth{
            patch.getDim(), static_cast<int>(GridLayoutT::options.particle_ghost_width)};

        forEachPhysicalBoundaryBox<dimension>(
            patch, ghostWidth, boundaryManager_, [&](boundary_box_type const& boundaryBox) {
                auto const condition = boundaryBox.boundary.getParticleCondition(population.name());
                if (!condition)
                    throw std::runtime_error("Particle boundary condition not found.");

                fn(*condition, boundaryBox);
            });
    }

    BoundaryManagerT& boundaryManager_; //!< a reference to the boundary manager
};

} // namespace PHARE::amr

#endif // PHARE_AMR_PARTICLE_BOUNDARY_PATCH_STRATEGY_HPP
