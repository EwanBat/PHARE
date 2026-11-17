#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

namespace PHARE::core
{
/**
 * @brief Particle boundary condition for an open (outflow) boundary.
 *
 * Particles that leave the domain through the boundary are kept in the ghost cells beyond it, and
 * deposited there, until they leave the ghost layer.
 *
 * Physics not implemented yet: the moments are chosen, the work is not written.
 */
template<typename ParticleArrayT, typename FieldT, typename GridLayoutT>
class ParticleBoundaryConditionOpen
    : public IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>
{
public:
    using Super     = IParticleBoundaryCondition<ParticleArrayT, FieldT, GridLayoutT>;
    using data_type = Super::data_type;

    /** @brief Open acts at PredictedMoments, Moments and Leavers. */
    bool actsAt(ParticleBoundaryStage const stage) const override
    {
        return stage == ParticleBoundaryStage::PredictedMoments
               or stage == ParticleBoundaryStage::Moments or stage == ParticleBoundaryStage::Leavers;
    }

    ParticleBoundaryConditionType getType() const override
    {
        return ParticleBoundaryConditionType::Open;
    }

    void apply(ParticleBoundaryStage const stage, data_type& data,
               BoundaryLocation const /*location*/, Box<int, Super::dimension> const& ghostBox,
               GridLayoutT const& gridLayout, double const /*time*/) override
    {
        switch (stage)
        {
            case ParticleBoundaryStage::PredictedMoments:
                // TODO: what to deposit after the push of copies is not decided yet: the real
                // stock has not moved, and the leavers of the copies were dropped. See the note
                // "Issues for later - Hybrid BC", "Open BC: what to deposit during domain_only".
                break;
            case ParticleBoundaryStage::Moments:
                // TODO: deposit the stock (the particles that left in earlier steps).
                break;
            case ParticleBoundaryStage::Leavers:
            {
                // strays: stock particles the real push took into a neighbour's cells. This
                // patch holds the only copy (a leaver is not drawn), so it is handed to the
                // neighbour: the particle exchange sends the patch ghosts in a neighbour's cells
                // to that neighbour's domain (ParticlesData::copy_from_ghost).
                auto& stock = data.levelGhostParticles;
                for (auto const& strays : strayBoxes(ghostBox, gridLayout.AMRBox()))
                {
                    stock.export_particles(strays, data.patchGhostParticles);
                    stock.erase_in(strays);
                }
                // TODO: move the particles that just left through the face from
                // patchGhostParticles to levelGhostParticles (the stock). The strays above are
                // in the neighbour's cells, not beyond the face, so they are not taken back.
                break;
            }
            case ParticleBoundaryStage::Init:
            case ParticleBoundaryStage::Refill: break; // open does not act then
        }
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_PARTICLE_BOUNDARY_CONDITION_OPEN_HPP
