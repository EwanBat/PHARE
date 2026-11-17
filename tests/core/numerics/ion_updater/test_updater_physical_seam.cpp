// A stock particle near the seam of two patches that both touch a physical face.
//
// Root level, 2D, order 1 (particle ghost width 1). Physical face x-lower at x = 0.
// Two patches, one above the other, share a seam at y = 4|5:
//
//             physical │
//              ghost   │  domain
//     Q  y=5  [-1, 5]  │ [0, 5] ...       Q: cells x in [0, 9], y in [5, 9]
//     ── seam ─────────┼─────────────
//     P  y=4  [-1, 4]  │ [0, 4] ...       P: cells x in [0, 9], y in [0, 4]
//
// The ghost box of a face also covers the cells in front of the same-level neighbour. So the
// stock cell [-1, 5] (in front of Q) is also in P's ghost box, and both patches hold the same
// stock particle there (draws are keyed by cell). Same for [-1, 4].
//
// After the real push the updater keeps, in levelGhostParticles, a stock particle that landed in a
// neighbour's cells (a "stray", see strayBoxes). At the Leavers moment the inflow rule erases its
// strays (the neighbour has the same particle) and the open rule hands them to the neighbour
// through patchGhostParticles (it holds the only copy).
//
// The test pushes with the real IonUpdater (E = B = 0, straight lines), calls the rules at the
// Leavers moment as SolverPPC does, and stands in for the SAMRAI particle exchange: what a patch
// gets from its neighbour is what the neighbour put in patchGhostParticles in the patch's cells.
// A particle is tracked by its unique weight.

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_inflow.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_open.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/utilities/box/box.hpp"

#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace PHARE::core;

namespace
{
constexpr std::size_t dim          = 2;
constexpr std::size_t interp_order = 1;
constexpr PHARE::SimOpts opts{dim, interp_order};

using PHARETypes    = PHARE::core::PHARE_Types<opts>;
using Ions          = PHARETypes::Hybrid::Ions_t;
using Electromag    = PHARETypes::Hybrid::Electromag_t;
using GridLayout    = PHARETypes::Hybrid::GridLayout_t;
using ParticleArray = PHARETypes::Hybrid::ParticleArray_t;
using Grid          = PHARETypes::Hybrid::Grid_t;
using IonUpdater_t  = IonUpdater<Ions, Electromag, GridLayout>;
using Boxing_t      = UpdaterSelectionBoxing<IonUpdater_t, GridLayout>;
using Box_t         = Box<int, dim>;
using Particle_t    = ParticleArray::value_type;
using Field_t       = Grid::field_type;
using Inflow_t      = ParticleBoundaryConditionInflow<ParticleArray, Field_t, GridLayout>;
using Open_t        = ParticleBoundaryConditionOpen<ParticleArray, Field_t, GridLayout>;
using Stage         = ParticleBoundaryStage;

constexpr double dx = 0.1;
constexpr double dt = 0.01;
// half a cell per step along each axis
constexpr double v = 0.5 * dx / dt;

// the particle we follow: no other particle has this weight
constexpr double tracked = 0.123;

Box_t const boxP{{0, 0}, {9, 4}};
Box_t const boxQ{{0, 5}, {9, 9}};


PHARE::initializer::PHAREDict makeDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    dict["ions"]["nbrPopulations"]                               = std::size_t{1};
    dict["ions"]["pop0"]["name"]                                 = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                                 = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["charge"]       = 1.;
    return dict;
}
auto const dict = makeDict();


// one patch: its layout, fields, one population, and the boxes the solver would build
struct Patch
{
    static constexpr int mapMargin = 2; // cell map a bit larger than the ghost box

    GridLayout layout;
    UsableVecField<dim> E, B;
    Electromag em{"EM"};

    Grid density, chargeDensity, ionChargeDensity, ionMassDensity;
    UsableVecField<dim> flux, Vi, chargeFlux;
    UsableTensorField<dim> M, popM;
    ParticleArray domain, patchGhost, levelGhost, levelGhostOld, levelGhostNew;
    ParticlesPack<ParticleArray> pack;
    Ions ions{dict["ions"]};

    Boxing_t boxing;

    Patch(Box_t const& box, Box_t const& neighbour)
        : layout{{dx, dx},
                 {static_cast<std::uint32_t>(box.shape()[0]),
                  static_cast<std::uint32_t>(box.shape()[1])},
                 {box.lower[0] * dx, box.lower[1] * dx},
                 box}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
        , B{"EM_B", layout, HybridQuantity::Vector::B}
        , density{"protons_particleDensity", HybridQuantity::Scalar::rho,
                  layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , chargeDensity{"protons_chargeDensity", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , ionChargeDensity{"chargeDensity", HybridQuantity::Scalar::rho,
                           layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , flux{"protons_flux", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , chargeFlux{"chargeFlux", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , popM{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , domain{grow(box, mapMargin)}
        , patchGhost{grow(box, mapMargin)}
        , levelGhost{grow(box, mapMargin)}
        , levelGhostOld{grow(box, mapMargin)}
        , levelGhostNew{grow(box, mapMargin)}
        , pack{"protons", &domain, &patchGhost, &levelGhost, &levelGhostOld, &levelGhostNew}
        // what makeNonLevelGhostBoxFor gives: the patch, plus its particle ghost box cut to the
        // same-level neighbour. Nothing beyond the physical face.
        , boxing{layout,
                 {box, *(grow(box, GridLayout::options.particle_ghost_width) * neighbour)},
                 {}}
    {
        E.set_on(em.E);
        B.set_on(em.B);

        auto const& [V, m, cd, md, cf] = ions.getCompileTimeResourcesViewList();
        Vi.set_on(V);
        chargeFlux.set_on(cf);
        M.set_on(m);
        cd.setBuffer(&ionChargeDensity);
        md.setBuffer(&ionMassDensity);

        auto& pops                           = ions.getRunTimeResourcesViewList();
        auto const& [F, popm, d, c, particles] = pops[0].getCompileTimeResourcesViewList();
        d.setBuffer(&density);
        c.setBuffer(&chargeDensity);
        popM.set_on(popm);
        flux.set_on(F);
        particles.setBuffer(&pack);
    }

    void push()
    {
        IonUpdater_t updater{dict["simulation"]["algo"]["ion_updater"]};
        updater.updatePopulations(ions, em, boxing, dt, UpdaterMode::all);
    }

    // the cells of this patch's x-lower face ghost box: what Refill empties (erase_in)
    Box_t faceGhostBox() const
    {
        return Box_t{{-1, layout.AMRBox().lower[1] - 1}, {-1, layout.AMRBox().upper[1] + 1}};
    }

    // what SolverPPC does after the real push, before the exchange
    void leavers(auto& condition)
    {
        typename Inflow_t::data_type data{domain, patchGhost, levelGhost, *(&density), flux};
        condition.apply(Stage::Leavers, data, BoundaryLocation::XLower, faceGhostBox(), layout,
                        0.);
    }

    // what the exchange does after sending
    void exchanged() { patchGhost.clear(); }
};

Inflow_t makeInflow()
{
    return Inflow_t{1., {0., 0., 0.}, {.1, .1, .1}, ParticleBoundaryInfo{"protons", 1, 1}};
}


std::size_t countTracked(auto const& particles, auto&& where)
{
    std::size_t n = 0;
    for (auto const& p : particles)
        if (p.weight == tracked and where(p.iCell))
            ++n;
    return n;
}
std::size_t countTracked(auto const& particles)
{
    return countTracked(particles, [](auto const&) { return true; });
}

// the stand-in for the particle exchange: copies of `from`'s patch ghosts that land in `to`
std::size_t exchanged(Patch const& from, Box_t const& to)
{
    return countTracked(from.patchGhost, [&](auto const& cell) { return isIn(cell, to); });
}

void flipVy(ParticleArray& particles)
{
    for (auto& p : particles)
        if (p.weight == tracked)
            p.v[1] = -p.v[1];
}

} // namespace



// The stock particle is drawn in front of P and enters Q's cells diagonally. It must end up in
// one domain only, once.
TEST(StockParticleNearASeam, entersTheNeighbourOnce)
{
    Patch P{boxP, boxQ};
    Patch Q{boxQ, boxP};

    // cell [-1, 4]: in P's ownBox and in Q's overhang. (-1, 4) + (0.8, 0.8) + (0.5, 0.5) = (0.3, 5.3)
    Particle_t const s{tracked, {-1, 4}, {0.8, 0.8}, {v, v, 0.}};
    P.levelGhost.push_back(s);
    Q.levelGhost.push_back(s);

    P.push();
    Q.push();

    // Q owns cell [0, 5]: it takes the particle into its domain
    EXPECT_EQ(1u, countTracked(Q.domain));
    EXPECT_EQ(0u, countTracked(P.domain));
    EXPECT_EQ(0u, exchanged(P, boxQ)); // P did not send it: it was never P's domain particle

    // but the updater keeps P's copy: [0, 5] is in P's ghost layer. Refill would not erase it
    // there: it is not in the face ghost box
    EXPECT_EQ(1u, countTracked(P.levelGhost));
    for (auto const& p : P.levelGhost)
        if (p.weight == tracked)
            EXPECT_FALSE(isIn(p.iCell, P.faceGhostBox()));

    // the inflow rule erases it at the Leavers moment
    auto inflow = makeInflow();
    P.leavers(inflow);
    Q.leavers(inflow);
    EXPECT_EQ(0u, countTracked(P.levelGhost));
    EXPECT_EQ(1u, countTracked(Q.domain));
}



// Same start, then the particle turns back into P (here: vy flips for every copy, as a field
// would do it to all of them). P then gets the particle twice: from Q through the exchange, and
// from its own stale stock copy, exported to its domain by the updater.
TEST(StockParticleNearASeam, turnsBackIntoThePatchItCameFromAndIsCountedOnce)
{
    Patch P{boxP, boxQ};
    Patch Q{boxQ, boxP};

    Particle_t const s{tracked, {-1, 4}, {0.8, 0.8}, {v, v, 0.}};
    P.levelGhost.push_back(s);
    Q.levelGhost.push_back(s);

    // step 1: into Q's cell [0, 5]. Then the Leavers moment, then the exchange
    P.push();
    Q.push();
    auto inflow = makeInflow();
    P.leavers(inflow);
    Q.leavers(inflow);
    P.exchanged();
    Q.exchanged();

    // step 2: back down. (0.3, 5.3) + (0.5, -0.5) = (0.8, 4.8): cell [0, 4], owned by P
    flipVy(P.levelGhost);
    flipVy(Q.domain);
    P.push();
    Q.push();

    auto const inP = countTracked(P.domain) + exchanged(Q, boxP);
    auto const inQ = countTracked(Q.domain) + exchanged(P, boxQ);

    EXPECT_EQ(1u, inP + inQ) << "the particle exists " << inP + inQ << " times after the exchange"
                             << " (P: " << countTracked(P.domain) << " from its stale stock copy + "
                             << exchanged(Q, boxP) << " from Q)";
}


// Open face: the stock is made of leavers, held by one patch only. A leaver drawn back across the
// face into the neighbour's cells must reach the neighbour, once.
TEST(StockParticleNearASeam, openLeaverEnteringTheNeighbourIsHandedToIt)
{
    Patch P{boxP, boxQ};
    Patch Q{boxQ, boxP};

    // a leaver in P's stock, in front of P, moving back in and up: (0.3, 5.3), Q's cell [0, 5]
    P.levelGhost.push_back(Particle_t{tracked, {-1, 4}, {0.8, 0.8}, {v, v, 0.}});

    P.push();
    Q.push();
    EXPECT_EQ(1u, countTracked(P.levelGhost)); // the updater kept it: a stray

    Open_t open;
    P.leavers(open);
    Q.leavers(open);

    EXPECT_EQ(0u, countTracked(P.levelGhost));
    EXPECT_EQ(1u, exchanged(P, boxQ)); // the exchange gives it to Q
    EXPECT_EQ(0u, countTracked(P.domain) + countTracked(Q.domain) + exchanged(Q, boxP));
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
