// How much plasma an inflow face lets in, and with which velocities. 1D, real IonUpdater.
//
// Root level, order 1 (particle ghost width 1). One patch, cells [0, 59], empty at the start.
// Inflow face x-lower:
//
//        stock │ domain (empty at t = 0)
//         [-1] │ [0] [1] ... [59]
//              x = 0
//
// Each step calls the inflow rule as SolverPPC does: PredictedMoments (after the push of copies,
// which draws the stock of the new time and keeps it), then the real push, then Refill (which
// reuses that kept stock), Moments, Leavers. E = B = 0, so the particles go in straight lines,
// and none leaves the domain during the run.
//
// The physics: through a face, a drifting Maxwellian (density n, drift u along the normal,
// thermal velocity vth along each axis) sends per unit time and unit area
//
//     Gamma = n [ vth phi(u/vth) + u Phi(u/vth) ]        phi, Phi: normal density and its integral
//
// which is more than n u: the thermal motion carries particles in too. The particles that come in
// are not a Maxwellian: a fast particle crosses more often, so the normal velocities of those that
// enter follow v f(v) for v > 0 ("flux-weighted"), with mean
//
//     <v_n> = [ (u^2 + vth^2) Phi + u vth phi ] / [ u Phi + vth phi ]
//
// The stock is one cell thick: a particle needs |v| dt < dx to cross it in one step. Here
// dx / dt = 10 vth, so the missing tail is far below the statistical error.

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_inflow.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/utilities/box/box.hpp"

#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <string>
#include <vector>

using namespace PHARE::core;

namespace
{
constexpr std::size_t dim          = 1;
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
using Field_t       = Grid::field_type;
using Inflow_t      = ParticleBoundaryConditionInflow<ParticleArray, Field_t, GridLayout>;
using Stage         = ParticleBoundaryStage;

constexpr double dx = 0.1;
constexpr double dt = 0.01;

// the wind
constexpr double n   = 1.;
constexpr double u   = 0.5; // drift along x, half the thermal velocity: the thermal part matters
constexpr double vth = 1.;
constexpr std::uint32_t ppc = 10000;

constexpr std::size_t nSteps = 20;

Box_t const patchBox{{0}, {59}};
Box_t const stockBox{{-1}, {-1}};


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


// normal density and its integral
double phi(double const z)
{
    return std::exp(-0.5 * z * z) / std::sqrt(2. * std::numbers::pi);
}
double Phi(double const z)
{
    return 0.5 * (1. + std::erf(z / std::numbers::sqrt2));
}

// particles per unit time per unit area through the face
double expectedFlux()
{
    return n * (vth * phi(u / vth) + u * Phi(u / vth));
}

// mean normal velocity of the particles that come in
double expectedMeanVx()
{
    double const z = u / vth;
    return ((u * u + vth * vth) * Phi(z) + u * vth * phi(z)) / (u * Phi(z) + vth * phi(z));
}


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
    Inflow_t inflow{n, {u, 0., 0.}, {vth, vth, vth}, ParticleBoundaryInfo{"protons", ppc, 42}};

    Patch()
        : layout{{dx}, {static_cast<std::uint32_t>(patchBox.shape()[0])}, {0.}, patchBox}
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
        , domain{grow(patchBox, mapMargin)}
        , patchGhost{grow(patchBox, mapMargin)}
        , levelGhost{grow(patchBox, mapMargin)}
        , levelGhostOld{grow(patchBox, mapMargin)}
        , levelGhostNew{grow(patchBox, mapMargin)}
        , pack{"protons", &domain, &patchGhost, &levelGhost, &levelGhostOld, &levelGhostNew}
        // no neighbour: nonLevelGhostBox is the patch itself
        , boxing{layout, {patchBox}, {}}
    {
        E.set_on(em.E);
        B.set_on(em.B);

        auto const& [V, m, cd, md, cf] = ions.getCompileTimeResourcesViewList();
        Vi.set_on(V);
        chargeFlux.set_on(cf);
        M.set_on(m);
        cd.setBuffer(&ionChargeDensity);
        md.setBuffer(&ionMassDensity);

        auto& pops                             = ions.getRunTimeResourcesViewList();
        auto const& [F, popm, d, c, particles] = pops[0].getCompileTimeResourcesViewList();
        d.setBuffer(&density);
        c.setBuffer(&chargeDensity);
        popM.set_on(popm);
        flux.set_on(F);
        particles.setBuffer(&pack);
    }

    void rule(Stage const stage, double const time)
    {
        Inflow_t::data_type data{domain, patchGhost, levelGhost, *density, flux};
        inflow.apply(stage, data, BoundaryLocation::XLower, stockBox, layout, time);
    }

    // one step as SolverPPC does it, without the exchange (no neighbour). The push of copies
    // (domain_only) only deposits moments: it is left out, but its PredictedRefill is kept
    // because it draws the stock that Refill then uses.
    void realStep(double const newTime)
    {
        rule(Stage::PredictedRefill, newTime);
        IonUpdater_t updater{dict["simulation"]["algo"]["ion_updater"]};
        updater.updatePopulations(ions, em, boxing, dt, UpdaterMode::all);
        rule(Stage::Refill, newTime);
        rule(Stage::Moments, newTime);
        rule(Stage::Leavers, newTime);
        patchGhost.clear(); // the exchange empties it
    }
};

} // namespace



TEST(InflowFlux, theFaceLetsInTheThermalFluxOfADriftingMaxwellian)
{
    Patch P;
    P.rule(Stage::Init, 0.);

    std::vector<std::size_t> enteredPerStep;
    for (std::size_t step = 1; step <= nSteps; ++step)
    {
        auto const before = P.domain.size();
        P.realStep(step * dt);
        enteredPerStep.push_back(P.domain.size() - before);
    }

    // with E = B = 0 nothing leaves: every domain particle came in through the face
    double weight = 0.;
    std::array<double, 3> mean{}, var{};
    for (auto const& p : P.domain)
    {
        weight += p.weight;
        for (std::size_t i = 0; i < 3; ++i)
            mean[i] += p.v[i] / P.domain.size();
    }
    for (auto const& p : P.domain)
        for (std::size_t i = 0; i < 3; ++i)
            var[i] += std::pow(p.v[i] - mean[i], 2) / P.domain.size();

    // a weight is a density times a cell: weight * dx is a number of physical particles (per
    // unit area). About 14000 particles came in: the statistical error is about 1%.
    double const entered  = weight * dx;
    double const expected = expectedFlux() * nSteps * dt;
    EXPECT_NEAR(expected, entered, 0.04 * expected)
        << "the face lets in " << entered / (nSteps * dt) << " per unit time, expected "
        << expectedFlux() << " (n u alone would give " << n * u << ")";

    // the normal velocities are flux-weighted: much faster on average than the drift u
    EXPECT_NEAR(expectedMeanVx(), mean[0], 0.03)
        << "drift " << u << ", mean over the half Maxwellian v > 0 "
        << u + vth * phi(u / vth) / Phi(u / vth);

    // the transverse velocities are not changed by crossing the face
    for (std::size_t i = 1; i < 3; ++i)
    {
        EXPECT_NEAR(0., mean[i], 0.03) << "component " << i;
        EXPECT_NEAR(vth, std::sqrt(var[i]), 0.03) << "component " << i;
    }

    // the stock is new at each step (the time is in the key of the draw): the same stock would
    // let in the same number of particles at every step
    std::size_t sameAsFirst = 0;
    for (auto const count : enteredPerStep)
        sameAsFirst += (count == enteredPerStep[0]);
    EXPECT_LT(sameAsFirst, nSteps);
}


// With E = B = 0 each particle moves freely. A particle of velocity v > 0 that came in at any
// time since t = 0 is now in (0, v t], and the face lets them in uniformly in time, so at x only
// those with v > x / t are there:
//
//     n(x, t) = n * P(v > x / t) = n Phi((u - x / t) / vth)
//
//     n ┤ n Phi(u/vth) at the face (only v > 0 come in: less than n)
//       │‾‾‾‾‾‾╲
//       │       ╲___
//       └──────────────── x        a front of width ~ vth t
//
// This holds step by step, not only for dt -> 0: a particle that crosses at step k is uniform in
// (v (t - t_k), v (t - t_k-1)], and these intervals tile (0, v t].
TEST(InflowFlux, theDomainFillsWithTheFreeStreamingProfile)
{
    Patch P;
    P.rule(Stage::Init, 0.);
    for (std::size_t step = 1; step <= nSteps; ++step)
        P.realStep(step * dt);
    double const t = nSteps * dt;

    // the domain particles only (the updater also deposited the stock, beyond x = 0). A weight is
    // a density times a cell, so the weights of a cell add up to its mean density.
    std::size_t constexpr nCells = 8;
    std::vector<double> perCell(nCells, 0.);
    for (auto const& p : P.domain)
        if (static_cast<std::size_t>(p.iCell[0]) < nCells)
            perCell[p.iCell[0]] += p.weight;

    double const weight = n / ppc;
    for (std::size_t i = 0; i < nCells; ++i)
    {
        // mean of the profile over the cell
        std::size_t constexpr nSub = 100;
        double expected            = 0.;
        for (std::size_t k = 0; k < nSub; ++k)
        {
            double const x = (i + (k + 0.5) / nSub) * dx;
            expected += n * Phi((u - x / t) / vth) / nSub;
        }
        // counting noise: sqrt(number of particles) particles of weight `weight`
        double const tolerance = 5. * std::sqrt(expected * weight) + 1e-3;
        EXPECT_NEAR(expected, perCell[i], tolerance) << "cell " << i;
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
