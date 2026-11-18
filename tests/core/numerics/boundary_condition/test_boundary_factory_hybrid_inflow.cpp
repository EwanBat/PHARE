// What the BoundaryFactory builds for a hybrid inflow face.
//
//     fields      rho  -> None            (the particles give the moments)
//                 B    -> DivergenceFreeTransverseDirichlet, with the B of the face
//                 E    -> Dirichlet, E = -V_q x B: B of the face, V_q the charge-weighted
//                         velocity of the wind, sum q n V / sum q n (windChargeVelocity)
//     particles   a population listed in "populations" -> Inflow, with its own n, V, vth
//                 a population not listed               -> Open
//
// The dict is built as pyphare/pharein/boundary.py writes it (HybridInflowBoundary.populate_dict).

#include "core/boundary/boundary_factory.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_inflow.hpp"
#include "initializer/data_provider.hpp"

#include "hybrid_bc_test_fixtures.hpp"

#include "gtest/gtest.h"

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace PHARE::core;
using PHARE::initializer::PHAREDict;

namespace
{
using ParticleArrayT = ParticleArray<1>;
using Factory        = BoundaryFactory<GridLayout1D, Field1D, ParticleArrayT>;
using Inflow_t       = ParticleBoundaryConditionInflow<ParticleArrayT, Field1D, GridLayout1D>;
using Scalar         = HybridQuantity::Scalar;
using Vector         = HybridQuantity::Vector;
using FBC            = FieldBoundaryConditionType;
using PBC            = ParticleBoundaryConditionType;

// what HybridModel registers
std::vector<Scalar> const scalars{Scalar::rho};
std::vector<Vector> const vectors{Vector::B, Vector::E};

// all components non-zero, so a wrong sign or a swapped component shows
std::array<double, 3> const faceB{1.0, 0.4, 0.7};

// two wind populations with different velocities: V_q is neither of them
double const protonN = 1.0, alphaN = 0.04;
double const protonQ = 1.0, alphaQ = 2.0;
std::array<double, 3> const protonV{2.0, 0.5, -0.3};
std::array<double, 3> const protonVth{0.3, 0.2, 0.1};
std::array<double, 3> const alphaV{2.2, -0.4, 0.6};
std::array<double, 3> const alphaVth{0.1, 0.1, 0.1};

// the populations of the model: the wind is protons + alphas, the beam is not in the wind
std::vector<ParticleBoundaryInfo> const populations{
    {.name = "protons", .nbrParticlesPerCell = 100, .seed = 1, .charge = protonQ},
    {.name = "alpha", .nbrParticlesPerCell = 50, .seed = 2, .charge = alphaQ},
    {.name = "beam", .nbrParticlesPerCell = 20, .seed = 3, .charge = 1.0}};


void addVector(PHAREDict& dict, std::string const& key, std::array<double, 3> const& v)
{
    dict[key]["x"] = v[0];
    dict[key]["y"] = v[1];
    dict[key]["z"] = v[2];
}

void addPopulation(PHAREDict& dict, std::string const& name, double const density,
                   std::array<double, 3> const& v, std::array<double, 3> const& vth)
{
    auto& pop      = dict["populations"][name];
    pop["density"] = density;
    addVector(pop, "velocity", v);
    addVector(pop, "vth", vth);
}

PHAREDict inflowDict()
{
    PHAREDict dict;
    dict["type"] = static_cast<int>(BoundaryType::HybridInflow);
    addVector(dict, "B", faceB); // pyphare writes no face velocity: one per population
    addPopulation(dict, "protons", protonN, protonV, protonVth);
    addPopulation(dict, "alpha", alphaN, alphaV, alphaVth);
    return dict;
}

auto create(PHAREDict const& dict)
{
    return Factory::create(BoundaryLocation::XLower, dict, scalars, vectors, 0.0, populations);
}

// V_q = sum q n V / sum q n, written here from its definition
std::array<double, 3> chargeVelocity()
{
    double const qp = protonQ * protonN, qa = alphaQ * alphaN;
    std::array<double, 3> v;
    for (std::size_t d = 0; d < 3; ++d)
        v[d] = (qp * protonV[d] + qa * alphaV[d]) / (qp + qa);
    return v;
}

// -V x B
std::array<double, 3> motionalE(std::array<double, 3> const& v, std::array<double, 3> const& B)
{
    return {-(v[1] * B[2] - v[2] * B[1]), -(v[2] * B[0] - v[0] * B[2]),
            -(v[0] * B[1] - v[1] * B[0])};
}

} // namespace



TEST(HybridInflowFactory, registersTheFieldConditions)
{
    auto boundary = create(inflowDict());
    ASSERT_NE(boundary, nullptr);
    EXPECT_EQ(BoundaryType::HybridInflow, boundary->getType());

    auto rho = boundary->getFieldCondition(Scalar::rho);
    auto B   = boundary->getFieldCondition(Vector::B);
    auto E   = boundary->getFieldCondition(Vector::E);
    ASSERT_NE(rho, nullptr);
    ASSERT_NE(B, nullptr);
    ASSERT_NE(E, nullptr);
    EXPECT_EQ(FBC::None, rho->getType());
    EXPECT_EQ(FBC::DivergenceFreeTransverseDirichlet, B->getType());
    EXPECT_EQ(FBC::Dirichlet, E->getType());
}


TEST(HybridInflowFactory, imposesTheMotionalElectricFieldInTheGhosts)
{
    auto boundary = create(inflowDict());
    auto bc       = boundary->getFieldCondition(Vector::E);
    ASSERT_NE(bc, nullptr);

    GridLayout1D layout{{0.1}, {nCells}, {0.0}};
    UsableTensorField<1, 1> E{"E", layout, Vector::E};
    for (std::size_t c = 0; c < 3; ++c)
        for (std::uint32_t i = 0; i < E[c].shape()[0]; ++i)
            E[c](i) = ghostSentinel;

    bc->apply(E.super(), BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);

    auto const expected = motionalE(chargeVelocity(), faceB);
    for (std::size_t c = 0; c < 3; ++c)
    {
        auto const qty   = HybridQuantity::componentsQuantities(Vector::E)[c];
        auto const start = layout.physicalStartIndex(qty, Direction::X);
        for (std::uint32_t i = 0; i < start; ++i)
            EXPECT_DOUBLE_EQ(expected[c], E[c](i)) << "component " << c << ", ghost node " << i;
    }
}


TEST(HybridInflowFactory, windPopulationsGetInflowWithTheirOwnPlasma)
{
    auto boundary = create(inflowDict());

    struct Wind
    {
        std::string name;
        double density;
        std::array<double, 3> v, vth;
        std::uint32_t ppc;
    };
    for (auto const& wind : {Wind{"protons", protonN, protonV, protonVth, 100},
                             Wind{"alpha", alphaN, alphaV, alphaVth, 50}})
    {
        auto bc = boundary->getParticleCondition(wind.name);
        ASSERT_NE(bc, nullptr) << wind.name;
        EXPECT_EQ(PBC::Inflow, bc->getType()) << wind.name;

        auto inflow = std::dynamic_pointer_cast<Inflow_t>(bc);
        ASSERT_NE(inflow, nullptr) << wind.name;
        EXPECT_DOUBLE_EQ(wind.density, inflow->density()) << wind.name;
        // each population is drawn at its own velocity, not at V_q
        EXPECT_EQ(wind.v, inflow->velocity()) << wind.name;
        EXPECT_EQ(wind.vth, inflow->thermalVelocity()) << wind.name;
        EXPECT_EQ(wind.ppc, inflow->nbrParticlesPerCell()) << wind.name; // from the population
    }
}


TEST(HybridInflowFactory, populationsNotInTheWindGetOpen)
{
    auto boundary = create(inflowDict());
    auto bc       = boundary->getParticleCondition("beam");
    ASSERT_NE(bc, nullptr);
    EXPECT_EQ(PBC::Open, bc->getType());
}


TEST(HybridInflowFactory, refusesAWindPopulationTheModelDoesNotHave)
{
    auto dict = inflowDict();
    addPopulation(dict, "protns", 1.0, protonV, protonVth); // a typo
    EXPECT_THROW(create(dict), std::runtime_error);
}


TEST(HybridInflowFactory, refusesAMissingKey)
{
    for (auto const key : {"B", "populations"})
    {
        PHAREDict dict;
        dict["type"] = static_cast<int>(BoundaryType::HybridInflow);
        if (std::string{key} != "B")
            addVector(dict, "B", faceB);
        if (std::string{key} != "populations")
            addPopulation(dict, "protons", 1.0, protonV, protonVth);
        EXPECT_THROW(create(dict), std::runtime_error) << "missing " << key;
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
