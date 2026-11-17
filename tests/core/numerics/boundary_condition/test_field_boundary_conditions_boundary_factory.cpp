#include "gtest/gtest.h"

#include "core/boundary/boundary_factory.hpp"
#include "core/numerics/hybrid_inflow/hybrid_inflow_wind.hpp"
#include "initializer/data_provider.hpp"

#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

using namespace PHARE::core;
using PHARE::initializer::PHAREDict;


namespace
{
using Factory = BoundaryFactory<GridLayoutMHD1D, FieldMHD<1>>;
using Scalar  = MHDQuantity::Scalar;
using Vector  = MHDQuantity::Vector;
using FBC     = FieldBoundaryConditionType;

double constexpr heatCapacityRatio = 5.0 / 3.0;

std::vector<Scalar> const mhdScalars{Scalar::rho, Scalar::Etot};
std::vector<Vector> const mhdVectors{Vector::B, Vector::E, Vector::rhoV};

void fillInflowData(PHAREDict& dict)
{
    dict["density"]  = 1.0;
    dict["pressure"] = 2.0;

    dict["velocity"]["x"] = 3.0;
    dict["velocity"]["y"] = 0.0;
    dict["velocity"]["z"] = 0.0;

    dict["B"]["x"] = 0.75;
    dict["B"]["y"] = 1.0;
    dict["B"]["z"] = 0.0;
}

PHAREDict dictFor(BoundaryType const type)
{
    PHAREDict dict;
    dict["type"] = static_cast<int>(type);
    if (type == BoundaryType::SuperMagnetofastInflow)
        fillInflowData(dict);
    return dict;
}

struct ExpectedDispatch
{
    std::map<Scalar, FBC> scalars;
    std::map<Vector, FBC> vectors;
};

void checkDispatch(BoundaryType const type, ExpectedDispatch const& expected)
{
    auto boundary = Factory::create(BoundaryLocation::XLower, dictFor(type), mhdScalars, mhdVectors,
                                    heatCapacityRatio);
    ASSERT_NE(boundary, nullptr) << "factory returned null for type " << static_cast<int>(type);

    for (auto const& [qty, expectedType] : expected.scalars)
    {
        auto bc = boundary->getFieldCondition(qty);
        ASSERT_NE(bc, nullptr) << "type " << static_cast<int>(type) << " left scalar quantity "
                               << static_cast<int>(qty) << " with no condition";
        EXPECT_EQ(bc->getType(), expectedType)
            << "type " << static_cast<int>(type) << ", scalar quantity " << static_cast<int>(qty)
            << ": got " << static_cast<int>(bc->getType()) << ", expected "
            << static_cast<int>(expectedType);
    }
    for (auto const& [qty, expectedType] : expected.vectors)
    {
        auto bc = boundary->getFieldCondition(qty);
        ASSERT_NE(bc, nullptr) << "type " << static_cast<int>(type) << " left vector quantity "
                               << static_cast<int>(qty) << " with no condition";
        EXPECT_EQ(bc->getType(), expectedType)
            << "type " << static_cast<int>(type) << ", vector quantity " << static_cast<int>(qty)
            << ": got " << static_cast<int>(bc->getType()) << ", expected "
            << static_cast<int>(expectedType);
    }
}
} // namespace

TEST(BoundaryFactory, NoneLeavesEveryQuantityUntouched)
{
    checkDispatch(BoundaryType::None,
                  {{{Scalar::rho, FBC::None}, {Scalar::Etot, FBC::None}},
                   {{Vector::B, FBC::None}, {Vector::E, FBC::None}, {Vector::rhoV, FBC::None}}});
}

TEST(BoundaryFactory, Reflective)
{
    checkDispatch(BoundaryType::Reflective,
                  {{{Scalar::rho, FBC::Neumann}, {Scalar::Etot, FBC::Neumann}},
                   {{Vector::B, FBC::DivergenceFreeTransverseNeumann},
                    {Vector::E, FBC::AntiSymmetric},
                    {Vector::rhoV, FBC::Symmetric}}});
}

TEST(BoundaryFactory, Open)
{
    checkDispatch(BoundaryType::Open, {{{Scalar::rho, FBC::Neumann}, {Scalar::Etot, FBC::Neumann}},
                                       {{Vector::B, FBC::DivergenceFreeTransverseNeumann},
                                        {Vector::E, FBC::None},
                                        {Vector::rhoV, FBC::Neumann}}});
}

TEST(BoundaryFactory, SuperMagnetofastInflow)
{
    checkDispatch(BoundaryType::SuperMagnetofastInflow,
                  {{{Scalar::rho, FBC::Dirichlet}, {Scalar::Etot, FBC::Dirichlet}},
                   {{Vector::B, FBC::DivergenceFreeTransverseDirichlet},
                    {Vector::E, FBC::Dirichlet},
                    {Vector::rhoV, FBC::Dirichlet}}});
}

TEST(BoundaryFactory, SuperMagnetofastInflowImposesTotalEnergyFromInflowState)
{
    auto boundary
        = Factory::create(BoundaryLocation::XLower, dictFor(BoundaryType::SuperMagnetofastInflow),
                          mhdScalars, mhdVectors, heatCapacityRatio);
    auto bc = boundary->getFieldCondition(Scalar::Etot);
    ASSERT_NE(bc, nullptr);

    double const rho = 1.0, P = 2.0, vx = 3.0, Bx = 0.75, By = 1.0;
    double const expected
        = P / (heatCapacityRatio - 1.0) + 0.5 * rho * vx * vx + 0.5 * (Bx * Bx + By * By);

    GridLayoutMHD1D layout{{0.1}, {nCellsMHD}, {0.0}};
    GridMHD1D grid{"Etot", Scalar::Etot, layout.allocSize(Scalar::Etot)};
    FieldMHD<1>& Etot{*(&grid)};
    for (std::uint32_t i = 0; i < Etot.shape()[0]; ++i)
        Etot(i) = expected;
    for (std::uint32_t i = 0; i < mhdGhostWidth; ++i)
        Etot(i) = -1.0;

    bc->apply(Etot, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, 0.0);

    for (std::uint32_t i = 0; i < mhdGhostWidth; ++i)
        EXPECT_DOUBLE_EQ(Etot(i), expected);
}

TEST(BoundaryFactory, SuperMagnetofastInflowImposesMotionalElectricField)
{
    auto boundary
        = Factory::create(BoundaryLocation::XLower, dictFor(BoundaryType::SuperMagnetofastInflow),
                          mhdScalars, mhdVectors, heatCapacityRatio);
    auto bc = boundary->getFieldCondition(Vector::E);
    ASSERT_NE(bc, nullptr);

    double const vx = 3.0, Bx = 0.75, By = 1.0;
    std::array<double, 3> const expected{0.0, 0.0, Bx * 0.0 - By * vx};

    GridLayoutMHD1D layout{{0.1}, {nCellsMHD}, {0.0}};
    UsableVecFieldMHD<1> Evec{"bc_test_E", layout, Vector::E};
    auto& E = Evec.super();

    for (std::size_t c = 0; c < 3; ++c)
    {
        auto& Ec = E[c];
        for (std::uint32_t i = 0; i < Ec.shape()[0]; ++i)
            Ec(i) = (i < mhdGhostWidth) ? -1.0 : expected[c];
    }

    bc->apply(E, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, 0.0);

    for (std::size_t c = 0; c < 3; ++c)
        for (std::uint32_t i = 0; i <= mhdGhostWidth; ++i)
            EXPECT_DOUBLE_EQ(E[c](i), expected[c]) << "component " << c << " index " << i;
}


namespace
{
void addWindPopulation(PHAREDict& wind, std::string const& name, double const density,
                       std::array<double, 3> const& v)
{
    wind[name]["density"]       = density;
    wind[name]["velocity"]["x"] = v[0];
    wind[name]["velocity"]["y"] = v[1];
    wind[name]["velocity"]["z"] = v[2];
}

ParticleBoundaryInfo population(std::string const& name, double const charge)
{
    return {.name = name, .nbrParticlesPerCell = 10, .seed = 1, .charge = charge};
}
} // namespace

// E at a hybrid inflow face uses V_q = sum q n V / sum q n, the electron velocity of Ohm's law
// in a uniform wind (J = 0)
TEST(WindChargeVelocity, isTheChargeWeightedVelocityOfTheWind)
{
    PHAREDict wind;
    addWindPopulation(wind, "protons", 1.0, {2.0, 0.0, 0.0});
    addWindPopulation(wind, "alpha", 0.05, {2.5, 0.1, 0.0});
    std::vector<ParticleBoundaryInfo> const pops{population("protons", 1.0),
                                                 population("alpha", 2.0)};

    auto const v = windChargeVelocity(wind, pops);

    double const qn = 1.0 * 1.0 + 2.0 * 0.05;
    EXPECT_DOUBLE_EQ((1.0 * 1.0 * 2.0 + 2.0 * 0.05 * 2.5) / qn, v[0]);
    EXPECT_DOUBLE_EQ((2.0 * 0.05 * 0.1) / qn, v[1]);
    EXPECT_DOUBLE_EQ(0.0, v[2]);
}

TEST(WindChargeVelocity, isTheCommonVelocityWhenAllPopulationsShareIt)
{
    PHAREDict wind;
    addWindPopulation(wind, "protons", 1.0, {2.0, 0.3, 0.0});
    addWindPopulation(wind, "alpha", 0.05, {2.0, 0.3, 0.0});
    std::vector<ParticleBoundaryInfo> const pops{population("protons", 1.0),
                                                 population("alpha", 2.0)};

    auto const v = windChargeVelocity(wind, pops);

    EXPECT_DOUBLE_EQ(2.0, v[0]);
    EXPECT_DOUBLE_EQ(0.3, v[1]);
    EXPECT_DOUBLE_EQ(0.0, v[2]);
}

TEST(WindChargeVelocity, ignoresThePopulationsNotInTheWind)
{
    PHAREDict wind;
    addWindPopulation(wind, "protons", 1.0, {2.0, 0.0, 0.0});
    std::vector<ParticleBoundaryInfo> const pops{population("protons", 1.0),
                                                 population("beam", 1.0)};

    EXPECT_DOUBLE_EQ(2.0, windChargeVelocity(wind, pops)[0]);
}

TEST(WindChargeVelocity, refusesAWindWithoutCharge)
{
    PHAREDict wind;
    addWindPopulation(wind, "protons", 1.0, {2.0, 0.0, 0.0});
    std::vector<ParticleBoundaryInfo> const pops{population("protons", 0.0)};

    EXPECT_THROW(windChargeVelocity(wind, pops), std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
