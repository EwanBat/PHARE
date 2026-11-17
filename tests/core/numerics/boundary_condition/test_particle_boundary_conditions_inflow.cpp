// Tests of the inflow particle boundary condition: the Moments and PredictedMoments stages.

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"
#include "core/numerics/boundary_condition/particle_boundary_condition_inflow.hpp"

#include "hybrid_bc_test_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>

using namespace PHARE::core;


TEST(FrontOfOwnCells, cutsAFaceBoxToThePatchAlongTheFace)
{
    Box<int, 2> const domain{{0, 0}, {9, 7}};
    // x-lower face, ghost width 1: it also covers the cells in front of the neighbours (y=-1, 8)
    Box<int, 2> const faceBox{{-1, -1}, {-1, 8}};

    EXPECT_EQ((Box<int, 2>{{-1, 0}, {-1, 7}}), frontOfOwnCells(faceBox, domain));
}

TEST(FrontOfOwnCells, keepsACornerBoxWhole)
{
    Box<int, 2> const domain{{0, 0}, {9, 7}};
    Box<int, 2> const cornerBox{{-2, -2}, {-1, -1}}; // beyond x-lower and y-lower

    EXPECT_EQ(cornerBox, frontOfOwnCells(cornerBox, domain));
}



TEST(BoxMinus, givesTheCellsOutsideTheHoleOnce)
{
    Box<int, 2> const box{{0, 0}, {9, 9}};
    Box<int, 2> const hole{{3, 4}, {5, 6}};

    auto const parts = boxMinus(box, hole);

    std::size_t cells = 0;
    for (auto const& part : parts)
    {
        cells += part.size();
        EXPECT_FALSE(part * hole); // no cell of the hole
        for (auto const& other : parts)
            if (&other != &part)
                EXPECT_FALSE(part * other); // no cell twice
    }
    EXPECT_EQ(box.size() - hole.size(), cells);
}

TEST(BoxMinus, keepsTheBoxWhenTheHoleIsElsewhere)
{
    Box<int, 2> const box{{0, 0}, {3, 3}};
    auto const parts = boxMinus(box, Box<int, 2>{{10, 10}, {12, 12}});
    ASSERT_EQ(1u, parts.size());
    EXPECT_EQ(box, parts[0]);
}



TEST(StrayBoxes, areTheCellsInsideInFrontOfTheNeighbours)
{
    Box<int, 2> const domain{{0, 0}, {9, 7}};

    // x-lower face, ghost width 2: x in [-2, -1], y in front of the patch and 2 cells beyond
    auto const lower = strayBoxes(Box<int, 2>{{-2, -2}, {-1, 9}}, domain);
    ASSERT_EQ(2u, lower.size());
    EXPECT_EQ((Box<int, 2>{{0, -2}, {1, -1}}), lower[0]);
    EXPECT_EQ((Box<int, 2>{{0, 8}, {1, 9}}), lower[1]);

    // x-upper face, ghost width 1
    auto const upper = strayBoxes(Box<int, 2>{{10, -1}, {10, 8}}, domain);
    ASSERT_EQ(2u, upper.size());
    EXPECT_EQ((Box<int, 2>{{9, -1}, {9, -1}}), upper[0]);
    EXPECT_EQ((Box<int, 2>{{9, 8}, {9, 8}}), upper[1]);
}

TEST(StrayBoxes, noneForACornerBox)
{
    Box<int, 2> const domain{{0, 0}, {9, 7}};
    EXPECT_TRUE(strayBoxes(Box<int, 2>{{-2, -2}, {-1, -1}}, domain).empty());
}

TEST(StrayBoxes, anEdgeBoxIn3DCanHaveStraysAlongItsLength)
{
    Box<int, 3> const domain{{0, 0, 0}, {9, 7, 5}};
    // beyond x-lower and y-lower (codim 2), along z also in front of the z-neighbours
    auto const strays = strayBoxes(Box<int, 3>{{-1, -1, -1}, {-1, -1, 6}}, domain);
    ASSERT_EQ(2u, strays.size());
    EXPECT_EQ((Box<int, 3>{{0, 0, -1}, {0, 0, -1}}), strays[0]);
    EXPECT_EQ((Box<int, 3>{{0, 0, 6}, {0, 0, 6}}), strays[1]);
}



struct InflowMoments2D : testing::Test
{
    using ParticleArrayT = ParticleArray<2>;
    using InflowT        = ParticleBoundaryConditionInflow<ParticleArrayT, Field2D, GridLayout2D>;
    using data_type      = InflowT::data_type;

    static constexpr double density = 2.;
    static constexpr std::uint32_t ppc = 50;
    static constexpr double time       = 0.37;

    GridLayout2D layout{{0.1, 0.1}, {nCellsX2D, nCellsY2D}, {0.0, 0.0}};
    Box<int, 2> const particleBox = grow(layout.AMRBox(), 1);

    // x-lower face, with the cells in front of the neighbours along y
    Box<int, 2> const ghostBox{{-1, -1}, {-1, static_cast<int>(nCellsY2D)}};

    InflowT inflow{density, {1., 0., 0.}, {.3, .3, .3}, ParticleBoundaryInfo{"protons", ppc, 42}};

    struct Moments
    {
        Moments(GridLayout2D const& layout)
            : rho{"rho", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)}
            , F{"F", layout, HybridQuantity::Vector::V}
        {
            rho.zero();
            for (std::size_t c = 0; c < 3; ++c)
                F[c].zero();
        }
        Grid2D rho;
        UsableTensorField<2, 1> F;
    };

    struct Arrays
    {
        Arrays(Box<int, 2> const& box)
            : domain{box}
            , patchGhost{box}
            , stock{box}
        {
        }
        ParticleArrayT domain, patchGhost, stock;
    };

    data_type dataOf(Arrays& a, Moments& m)
    {
        return data_type{a.domain, a.patchGhost, a.stock, *(&m.rho), m.F};
    }

    void apply(ParticleBoundaryStage stage, Arrays& a, Moments& m)
    {
        auto data = dataOf(a, m);
        inflow.apply(stage, data, BoundaryLocation::XLower, ghostBox, layout, time);
    }
};


TEST_F(InflowMoments2D, predictedMomentsDepositsWhatMomentsDepositsAfterARefill)
{
    // all push: the stock is drawn at `time` (here by Init, as Refill would), then deposited
    Arrays refilled{particleBox};
    Moments fromStock{layout};
    apply(ParticleBoundaryStage::Init, refilled, fromStock);
    apply(ParticleBoundaryStage::Moments, refilled, fromStock);

    // push of copies: the real stock is empty here (it would be at the start of the step)
    Arrays notMoved{particleBox};
    Moments fromTempDraw{layout};
    apply(ParticleBoundaryStage::PredictedMoments, notMoved, fromTempDraw);

    EXPECT_EQ(0u, notMoved.stock.size()); // the temp draw is not kept

    auto const& a = fromStock.rho;
    auto const& b = fromTempDraw.rho;
    for (std::size_t i = 0; i < a.size(); ++i)
        EXPECT_DOUBLE_EQ(a.data()[i], b.data()[i]);
    for (std::size_t c = 0; c < 3; ++c)
        for (std::size_t i = 0; i < a.size(); ++i)
            EXPECT_DOUBLE_EQ(fromStock.F[c].data()[i], fromTempDraw.F[c].data()[i]);
}


TEST_F(InflowMoments2D, depositsOnlyTheCellsInFrontOfThePatch)
{
    for (auto const stage : {ParticleBoundaryStage::Moments, ParticleBoundaryStage::PredictedMoments})
    {
        Arrays arrays{particleBox};
        Moments moments{layout};
        apply(ParticleBoundaryStage::Init, arrays, moments);
        apply(stage, arrays, moments);

        // the stock covers the whole ghost box (nCellsY2D + 2 cells)...
        EXPECT_EQ((nCellsY2D + 2) * ppc, arrays.stock.size());

        // ...but only the nCellsY2D cells in front of the patch are deposited. The weights of a
        // particle sum to 1 over the nodes, so the sum of rho is the total weight deposited.
        double total = 0.;
        for (std::size_t i = 0; i < moments.rho.size(); ++i)
            total += moments.rho.data()[i];
        EXPECT_NEAR(density * nCellsY2D, total, 1e-10);
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
