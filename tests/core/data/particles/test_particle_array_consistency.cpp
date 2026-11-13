#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cmath>
#include <map>


namespace PHARE::core
{
std::size_t static constexpr cells = 3;
std::size_t static constexpr ppc   = 10;


template<std::size_t dim, typename ICell>
PHARE::core::Particle<dim> particle(ICell const& icell)
{
    return {/*.weight = */ 0,
            /*.iCell  = */ icell,
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
            /*.v      = */ {{.00001, .00001, .00001}}};
}

template<typename ParticleArray_t, typename Box_t>
void add_particles_in(ParticleArray_t& particles, Box_t const& box)
{
    for_each_point(box, [&](auto const& amr_idx) {
        for (std::size_t i = 0; i < ppc; ++i)
            particles.emplace_back(particle<ParticleArray_t::dimension>(*amr_idx));
    });
}



template<typename ParticleArray_>
struct ParticleArrayConsistencyTest : public ::testing::Test
{
    auto constexpr static dim    = ParticleArray_::dimension;
    auto constexpr static interp = 1;

    using GridLayout_t
        = TestGridLayout<typename PHARE_Types<SimOpts{dim, interp}>::Hybrid::GridLayout_t>;
    using ParticleArray_t = ParticleArray_;

    GridLayout_t layout{cells};
};



using Permutations_t = testing::Types<ParticleArray<1>>;


TYPED_TEST_SUITE(ParticleArrayConsistencyTest, Permutations_t, );



TYPED_TEST(ParticleArrayConsistencyTest, test_is_consistent_after_swap_copy)
{
    using ParticleArray_t     = TestFixture::ParticleArray_t;
    auto static constexpr dim = ParticleArray_t::dimension;

    auto levelGhostParticles = ParticleArray<dim>{this->layout.AMRBox()};
    add_particles_in(levelGhostParticles, this->layout.AMRBox());

    auto levelGhostParticlesNew = ParticleArray<dim>{this->layout.AMRBox()};
    add_particles_in(levelGhostParticlesNew, this->layout.AMRBox());

    auto levelGhostParticlesOld = ParticleArray<dim>{this->layout.AMRBox()};
    add_particles_in(levelGhostParticlesOld, this->layout.AMRBox());

    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
    levelGhostParticlesNew.clear();
    levelGhostParticles = levelGhostParticlesOld;

    EXPECT_EQ(levelGhostParticlesNew.size(), 0);
    EXPECT_EQ(levelGhostParticlesOld.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(levelGhostParticles.size(), ppc * std::pow(cells, TestFixture::dim));

    EXPECT_TRUE(levelGhostParticlesNew.is_consistent());
    EXPECT_TRUE(levelGhostParticlesOld.is_consistent());
    EXPECT_TRUE(levelGhostParticles.is_consistent());
}


// --- erase_in -------------------------------------------------------------------------------

// a 2D array over [0, 5] x [0, 5]. The particles are added cell after cell, ppc times over, so
// the particles of a cell are spread through the array, not next to each other.
ParticleArray<2> interleavedParticles()
{
    Box<int, 2> const box{Point{0, 0}, Point{5, 5}};
    ParticleArray<2> particles{box};
    for (std::size_t i = 0; i < ppc; ++i)
        for_each_point(box, [&](auto const& cell) {
            particles.emplace_back(particle<2>(std::array<int, 2>{cell[0], cell[1]}));
        });
    return particles;
}

// number of particles in each cell of the array's box, by cell
std::map<std::array<int, 2>, std::size_t> countPerCell(ParticleArray<2> const& particles)
{
    std::map<std::array<int, 2>, std::size_t> counts;
    for (auto const& particle : particles)
        ++counts[particle.iCell];
    return counts;
}

TEST(ParticleArrayEraseIn, erasesOnlyTheParticlesInTheBoxAndKeepsTheMapConsistent)
{
    auto particles    = interleavedParticles();
    auto const before = countPerCell(particles);
    Box<int, 2> const box{Point{0, 0}, Point{1, 5}}; // 2 columns of 6 cells

    particles.erase_in(box);

    EXPECT_EQ(particles.size(), (36 - 12) * ppc);
    EXPECT_EQ(particles.nbr_particles_in(box), 0u);
    EXPECT_TRUE(particles.is_consistent());

    auto const after = countPerCell(particles);
    for (auto const& [cell, count] : before)
    {
        if (isIn(Point{cell}, box))
            EXPECT_EQ(after.count(cell), 0u);
        else
            EXPECT_EQ(after.at(cell), count); // untouched
    }
}

TEST(ParticleArrayEraseIn, worksWithABoxPartlyOutsideTheArray)
{
    auto particles = interleavedParticles();
    particles.erase_in(Box<int, 2>{Point{-3, 4}, Point{2, 9}}); // overlaps [0, 2] x [4, 5]

    EXPECT_EQ(particles.size(), (36 - 6) * ppc);
    EXPECT_TRUE(particles.is_consistent());
}

TEST(ParticleArrayEraseIn, doesNothingWithABoxOutsideTheArray)
{
    auto particles = interleavedParticles();
    particles.erase_in(Box<int, 2>{Point{-3, -3}, Point{-1, -1}});

    EXPECT_EQ(particles.size(), 36 * ppc);
    EXPECT_TRUE(particles.is_consistent());
}

TEST(ParticleArrayEraseIn, canEmptyTheWholeArray)
{
    auto particles = interleavedParticles();
    particles.erase_in(particles.box());

    EXPECT_EQ(particles.size(), 0u);
    EXPECT_TRUE(particles.is_consistent());
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
