#ifndef PHARE_CORE_NUMERICS_HYBRID_INFLOW_HYBRID_INFLOW_WIND_HPP
#define PHARE_CORE_NUMERICS_HYBRID_INFLOW_HYBRID_INFLOW_WIND_HPP

#include "core/numerics/boundary_condition/particle_boundary_condition.hpp"

#include "initializer/data_provider.hpp"
#include "initializer/dict_utils.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

/** @file What a hybrid inflow face derives from its wind (the populations that enter), kept out
 * of BoundaryFactory so that the factory stays the same for MHD and hybrid faces. */

namespace PHARE::core
{
/**
 * @brief Charge-weighted velocity of an inflow wind: V_q = sum_s q_s n_s V_s / sum_s q_s n_s.
 *
 * Ohm's law takes the electron velocity Ve = (sum_s q_s F_s - J) / Ne (Electrons). In the uniform
 * wind J = 0, so Ve = V_q, and E = -V_q x B is what Ohm's law gives on the other side of the face.
 * With a single velocity for every population, V_q is that velocity.
 *
 * @param wind {population name: {density, velocity, vth}}, as written by pyphare.
 * @param populations the populations of the model, for their charge. A wind population that the
 * model does not have is skipped here (the particle registration refuses it).
 */
inline std::array<double, 3> windChargeVelocity(initializer::PHAREDict const& wind,
                                                std::vector<ParticleBoundaryInfo> const& populations)
{
    std::array<double, 3> qnv{};
    double qn = 0.;
    for (auto const& pop : populations)
        if (wind.contains(pop.name))
        {
            auto const& w   = wind[pop.name];
            double const n  = w["density"].template to<double>();
            auto const v    = initializer::parseDimXYZType<double, 3>(w, "velocity");
            double const qw = pop.charge * n;
            qn += qw;
            for (std::size_t d = 0; d < 3; ++d)
                qnv[d] += qw * v[d];
        }
    if (!(qn > 0.))
        throw std::runtime_error("BoundaryFactory: HybridInflow: the wind has no positive charge "
                                 "density (sum of charge x density), cannot compute its velocity.");
    for (auto& c : qnv)
        c /= qn;
    return qnv;
}
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_HYBRID_INFLOW_HYBRID_INFLOW_WIND_HPP
