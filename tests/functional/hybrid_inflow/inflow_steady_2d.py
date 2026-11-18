#!/usr/bin/env python3
"""
Hybrid inflow, steady state. 2D, root level only (no AMR), run with MPI.

The domain holds a uniform plasma flowing along B: n = 1, V = (V0, 0, 0), B = (1, 0, 0).
V is parallel to B, so E = -V x B = 0 and the uniform flow is an exact steady state.
The x-lower face injects the same plasma (hybrid-inflow); the x-upper face is open.
y is periodic. Patches of 16 cells: two patches touch the inflow face, so the seam between
them is crossed by the stock particles too.

        y ┤ periodic
          │ ─► n = 1, V0 ─►  ─►  ─►  ─►
 inflow   │                              open
          │ ─►  ─►  ─►  ─►  ─►  ─►  ─►
          └──────────────────────────── x
            checked: the 6 nodes next to the face, and the mean over the left half

If the inflow lets in the right plasma, the density and the bulk velocity stay at their
initial values next to the face, for the whole run. Too few particles (wrong flux) or too many
(a stock particle counted twice at a seam) show as a density step at the face.

The x-upper side is not checked: on this branch the particle part of the open rule is not
written yet, so the density drops there. The left half is far enough from it for the run time.

Usage:  mpirun -n 4 python3 inflow_steady_2d.py
"""

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

V0 = 0.5
VTH = 0.3
CELLS = (40, 32)
DL = (0.2, 0.2)
PPC = 100
DIAG_DIR = "phare_outputs/hybrid_inflow_steady_2d"
TIME_STEP = 0.005
FINAL_TIME = 2.0
DUMP_EVERY = 40  # steps

NODES_CHECKED = 6

# noise of a mean over y: about 1 / sqrt(PPC * CELLS[1]) = 0.018 on n, VTH times that on V.
# The bounds are about 4 times the noise (a run on kaa, 2026-11-17: |n - 1| <= 0.025).
N_TOLERANCE = 0.08
V_TOLERANCE = 0.03
MEAN_N_TOLERANCE = 0.02  # mean over 20 x nodes: about 4 times less noise


def config():
    sim = ph.Simulation(
        smallest_patch_size=16,
        largest_patch_size=16,
        time_step=TIME_STEP,
        final_time=FINAL_TIME,
        cells=CELLS,
        dl=DL,
        interp_order=1,
        hyper_resistivity=0.01,
        domain_boundaries={
            "xlower": {
                "type": "hybrid-inflow",
                "velocity": [V0, 0.0, 0.0],
                "B": [1.0, 0.0, 0.0],
                "populations": {"protons": {"density": 1.0, "vth": [VTH, VTH, VTH]}},
            },
            "xupper": {"type": "open"},
        },
        diag_options={
            "format": "phareh5",
            "options": {"dir": DIAG_DIR, "mode": "overwrite"},
        },
    )

    def one(x, y):
        return 1.0

    def zero(x, y):
        return 0.0

    def vx(x, y):
        return V0

    def vth(x, y):
        return VTH

    ph.MaxwellianFluidModel(
        bx=one,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": one,
            "vbulkx": vx,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": PPC,
            "init": {"seed": 12334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.1)

    n_dumps = int(round(FINAL_TIME / (TIME_STEP * DUMP_EVERY))) + 1
    timestamps = TIME_STEP * DUMP_EVERY * np.arange(n_dumps)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim, timestamps


def profile_over_y(hier, name):
    """Mean over y of a primal quantity, on each x node of the domain (level 0)."""
    sums, counts = {}, {}
    for patch in hier.level(0).patches:
        values = np.asarray(patch.patch_datas[name][patch.box])  # domain nodes, +1 per direction
        x0 = patch.box.lower[0]
        for ix in range(values.shape[0]):
            i = x0 + ix
            # the last y node is the first one of the next patch (or the periodic image)
            sums[i] = sums.get(i, 0.0) + float(np.sum(values[ix, :-1]))
            counts[i] = counts.get(i, 0) + values.shape[1] - 1
    nodes = np.array(sorted(sums))
    return np.array([sums[i] / counts[i] for i in nodes])


def check(timestamps):
    from pyphare.pharesee.run import Run

    run = Run(DIAG_DIR)
    print("\nmean over y, at the x nodes next to the inflow face (node 0 is on the face)")
    print("   t     n: " + "  ".join(f"x{i:<4d}" for i in range(NODES_CHECKED)) + "   Vx: worst")

    failures = []
    for t in timestamps:
        n = profile_over_y(run.GetNi(t), "value")
        vx = profile_over_y(run.GetVi(t), "x")
        near_n, near_vx = n[:NODES_CHECKED], vx[:NODES_CHECKED]
        left = n[: len(n) // 2]

        worst_vx = float(np.max(np.abs(near_vx - V0)))
        print(
            f"{t:6.3f}     "
            + "  ".join(f"{v:5.3f}" for v in near_n)
            + f"   {worst_vx:.3f}   mean n (left half) {left.mean():.4f}"
        )

        if not np.all(np.isfinite(n)) or not np.all(np.isfinite(vx)):
            failures.append(f"t={t}: non-finite n or Vx")
            continue
        if np.max(np.abs(near_n - 1.0)) > N_TOLERANCE:
            failures.append(f"t={t}: n next to the face = {near_n}, expected 1 +- {N_TOLERANCE}")
        if worst_vx > V_TOLERANCE:
            failures.append(f"t={t}: Vx next to the face = {near_vx}, expected {V0} +- {V_TOLERANCE}")
        if abs(left.mean() - 1.0) > MEAN_N_TOLERANCE:
            failures.append(f"t={t}: mean n over the left half = {left.mean()}, expected 1")

    assert not failures, "the inflow does not hold the steady state:\n" + "\n".join(failures)


def main():
    sim, timestamps = config()
    Simulator(sim).run()
    if cpp.mpi_rank() == 0:
        check(timestamps)


if __name__ == "__main__":
    startMPI()
    main()
