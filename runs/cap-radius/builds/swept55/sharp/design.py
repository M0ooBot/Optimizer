"""Choose the profile knob x0(q) so the lip meets the printable radius rule.

Rule (skill step 5), judged against the parent's attachment thickness t(q) at
the same chord position:
  t >= 0.8 mm : every radius >= 0.4 mm  -> design lip radius R_FLOOR = 0.5 mm (25 % margin)
  t <  0.8 mm : lip a full half-round of radius t/2
Design target R(q) = smooth-min(R_FLOOR, 1.03 t/2). For each q the lip's principal
radius is found as a function of x0 (bisection on the analytic jets), and x0
is capped where a larger value would put a tighter spot (shoulder) inboard of
the lip. x0(q) = smooth-max(0, x0_needed), slope-limited and smoothed along the chord
(soften), then held to the shoulder-free cap.
Writes design_table.json (q, x0, target, achieved lip radius, shoulder-free cap).
Run from the project root:
  .venv/bin/python runs/cap-radius/builds/swept55/sharp/design.py
"""
import json
import pathlib
import sys
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import cap_model as cm  # noqa: E402

warnings.simplefilter("ignore")
THIN_FACTOR = 1.03   # aim 3 % above t/2 so the 0.4 mm floor holds with margin where t is just over 0.8 mm
W_MAX = 0.03     # smoothing of max(0, x0)
TH = np.concatenate([np.linspace(1e-3, np.pi / 2 - 1e-3, 3000), [np.pi / 2]])


def smin(a, b, p=20):
    """Smooth min that stays proportional to b when b << a (needed as t -> 0 at the trailing edge)."""
    return (a ** -p + b ** -p) ** (-1.0 / p)


def smax(a, b, w):
    return 0.5 * (a + b + np.sqrt((a - b) ** 2 + w * w))


SLOPE = 3.0                      # max |dx0/dq| when x0 starts to rise (avoids a concave crease)
SIGMA = (0.03, 0.72, 0.84)       # Gaussian smoothing width in q, faded out between 0.72 and 0.84


def soften(q, x0, cap):
    """Slope-limited upper envelope, then Gaussian smoothing (only where x0 is not at the
    trailing-edge limit), then held under the shoulder-free cap."""
    env = np.max(x0[None, :] - SLOPE * np.abs(q[:, None] - q[None, :]), axis=1)
    s0, qa, qb = SIGMA
    sig = s0 * np.clip((qb - q) / (qb - qa), 0, 1)
    out = env.copy()
    for i, (qq, s) in enumerate(zip(q, sig)):
        if s > 1e-4:
            wgt = np.exp(-0.5 * ((q - qq) / s) ** 2)
            # mirror-free: renormalise over the available rows
            out[i] = (wgt * env).sum() / wgt.sum()
    return np.minimum(np.maximum(out, 0.0), cap)


def radii(q, x0):
    P = cm.principal(cm.surface_jet(np.full_like(TH, q), TH, x0_const=x0, b_const=float(cm.b_of_q(q)), ku_const=float(cm.ku_of_q(q))))
    r = 1 / P["kmax"]
    return r[-1], r.min()


def solve_q(q, target):
    lip0, _ = radii(q, 0.0)
    if lip0 >= target:
        x_need = target / lip0 - 1         # negative: how far below zero (for smooth max)
    else:
        lo, hi = 0.0, 0.97
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if radii(q, mid)[0] < target else (lo, mid)
        x_need = 0.5 * (lo + hi)
    # shoulder-free cap: largest x0 with min radius >= 0.999 lip radius
    lo, hi = 0.0, 0.97
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lip, mn = radii(q, mid)
        lo, hi = (mid, hi) if mn >= 0.999 * lip else (lo, mid)
    return x_need, lo


def main():
    q = np.unique(np.concatenate([np.linspace(0.02, 0.98, 97), 1 - np.geomspace(1e-3, 0.02, 8), np.geomspace(2e-3, 0.02, 8)]))
    t = cm.thickness_at_attachment(q)
    target = smin(cm.R_FLOOR, THIN_FACTOR * t / 2)
    need, cap = np.array([solve_q(qq, tt) for qq, tt in zip(q, target)]).T
    x0 = np.minimum(smax(0.0, need, W_MAX), cap)
    x0 = soften(q, x0, cap)
    lip = np.array([radii(qq, xx)[0] for qq, xx in zip(q, x0)])
    qn = np.concatenate([[0.0], q, [1.0]])
    # beyond the last row (q > 0.99) hold the last value: the pole limit of x0 (sizes there are microns)
    xn = np.concatenate([[x0[0]], x0, [x0[-1]]])
    out = {"q": qn.tolist(), "x0": xn.tolist(), "x0_needed": need.tolist(), "x0_shoulder_cap": cap.tolist(),
           "target_mm": target.tolist(), "lip_radius_mm": lip.tolist(), "q_rows": q.tolist(),
           "note": "x0(q) knob; target = design lip radius; lip_radius = analytic lip principal radius at that x0"}
    pathlib.Path(__file__).with_name("design_table.json").write_text(json.dumps(out, indent=1))
    for i in range(0, len(q), 4):
        print(f"q={q[i]:.4f} t={t[i]:.3f} target={target[i]:.3f} need={need[i]:.3f} cap={cap[i]:.3f} x0={x0[i]:.3f} lip={lip[i]:.4f}")


if __name__ == "__main__":
    main()
