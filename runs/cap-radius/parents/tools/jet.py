# ---------------------------------------------------------------------------
# Second-order forward differentiation (NumPy only).
#
# A Jet carries a value together with its exact first and second partial
# derivatives with respect to n seed variables. Derivative axes are trailing:
# d has shape value.shape + (n,), h has shape value.shape + (n, n).
# ---------------------------------------------------------------------------


class Jet:
    __array_priority__ = 1000
    __slots__ = ("v", "d", "h")

    def __init__(self, v, d, h):
        self.v, self.d, self.h = v, d, h

    @property
    def n(self):
        return self.d.shape[-1]

    @classmethod
    def variable(cls, value, index, n):
        v = np.array(value, dtype=float)
        d = np.zeros(v.shape + (n,))
        d[..., index] = 1.0
        return cls(v, d, np.zeros(v.shape + (n, n)))

    @classmethod
    def constant(cls, value, n):
        v = np.array(value, dtype=float)
        return cls(v, np.zeros(v.shape + (n,)), np.zeros(v.shape + (n, n)))

    def _lift(self, other):
        return other if isinstance(other, Jet) else Jet.constant(other, self.n)

    def chain(self, f0, f1, f2):
        """Compose with a scalar function given its value and first two derivatives."""
        f1 = np.asarray(f1, dtype=float)
        f2 = np.asarray(f2, dtype=float)
        d = self.d
        return Jet(np.asarray(f0, dtype=float), f1[..., None] * d,
                   f1[..., None, None] * self.h + f2[..., None, None] * d[..., :, None] * d[..., None, :])

    def __add__(self, other):
        other = self._lift(other)
        return Jet(self.v + other.v, self.d + other.d, self.h + other.h)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d, -self.h)

    def __sub__(self, other):
        return self + (-self._lift(other))

    def __rsub__(self, other):
        return self._lift(other) - self

    def __mul__(self, other):
        if not isinstance(other, Jet):
            c = np.asarray(other, dtype=float)
            return Jet(self.v * c, self.d * c[..., None], self.h * c[..., None, None])
        a, b = self, other
        return Jet(
            a.v * b.v,
            a.d * b.v[..., None] + b.d * a.v[..., None],
            a.h * b.v[..., None, None] + b.h * a.v[..., None, None]
            + a.d[..., :, None] * b.d[..., None, :] + b.d[..., :, None] * a.d[..., None, :],
        )

    __rmul__ = __mul__

    def reciprocal(self):
        v = self.v
        return self.chain(1.0 / v, -1.0 / v**2, 2.0 / v**3)

    def __truediv__(self, other):
        if not isinstance(other, Jet):
            return self * (1.0 / np.asarray(other, dtype=float))
        return self * other.reciprocal()

    def __rtruediv__(self, other):
        return self.reciprocal() * other

    def __pow__(self, p):
        if not isinstance(p, int) or p < 0:
            raise ValueError("Jet powers take non-negative integers; use jsqrt for roots")
        v = self.v
        f1 = p * v ** (p - 1) if p >= 1 else np.zeros_like(v)
        f2 = p * (p - 1) * v ** (p - 2) if p >= 2 else np.zeros_like(v)
        return self.chain(v**p, f1, f2)


def jsin(x):
    return x.chain(np.sin(x.v), np.cos(x.v), -np.sin(x.v))


def jcos(x):
    return x.chain(np.cos(x.v), -np.sin(x.v), -np.cos(x.v))


def jsqrt(x):
    r = np.sqrt(x.v)
    return x.chain(r, 0.5 / r, -0.25 / (r * x.v))


def jarccos(x):
    w = np.sqrt(1.0 - x.v**2)
    return x.chain(np.arccos(x.v), -1.0 / w, -x.v / w**3)


def jabs(x):
    """|x| with the derivative of the positive branch at x = 0 (one-sided at a mirror plane)."""
    sign = np.where(x.v < 0, -1.0, 1.0)
    return x * sign


def jwhere(condition, a, b):
    condition = np.asarray(condition)
    return Jet(np.where(condition, a.v, b.v), np.where(condition[..., None], a.d, b.d),
               np.where(condition[..., None, None], a.h, b.h))


def _seed(chord, span):
    """Seed the two surface variables: index 0 is the chord parameter, index 1 the span station."""
    chord, span = np.broadcast_arrays(np.asarray(chord, float), np.asarray(span, float))
    return Jet.variable(chord, 0, 2), Jet.variable(span, 1, 2)


def _const(value, like):
    return Jet.constant(np.broadcast_to(np.asarray(value, float), like.v.shape), like.n)


def _pack(xyz, chord_name="chord", span_name="span"):
    """Turn three coordinate jets into point and derivative arrays (last axis = XYZ)."""
    shape = np.broadcast_shapes(*(c.v.shape for c in xyz))
    xyz = [Jet(np.broadcast_to(c.v, shape), np.broadcast_to(c.d, shape + (2,)),
               np.broadcast_to(c.h, shape + (2, 2))) for c in xyz]
    stack = lambda f: np.stack([f(c) for c in xyz], axis=-1)
    return {
        "point": stack(lambda c: c.v),
        f"d_{chord_name}": stack(lambda c: c.d[..., 0]),
        f"d_{span_name}": stack(lambda c: c.d[..., 1]),
        f"d_{chord_name}_{chord_name}": stack(lambda c: c.h[..., 0, 0]),
        f"d_{chord_name}_{span_name}": stack(lambda c: c.h[..., 0, 1]),
        f"d_{span_name}_{span_name}": stack(lambda c: c.h[..., 1, 1]),
    }


def _pack1(xyz):
    """Point and derivatives of a curve seeded with one variable p."""
    shape = np.broadcast_shapes(*(c.v.shape for c in xyz))
    stack = lambda f: np.stack([np.broadcast_to(f(c), shape) for c in xyz], axis=-1)
    return {"point": stack(lambda c: c.v), "d_p": stack(lambda c: c.d[..., 0]),
            "d_p_p": stack(lambda c: c.h[..., 0, 0])}


def _bisect(function, lo, hi, iterations=200):
    """Vectorised bisection for a sign change of function between lo and hi."""
    lo = np.array(lo, dtype=float)
    hi = np.array(hi, dtype=float)
    lo, hi = np.broadcast_arrays(lo, hi)
    lo, hi = lo.copy(), hi.copy()
    f_lo = function(lo)
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        f_mid = function(mid)
        left = np.sign(f_mid) == np.sign(f_lo)
        lo = np.where(left, mid, lo)
        f_lo = np.where(left, f_mid, f_lo)
        hi = np.where(left, hi, mid)
        if np.all(hi - lo <= 4e-16 * np.maximum(1.0, np.abs(hi))):
            break
    return 0.5 * (lo + hi)


def _implicit(residual, roots, outer):
    """Exact first and second derivatives of implicitly defined unknowns.

    residual(unknowns, outer) returns k residual jets for k unknown jets. roots
    are the solved values. The unknowns depend on the outer jet only; the
    result is expressed in the outer jet's seed variables (implicit function
    theorem through second order).
    """
    k = len(roots)
    shape = np.broadcast_shapes(outer.v.shape, *(np.shape(r) for r in roots))
    unknowns = [Jet.variable(np.broadcast_to(r, shape), i, k + 1) for i, r in enumerate(roots)]
    seed = Jet.variable(np.broadcast_to(outer.v, shape), k, k + 1)
    res = residual(unknowns, seed)
    jac = np.stack([np.stack([r.d[..., j] for j in range(k)], -1) for r in res], -2)
    r_y = np.stack([r.d[..., k] for r in res], -1)
    r_uu = np.stack([r.h[..., :k, :k] for r in res], -3)
    r_uy = np.stack([r.h[..., :k, k] for r in res], -2)
    r_yy = np.stack([r.h[..., k, k] for r in res], -1)
    first = -np.linalg.solve(jac, r_y[..., None])[..., 0]
    rhs = r_yy + 2 * np.einsum("...ij,...j->...i", r_uy, first) \
        + np.einsum("...ijl,...j,...l->...i", r_uu, first, first)
    second = -np.linalg.solve(jac, rhs[..., None])[..., 0]
    od = np.broadcast_to(outer.d, shape + (outer.n,))
    oh = np.broadcast_to(outer.h, shape + (outer.n, outer.n))
    return [Jet(np.broadcast_to(np.asarray(r, float), shape).copy(),
                first[..., i, None] * od,
                first[..., i, None, None] * oh + second[..., i, None, None] * od[..., :, None] * od[..., None, :])
            for i, r in enumerate(roots)]


def _newton_polish(residual, roots, outer_value, steps=4):
    """Polish real roots of residual(unknowns, outer) with Newton steps (k unknowns)."""
    k = len(roots)
    roots = [np.array(r, dtype=float) for r in roots]
    for _ in range(steps):
        shape = np.broadcast_shapes(np.shape(outer_value), *(np.shape(r) for r in roots))
        unknowns = [Jet.variable(np.broadcast_to(r, shape), i, k) for i, r in enumerate(roots)]
        outer = Jet.constant(np.broadcast_to(outer_value, shape), k)
        res = residual(unknowns, outer)
        jac = np.stack([np.stack([r.d[..., j] for j in range(k)], -1) for r in res], -2)
        value = np.stack([r.v for r in res], -1)
        step = np.linalg.solve(jac, value[..., None])[..., 0]
        roots = [np.broadcast_to(r, shape) - step[..., i] for i, r in enumerate(roots)]
    return roots
