"""Build a symmetric 55-degree swept wing in millimetres with analytic wrap guides."""
from math import comb, radians, tan
from pathlib import Path

import cadquery as cq
import numpy as np
from numpy.polynomial import Polynomial
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.Geom import Geom_BezierSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.gp import gp_Pnt

HERE = Path(__file__).resolve().parent
JOIN = 390.
TAPER = 115./JOIN
SWEEP = tan(radians(55))
Q = Polynomial([0., 1.])
U = 3*Q**2-2*Q**3
L = 240*Q*(1-Q)
N = 14
M = 8


def thickness(q):
    u = q*q*(3-2*q)
    return .45*(.2969*q*np.sqrt(3-2*q)-.1260*u-.3516*u**2+.2843*u**3-.1036*u**4)


def bernstein(power, n):
    """Elevate a power-basis polynomial to Bernstein degree n."""
    return np.array([
        sum(a*comb(i,k)/comb(n,k) for k,a in enumerate(power.coef) if k <= i)
        for i in range(n+1)
    ])


def basis(t, n):
    return np.array([comb(n,k)*t**k*(1-t)**(n-k) for k in range(n+1)]).T


def times_poly(poles, power):
    """Multiply a Bernstein curve by a power-basis polynomial."""
    m = len(poles)-1
    n = power.degree()
    other = bernstein(power,n)
    result = np.zeros(m+n+1)
    for j, p in enumerate(poles):
        for k, a in enumerate(other):
            result[j+k] += p*a*comb(m,j)*comb(n,k)/comb(m+n,j+k)
    return result


samples = (1-np.cos(np.linspace(0,np.pi,401)))/2
fit = basis(samples,N)
T = np.zeros(N+1)
T[1] = .45*.2969*np.sqrt(3)/N
T[2:-2] = np.linalg.lstsq(
    fit[:,2:-2], thickness(samples)-fit@T, rcond=None,
)[0]
assert np.max(abs(fit@T-thickness(samples))) < 1e-9


def make_face(side, skin, cap):
    """Build a tensor-product patch, exact in the wrap parameter."""
    if cap:
        # f turns the span coordinate vertical through the lip; g starts with
        # zero first and second derivatives to match the unchanged skin.
        f = Polynomial([0, 1.6, 0, -.8, .2])
        g = Polynomial([1, 0, 0, -2, 1])
        base = JOIN
        extension = L
    else:
        f = Polynomial([0, 1])
        g = Polynomial([1])
        base = 0.
        extension = Polynomial([JOIN])

    chord = 210-TAPER*base
    n = N+extension.degree()
    net = np.zeros((n+1,M+1,3))
    net[:,:,0] = (
        bernstein(SWEEP*base+chord*U,n)[:,None]
        +np.outer(bernstein(extension*(SWEEP-TAPER*U),n),bernstein(f,M))
    )
    net[:,:,1] = side*(
        base+np.outer(bernstein(extension,n),bernstein(f,M))
    )
    net[:,:,2] = skin*(
        np.outer(times_poly(T,Polynomial([chord]+[0.]*extension.degree())),bernstein(g,M))
        -TAPER*np.outer(times_poly(T,extension),bernstein(f*g,M))
    )
    poles = TColgp_Array2OfPnt(1,n+1,1,M+1)
    for i in range(n+1):
        for j in range(M+1):
            poles.SetValue(i+1,j+1,gp_Pnt(*net[i,j]))
    return BRepBuilderAPI_MakeFace(Geom_BezierSurface(poles),1e-8).Face()


def build():
    sewing = BRepBuilderAPI_Sewing(1e-7)
    for side in (-1,1):
        for skin in (-1,1):
            for cap in (False,True):
                sewing.Add(make_face(side,skin,cap))
    sewing.Perform()
    shell = cq.Shape.cast(sewing.SewedShape())
    assert shell.ShapeType() == "Shell" and sewing.NbFreeEdges() == 0, (
        shell.ShapeType(),sewing.NbFreeEdges()
    )
    shell.wrapped.Closed(True)
    solid = cq.Solid.makeSolid(shell)
    assert solid.isValid() and solid.Volume() > 0
    cq.exporters.export(solid, str(HERE/"wing.step"))


if __name__ == "__main__":
    build()
