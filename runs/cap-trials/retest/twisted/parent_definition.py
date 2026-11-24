import numpy as np
CUT, TIP = 440.0, 500.0
DIHEDRAL = np.tan(np.deg2rad(5))
WASHOUT = np.deg2rad(3)
ANGLE = np.pi / 2
U = np.linspace(0., 1., 121)


def airfoil(u):
    """Closed-TE NACA 2412. Smooth the two camber parabolas over 0.36–0.44c."""
    x = np.sin(np.pi * u / 2)**2
    forward = 0.02 / 0.4**2 * (0.8 * x - x*x)
    aft = 0.02 / 0.6**2 * (0.2 + 0.8*x - x*x)
    fraction = np.clip((x - 0.36) / 0.08, 0, 1)
    blend = fraction**3 * (10 - 15*fraction + 6*fraction**2)
    blend_derivative = 30*fraction**2 * (1-fraction)**2 / 0.08
    camber = forward + (aft-forward)*blend
    slope = (0.8 - 2*x) * (0.02/0.4**2 * (1-blend) + 0.02/0.6**2 * blend)
    slope += (aft-forward)*blend_derivative
    thickness = 0.6 * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x**3 - 0.1036*x**4)
    thickness = np.where((u == 0) | (u == 1), 0, thickness)
    camber = np.where((u == 0) | (u == 1), 0, camber)
    normal = np.stack((-slope, np.ones_like(slope)), axis=-1) / np.sqrt(1 + slope*slope)[..., None]
    return np.stack((x, camber), axis=-1), thickness[..., None]*normal


def section(y, u, thickness_factor):
    """Apply taper and true twist about the rising leading-edge reference line."""
    xy, normal_offset = airfoil(u)
    x, z = (xy + np.asarray(thickness_factor)[..., None]*normal_offset).T
    theta = WASHOUT * y/CUT
    chord = 180 - 80*y/CUT
    return np.stack((
        90*y/CUT + chord*(x*np.cos(theta) - z*np.sin(theta)),
        np.broadcast_to(y, np.shape(x)),
        y*DIHEDRAL + chord*(x*np.sin(theta) + z*np.cos(theta)),
    ), axis=-1)


def span_jet(y, u, thickness_factor):
    """The first and second span derivatives of the unchanged parent skin."""
    position = section(y, u, thickness_factor)
    chord = 180 - 80*y/CUT
    theta_rate = WASHOUT/CUT
    chord_rate = -80/CUT
    baseline = np.array([90/CUT, 1., DIHEDRAL])
    offset = position - np.asarray(y)[..., None]*baseline
    rotated = np.stack((-offset[:, 2], np.zeros(len(u)), offset[:, 0]), axis=-1)
    return (
        baseline + np.asarray(chord_rate/chord)[..., None]*offset + theta_rate*rotated,
        2*theta_rate*np.asarray(chord_rate/chord)[..., None]*rotated - theta_rate**2*offset,
    )


def reach(u):
    """Outboard distance to the common 3D perimeter for each chord guide."""
    return np.where((u == 0) | (u == 1), 0, (TIP-CUT)*np.sin(np.pi*u))
