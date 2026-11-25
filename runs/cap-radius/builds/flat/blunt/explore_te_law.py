"""Design-iteration scan of the trailing-edge reach law (analytic surface, before CAD).

Usage (project root): .venv/bin/python runs/cap-radius/builds/flat/blunt/explore_te_law.py 1e-3,3e-3 2e-3,6e-3
Each argument is EPS0,EPS1 for design.LAW = 'local'. Prints the smallest principal radius against the
printable rule and the lip's plan deviation from the locked outline.
"""

import sys, numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design as D
def grid():
    a = np.concatenate([np.geomspace(1e-6, 0.05, 120), np.linspace(0.05, 0.95, 400)[1:-1], 1-np.geomspace(0.05, 1e-9, 200)])
    s = np.unique(np.concatenate([[1.0], a]))
    w = np.unique(np.concatenate([np.linspace(-1,1,401), 1-np.geomspace(1e-6,0.05,60), -(1-np.geomspace(1e-6,0.05,60)), np.geomspace(1e-5, 0.05, 40), -np.geomspace(1e-5,0.05,40)]))
    w = w[w>=0]  # symmetric
    return s, w
import itertools
A = 0.0
for spec in sys.argv[1:]:
    e0, e1 = map(float, spec.split(','))
    D.LAW = 'local'; D.EPS0 = e0; D.EPS1 = e1
    print('LAW local EPS0', e0, 'EPS1', e1)
    s, w = grid(); S, W = np.meshgrid(s, w, indexing='ij')
    j = D.surface_jet(S, W); pc = D.principal(j)
    r = 1/pc['kmax']
    h = D.half_thickness(s)[:,None]*np.ones_like(W)
    req = np.where(2*h >= 0.8, 0.4, h)
    ratio = r/req
    i = np.unravel_index(np.nanargmin(ratio), ratio.shape)
    print('A', A, 'K', D.K, 'min ratio', ratio[i], 'r', r[i], 'at s', S[i], 'w', W[i], 'X', j['point'][i][0], 'h', h[i])
    # min radius in thick zone
    thick = 2*h >= 0.8
    rr = np.where(thick, r, np.inf); i2 = np.unravel_index(np.argmin(rr), rr.shape)
    print('  thick-zone min r', rr[i2], 's', S[i2], 'w', W[i2], 'P', j['point'][i2])
    rr = np.where(~thick, r, np.inf); i3 = np.unravel_index(np.argmin(rr), rr.shape)
    print('  thin-zone min r', rr[i3], 's', S[i3], 'w', W[i3], 'P', j['point'][i3], 'h', h[i3])
    # lip radius along s
    for sv in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99, 0.999, 0.9999, 1.0]:
        k = np.argmin(abs(s-sv)); lip = r[k, 0]; mn = r[k].min(); wm = w[np.argmin(r[k])]
        print(f'   s={s[k]:.6f} X={j["point"][k,0,0]:.5f} R={j["point"][k,0,1]-350:.4f} h={h[k,0]:.4f} lip r={lip:.4f} min r={mn:.4f} at w={wm:.4f}')

    # plan deviation of the lip from the locked outline (distance to the curve)
    import parent as PP
    from scipy.spatial import cKDTree
    pp = np.unique(np.concatenate([np.linspace(0,1,200001), 1-np.geomspace(1e-9,1e-2,20000), np.geomspace(1e-9,1e-2,20000)]))
    oc = PP.tip_outline(pp)[:, :2]
    tree = cKDTree(oc)
    sl = np.unique(np.concatenate([np.linspace(0,1,20001), 1-np.geomspace(1e-10,1e-1,5000)]))
    lip = D.surface_jet(sl, 0.0)['point'][:, :2]
    dd, ii = tree.query(lip)
    # refine with local segment distance
    k = np.argmax(dd)
    print('  plan deviation max (kdtree) %.3g at s=%.8f X=%.6f Y=%.6f' % (dd[k], sl[k], lip[k,0], lip[k,1]))
    R = lip[:,1]-350
    print('  R monotone check: R increasing then decreasing? min dR in last 1e-2:', np.diff(R[sl>0.99]).max())
