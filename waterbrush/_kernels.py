"""Hot loops. Compiled with Numba when it is installed; plain Python otherwise.

Loops run in a fixed order with no parallel reductions, so results are
bit-for-bit repeatable.
"""
import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def thin_film_flow(w, G, M, z, mob, dx, dt):
    """One explicit step of thin-film flow with upwind pigment transport.

    w: water depth (H, W); G: suspended pigment (n, H, W); M: wet mask;
    z: surface elevation of the paper; mob: mobility coefficient so that the
    flux across a face is -mob * w_face^3 * d(w + z)/dx. Updates w and G in place.
    """
    H, W = w.shape
    n = G.shape[0]
    qx = np.zeros((H, W + 1))
    qy = np.zeros((H + 1, W))
    out = np.zeros((H, W))
    for i in range(H):
        for j in range(W - 1):
            if M[i, j] and M[i, j + 1]:
                wf = 0.5 * (w[i, j] + w[i, j + 1])
                q = -mob * wf * wf * wf * ((w[i, j + 1] + z[i, j + 1]) - (w[i, j] + z[i, j])) / dx
                qx[i, j + 1] = q
                if q > 0:
                    out[i, j] += q
                else:
                    out[i, j + 1] -= q
    for i in range(H - 1):
        for j in range(W):
            if M[i, j] and M[i + 1, j]:
                wf = 0.5 * (w[i, j] + w[i + 1, j])
                q = -mob * wf * wf * wf * ((w[i + 1, j] + z[i + 1, j]) - (w[i, j] + z[i, j])) / dx
                qy[i + 1, j] = q
                if q > 0:
                    out[i, j] += q
                else:
                    out[i + 1, j] -= q
    # Positivity limiter: a cell may not send out more water than it holds.
    k = dt / dx
    lim = np.ones((H, W))
    for i in range(H):
        for j in range(W):
            o = out[i, j] * k
            if o > w[i, j]:
                lim[i, j] = w[i, j] / o
    dw = np.zeros((H, W))
    dG = np.zeros((n, H, W))
    for i in range(H):
        for j in range(W - 1):
            q = qx[i, j + 1]
            if q == 0.0:
                continue
            src = j if q > 0 else j + 1
            q *= lim[i, src]
            dw[i, j] -= q
            dw[i, j + 1] += q
            ws = w[i, src]
            if ws > 1e-12:
                for p in range(n):
                    f = q * G[p, i, src] / ws
                    dG[p, i, j] -= f
                    dG[p, i, j + 1] += f
    for i in range(H - 1):
        for j in range(W):
            q = qy[i + 1, j]
            if q == 0.0:
                continue
            src = i if q > 0 else i + 1
            q *= lim[src, j]
            dw[i, j] -= q
            dw[i + 1, j] += q
            ws = w[src, j]
            if ws > 1e-12:
                for p in range(n):
                    f = q * G[p, src, j] / ws
                    dG[p, i, j] -= f
                    dG[p, i + 1, j] += f
    for i in range(H):
        for j in range(W):
            v = w[i, j] + k * dw[i, j]
            w[i, j] = v if v > 0.0 else 0.0
            for p in range(n):
                g = G[p, i, j] + k * dG[p, i, j]
                G[p, i, j] = g if g > 0.0 else 0.0
