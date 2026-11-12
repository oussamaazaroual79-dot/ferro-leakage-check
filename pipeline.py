"""Analysis pipeline used in the article (functions shared by the real-data and synthetic-data scripts).

Conventions: voltage in V, polarisation in uC/cm^2, current in A, time in s, frequency in Hz.
Loops are one period of a sinusoidal excitation sampled at N points, starting at V = 0 on the rising side.
"""
import numpy as np

EPS0 = 8.854e-12


def crossings(x, y):
    """Zero crossings of a closed loop: values of y at x = 0 (P_r) and of x at y = 0 (V_c)."""
    P, V = [], []
    for i in range(len(x) - 1):
        if x[i] * x[i + 1] < 0 or x[i] == 0:
            P.append(y[i] if x[i] == 0 else y[i] + (y[i + 1] - y[i]) * (-x[i]) / (x[i + 1] - x[i]))
        if y[i] * y[i + 1] < 0 or y[i] == 0:
            V.append(x[i] if y[i] == 0 else x[i] + (x[i + 1] - x[i]) * (-y[i]) / (y[i + 1] - y[i]))
    return np.array(P), np.array(V)


def pr_vc(x, y):
    """Half-differences of the crossings: P_r and V_c of a loop."""
    P, V = crossings(np.r_[x, x[:1]], np.r_[y, y[:1]])
    return (P.max() - P.min()) / 2, (V.max() - V.min()) / 2


def loop_area(x, y):
    """Area of the loop (energy density per cycle in uC cm^-2 V)."""
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def current_analysis(t, V, I, Pmax, vthr=13.5):
    """Linear capacitance from the high-field current and electrode area from the charge of one rising half-cycle.
    Pmax: half peak-to-peak polarisation of the loop measured simultaneously (uC/cm^2).
    Returns C_l (F), I0 (A), A (m^2)."""
    dVdt = np.gradient(V, t)
    m = np.abs(V) > vthr
    M = np.vstack([dVdt[m], np.ones(m.sum())]).T
    (C, I0), *_ = np.linalg.lstsq(M, I[m], rcond=None)
    closed = abs(V[-1] - V[0]) < 1e-3 * np.abs(V).max()     # last sample repeats the first one
    T = t[-1] - t[0] if closed else t[-1] - t[0] + (t[1] - t[0])
    tr = np.where(t >= 0.75 * T, t - T, t)
    rise = (t - t[0] <= 0.25 * T) | (t - t[0] >= 0.75 * T)
    o = np.argsort(tr[rise])
    Q = np.trapezoid(I[rise][o], tr[rise][o])
    A = Q / (2 * Pmax * 1e-2)
    return C, I0, A


def leak_fit(f, Pr, W):
    """Fit P_r = P_inf + k_P/f and W = W_inf + k_W/f. Returns dict with values and standard errors."""
    f = np.asarray(f, float)
    M = np.vstack([np.ones_like(f), 1 / f]).T
    out = {}
    for name, y in (("P", np.asarray(Pr)), ("W", np.asarray(W))):
        c, *_ = np.linalg.lstsq(M, y, rcond=None)
        r = y - M @ c
        dof = max(len(f) - 2, 1)
        s2 = (r @ r) / dof
        cov = s2 * np.linalg.inv(M.T @ M)
        out[name] = dict(inf=c[0], k=c[1], s_inf=np.sqrt(cov[0, 0]), s_k=np.sqrt(cov[1, 1]), rms=np.sqrt(np.mean(r ** 2)))
    return out


def leak_shape(x, n=1.0):
    """Normalised leakage-charge waveform for a leakage current J ~ s|s|^(n-1), s = V/V0, under sinusoidal drive.
    Returns w(phi) such that the leakage contribution to P is (k_P/f) * w, with w = -cos(phi) for n = 1
    (so that the half-difference at V = 0 equals k_P/f)."""
    N = len(x)
    phi = 2 * np.pi * np.arange(N) / N
    fine = np.linspace(0, 2 * np.pi, 20 * N, endpoint=False)
    s = np.sin(fine)
    J = s * np.abs(s) ** (n - 1)
    q = np.cumsum(J) * (fine[1] - fine[0])
    q -= q.mean()
    q /= (q.max() - q.min()) / 2           # half-difference between V = 0 crossings is 1
    return np.interp(phi, fine, q)


def compensate(x, y, kP, f, n=1.0):
    """Subtract the leakage charge (k_P/f) * w(phi) from a loop sampled on the standard phase grid."""
    return y - (kP / f) * leak_shape(x, n)


def resistivity(kP, V0, A, t):
    """Leakage resistance (ohm) and resistivity (ohm m) from k_P (uC cm^-2 Hz)."""
    R = V0 / (2 * np.pi * kP * 1e-2 * A)
    return R, R * A / t


def powerlaw_exponent(f, y):
    p = np.polyfit(np.log(f), np.log(y), 1)
    return p[0]
