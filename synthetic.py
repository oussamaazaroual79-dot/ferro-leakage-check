"""Validation of the pipeline on synthetic loops and currents with known parameters.

Model of the measured polarisation (uC/cm^2) under V = V0 sin(phi):
  P = Ps(f) * tanh((V -/+ Vc(f))/w)   (rising / falling branch)
    + chi * V                          (linear dielectric part)
    + (kP/f) * w_n(phi)                (leakage charge, J ~ s|s|^(n-1))
  Vc(f) = Vc0 (f/1 kHz)^beta,  Ps(f) = Ps + kS/f (slow-switching scenario only)
The switching current at 2 kHz is i = A dP/dt.  Noise is added to P and i.
"""
import json, numpy as np
from pipeline import *

N = 400
FREQ0 = [100, 500, 1000, 1300, 2000.]
BASE = dict(PT=dict(V0=16.0, Ps=17.5, Vc0=8.1, w=2.0, chi=0.95, kP=650., beta=0.015, A=0.0172e-6, t=650e-9),
            PLT7=dict(V0=15.0, Ps=10.9, Vc0=5.6, w=1.6, chi=1.05, kP=860., beta=0.05, A=0.0172e-6, t=650e-9),
            PLT14=dict(V0=15.0, Ps=11.9, Vc0=5.7, w=1.6, chi=1.08, kP=940., beta=0.03, A=0.0172e-6, t=650e-9,
                       freq=[100, 400, 800, 1000, 2000.]))
SIG_P, SIG_I = 0.02, 1.6e-6          # noise levels estimated from the measured loops and current


def synth(p, f, n_leak=1.0, kS=0.0, fine=8):
    M = N * fine
    phi = 2 * np.pi * np.arange(M) / M
    V = p["V0"] * np.sin(phi)
    rising = np.cos(phi) >= 0
    Vc = p["Vc0"] * (f / 1000) ** p["beta"]
    Ps = p["Ps"] + kS / f
    Psw = np.where(rising, Ps * np.tanh((V - Vc) / p["w"]), Ps * np.tanh((V + Vc) / p["w"]))
    P = Psw + p["chi"] * V + (p["kP"] / f) * leak_shape(V, n_leak) if p["kP"] else Psw + p["chi"] * V
    return phi, V, P


def run_once(p, rng, n_leak=1.0, kS=0.0, noise=1.0, n_comp=1.0):
    FREQ = np.array(p.get("freq", FREQ0), float)
    Pr, W, Vc_raw, loops = [], [], [], {}
    for f in FREQ:
        phi, V, P = synth(p, f, n_leak, kS)
        x, y = V[::8], P[::8] + rng.normal(0, SIG_P * noise, N)
        a, b = pr_vc(x, y); Pr.append(a); Vc_raw.append(b); W.append(loop_area(x, y)); loops[f] = (x, y)
    # switching current at 2 kHz
    phi, V, P = synth(p, 2000., n_leak, kS)
    T = 1 / 2000.
    tt = phi / (2 * np.pi) * T
    i = p["A"] * 1e-2 * np.gradient(P, tt)
    t, Vs, I = tt[::8], V[::8], i[::8] + rng.normal(0, SIG_I * noise, N)
    x2, y2 = loops[2000.]
    Pmax = (y2.max() - y2.min()) / 2
    C, I0, A = current_analysis(t, Vs, I, Pmax)
    lf = leak_fit(FREQ, Pr, W)
    kP = lf["P"]["k"]
    Vc_c = []
    for f in FREQ:
        x, y = loops[f]
        Vc_c.append(pr_vc(x, compensate(x, y, kP, f, n_comp))[1])
    beta = powerlaw_exponent(FREQ, np.array(Vc_c))
    R, rho = resistivity(kP, p["V0"], A, p["t"]) if kP > 0 else (np.inf, np.inf)
    return dict(A=A, C=C, Pinf=lf["P"]["inf"], kP=kP, ratio=lf["W"]["k"] / (np.pi * p["V0"] * kP),
                beta=beta, beta_raw=powerlaw_exponent(FREQ, np.array(Vc_raw)), rho=rho, rmsP=lf["P"]["rms"])


def truth(p):
    Pr = p["Ps"] * np.tanh(p["Vc0"] / p["w"])
    C = p["chi"] * 1e-2 * p["A"]
    R, rho = resistivity(p["kP"], p["V0"], p["A"], p["t"]) if p["kP"] else (np.inf, np.inf)
    return dict(A=p["A"], C=C, Pinf=Pr, kP=p["kP"], ratio=1.0, beta=p["beta"], rho=rho)


SCEN = {
    "S1 ohmic leakage": dict(),
    "S2 ohmic, noise x10": dict(noise=10.0),
    "S3 superlinear leakage (n = 3)": dict(n_leak=3.0),
    "S4 no leakage": dict(kP0=True),
    "S5 slow switching, no leakage": dict(kP0=True, kS=650.),
    "S6 ohmic leakage, beta = 0": dict(beta0=True),
}

if __name__ == "__main__":
    rng = np.random.default_rng(2013)
    import warnings; warnings.filterwarnings("ignore")
    res = {}
    for mat, p0 in BASE.items():
        res[mat] = {}
        for sname, opt in SCEN.items():
            p = dict(p0)
            if opt.get("kP0"): p["kP"] = 0.0
            if "kS" in opt: p["kS"] = opt["kS"]
            if opt.get("beta0"): p["beta"] = 0.0
            tr = truth(p)
            runs = [run_once(p, rng, opt.get("n_leak", 1.0), opt.get("kS", 0.0), opt.get("noise", 1.0)) for _ in range(200)]
            keys = runs[0].keys()
            stat = {k: (float(np.mean([r[k] for r in runs])), float(np.std([r[k] for r in runs]))) for k in keys}
            # for S3 also compensate with the right shape
            if sname.startswith("S3"):
                runs3 = [run_once(p, rng, 3.0, 0.0, 1.0, n_comp=3.0) for _ in range(100)]
                stat["beta_rightshape"] = (float(np.mean([r["beta"] for r in runs3])), float(np.std([r["beta"] for r in runs3])))
            res[mat][sname] = dict(truth={k: float(v) for k, v in tr.items()}, est=stat)
            e = stat
            print(f"{mat:5s} {sname:32s} A {e['A'][0]/tr['A']:.3f}±{e['A'][1]/tr['A']:.3f}  C {e['C'][0]/tr['C']:.3f}±{e['C'][1]/tr['C']:.3f}  "
                  f"Pinf {e['Pinf'][0]:.2f}±{e['Pinf'][1]:.2f} (true {tr['Pinf']:.2f})  kP {e['kP'][0]:.0f}±{e['kP'][1]:.0f} (true {tr['kP']:.0f})  "
                  f"ratio {e['ratio'][0]:.2f}±{e['ratio'][1]:.2f}  beta {e['beta'][0]:.3f}±{e['beta'][1]:.3f} (true {tr['beta']:.3f}; raw {e['beta_raw'][0]:.3f})  rmsP {e['rmsP'][0]:.3f}"
                  + (f"  beta(n=3 comp) {stat['beta_rightshape'][0]:.3f}" if 'beta_rightshape' in stat else ""))
    json.dump(res, open("synthetic.json", "w"), indent=1)
