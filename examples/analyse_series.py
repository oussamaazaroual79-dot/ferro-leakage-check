"""Analyse your own frequency series.
Input: one CSV file per frequency with two columns, V (volt) and P (uC/cm^2), one period of a sinusoidal
excitation sampled at N points starting at V = 0 on the rising side.
Usage: python analyse_series.py 100:loop_100Hz.csv 500:loop_500Hz.csv ... --area 1.72e-8 --thickness 650e-9
"""
import sys, argparse, numpy as np
sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
from pipeline import *
ap = argparse.ArgumentParser(); ap.add_argument("files", nargs="+"); ap.add_argument("--area", type=float, required=True)
ap.add_argument("--thickness", type=float, required=True); a = ap.parse_args()
f, Pr, Vc, W, L = [], [], [], [], {}
for item in a.files:
    fr, path = item.split(":", 1); fr = float(fr); x, y = np.loadtxt(path, delimiter=",", unpack=True)
    p, v = pr_vc(x, y); f.append(fr); Pr.append(p); Vc.append(v); W.append(loop_area(x, y)); L[fr] = (x, y)
f = np.array(f); lf = leak_fit(f, Pr, W); kP = lf["P"]["k"]; V0 = max(L[f.max()][0])
ratio = lf["W"]["k"] / (np.pi * V0 * kP)
R_, rho = resistivity(kP, V0, a.area, a.thickness)
Vcc = [pr_vc(L[fr][0], compensate(*L[fr], kP, fr))[1] for fr in f]
print(f"P_r(inf) = {lf['P']['inf']:.2f} +- {lf['P']['s_inf']:.2f} uC/cm2;  k_P = {kP:.0f} +- {lf['P']['s_k']:.0f} uC cm-2 Hz")
print(f"k_W/(pi V0 k_P) = {ratio:.2f}  (ohmic leakage ~0.9-1.1; slow switching << 1, see synthetic.py)")
print(f"R = {R_/1e6:.1f} MOhm, rho = {rho:.3g} Ohm m;  beta (compensated V_c) = {powerlaw_exponent(f, np.array(Vcc)):.3f}")
