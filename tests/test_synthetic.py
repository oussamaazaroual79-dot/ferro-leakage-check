import sys, numpy as np
sys.path.insert(0, __file__.rsplit("/tests", 1)[0])
import synthetic as S

def test_recovers_ohmic_parameters():
    p = dict(S.BASE["PLT7"]); r = S.run_once(p, np.random.default_rng(0))
    assert abs(r["A"] / p["A"] - 1) < 0.02
    assert abs(r["kP"] / p["kP"] - 1) < 0.02
    assert abs(r["beta"] - p["beta"]) < 0.005

def test_slow_switching_fails_area_test():
    p = dict(S.BASE["PLT7"]); p["kP"] = 0.0
    r = S.run_once(p, np.random.default_rng(0), kS=650.)
    assert r["ratio"] < 0.6
