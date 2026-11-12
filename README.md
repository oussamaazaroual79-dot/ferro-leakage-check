# ferro-leakage-check

Internal consistency tests for leakage in ferroelectric hysteresis loops measured at several frequencies,
as used in: O. Azaroual, *Separating switching, dielectric, and leakage contributions in sol–gel PbTiO3
and (Pb,La)TiO3 thin films* (submitted to Journal of Applied Physics).

## What it does
- `pipeline.py`: P_r and V_c of a loop, loop area, fit of P_r and W to `y_inf + k/f`,
  the parameter-free test `k_W / (pi V0 k_P)`, leakage compensation for ohmic or power-law leakage,
  electrode area and linear capacitance from the switching current, leakage resistivity.
- `synthetic.py`: generates loops and currents with known parameters (ohmic or superlinear leakage,
  slow switching, noise) and checks that the pipeline recovers them.
- `examples/analyse_series.py`: run the analysis on your own loops (CSV, one file per frequency).

## Quick start
```bash
pip install -r requirements.txt
python synthetic.py                 # validation table (writes synthetic.json)
pytest tests/                       # two quick checks
python examples/analyse_series.py 100:l100.csv 500:l500.csv 1000:l1k.csv 2000:l2k.csv --area 1.72e-8 --thickness 650e-9
```

## Interpreting the test
For a frequency-independent switching response plus resistive leakage, `k_W/(pi V0 k_P) = 1`.
A weak increase of V_c with frequency lowers it to about 0.92–0.95, superlinear leakage (J ~ V^3) raises it
to about 1.05–1.08, and a slow switching contribution without leakage gives 0.3–0.55 (see the article, Fig. 4).

## Data
The measured data of the article are available from the author on reasonable request.

## License
MIT
