"""Tasks B and C. Temperature-oxygen coupling at Lake Karapiro and simple bottom-water oxygen projections.

B. Daily means per profiler layer (T ~1 m, M ~13 m, B ~21.5 m). Slope of DO on T (OLS, mg/L per degC), split
   exactly into a solubility part (slope of saturation DO on T) and a deficit part (slope of -(sat - DO) on T),
   since DO = sat - deficit. Saturation: Benson & Krause (1984), freshwater, corrected to the barometric
   pressure at lake level (~53 m a.s.l., factor 0.994). 90% intervals by resampling whole ISO weeks (2000 x).
   CAVEAT: across this record, warm days are stratified days, so the deficit part of the slope mixes the effect
   of temperature on respiration with the length and strength of stratification; it is not a causal sensitivity.
   Sensor artefact removed: 22 readings where the DO channel repeats temperature (|DO - T| < 0.2 and DO > 13).

C. Bottom-layer (B) daily DO shifted by the projected change at each SSP (run_projection.warming, water
   warming relative to 2025, ~0.3 degC already realised):
     (i)   solubility only: DO' = DO + sat(T + dT) - sat(T)  (deficit unchanged)
     (ii)  solubility + respiration: deficit' = deficit x rR(T + dT)/rR(T), with rR from the MMRT form in
           unified.py at the literature ER activation energy (0.62 eV, Topt 35 degC); days with a negative
           deficit are left unscaled. Organic-matter supply and GPP are not changed.
     (iii) empirical bottom-water coupling from B: DO' = DO + b_B dT (upper bound; confounded with stratification).
   The observed daily series is shifted as a whole; this ignores any change in the timing or length of
   stratification. Writes oxygen_coupling.csv, oxygen_projection.csv, oxygen_daily_2100.csv.
"""
import sys
import numpy as np, pandas as pd

sys.path.insert(0, "../model")
import unified as u
from run_projection import warming

OUT = "../redox/"
rng = np.random.default_rng(7)
P_FACTOR = (1 - 2.25577e-5 * 53) ** 5.25588          # barometric pressure at ~53 m a.s.l. / 1 atm


def do_sat(T_C, p_factor=P_FACTOR):
    """Benson & Krause (1984) O2 saturation in fresh water (mg/L) at 1 atm, scaled to local pressure
    (water-vapour correction neglected at this small pressure difference)."""
    Tk = np.asarray(T_C, float) + 273.15
    ln = (-139.34411 + 1.575701e5 / Tk - 6.642308e7 / Tk ** 2 + 1.243800e10 / Tk ** 3 - 8.621949e11 / Tk ** 4)
    return np.exp(ln) * p_factor


import os
if os.path.exists("../karapiro/s2_readings.csv"):
    r = pd.read_csv("../karapiro/s2_readings.csv", parse_dates=["t"])
    art = ((r.DO - r["T"]).abs() < 0.2) & (r.DO > 13)
    print("artefact rows removed:", int(art.sum()))
    r = r[~art].dropna(subset=["T", "DO"]).copy()
    r["sat"] = do_sat(r["T"])
    r["day"] = r.t.dt.normalize()
    d = r.groupby(["layer", "day"]).agg(T=("T", "mean"), DO=("DO", "mean"), sat=("sat", "mean"), n=("DO", "size")).reset_index()
    d = d[d.n >= 3]
    d["deficit"] = d.sat - d.DO
    d["week"] = d.day.dt.to_period("W").astype(str)
    d.to_csv(OUT + "oxygen_daily_layers.csv", index=False)
else:   # profiler readings not deposited: start from the deposited daily layer means
    d = pd.read_csv(OUT + "oxygen_daily_layers.csv", parse_dates=["day"])


def slopes(g):
    x = g["T"].values - g["T"].mean()
    f = lambda y: float(x @ (y - y.mean()) / (x @ x))
    return f(g.DO.values), f(g.sat.values), f(-g.deficit.values)


rows = []
for L in ["T", "M", "B"]:
    g = d[d.layer == L].reset_index(drop=True)
    b, bs, bd = slopes(g)
    wk = [v.index.values for _, v in g.groupby("week")]
    bb = np.array([slopes(g.iloc[np.concatenate([wk[k] for k in rng.integers(0, len(wk), len(wk))])]) for _ in range(2000)])
    share = bb[:, 1] / bb[:, 0]
    rows.append(dict(layer=L, days=len(g), weeks=len(wk), T_range=f"{g['T'].min():.1f}-{g['T'].max():.1f}",
                     slope_DO=b, lo=np.percentile(bb[:, 0], 5), hi=np.percentile(bb[:, 0], 95),
                     slope_solubility=bs, sol_lo=np.percentile(bb[:, 1], 5), sol_hi=np.percentile(bb[:, 1], 95),
                     slope_deficit=bd, def_lo=np.percentile(bb[:, 2], 5), def_hi=np.percentile(bb[:, 2], 95),
                     share_solubility=bs / b, share_lo=np.percentile(share, 5), share_hi=np.percentile(share, 95),
                     r2=np.corrcoef(g["T"], g.DO)[0, 1] ** 2))
CPL = pd.DataFrame(rows)
CPL.to_csv(OUT + "oxygen_coupling.csv", index=False)
pd.set_option("display.width", 250)
print(CPL.round(3).to_string(index=False))

# ---------------- C. projections
bB = CPL.set_index("layer").loc["B"]
B = d[d.layer == "B"].copy()
M = u.METAB
rR = lambda T: u.mmrt_rel(T, M["Ea_R"], M["Topt_R"], M["Tref"])


def project(dT, case, slope=bB.slope_DO):
    T, DO, sat, df = B["T"].values, B.DO.values, B.sat.values, B.deficit.values
    sat2 = do_sat(T + dT)
    if case == "observed":
        return DO
    if case == "(i) solubility":
        return sat2 - df
    if case == "(ii) solubility + respiration":
        f = rR(T + dT) / rR(T)
        return sat2 - np.where(df > 0, df * f, df)
    if case == "(iii) empirical coupling":
        return DO + slope * dT
    raise ValueError(case)


def stats_(x):
    x = np.clip(x, 0, None)
    return dict(mean=x.mean(), p5=np.percentile(x, 5), min=x.min(), frac_lt2=(x < 2).mean(), frac_lt4=(x < 4).mean(),
                days_lt4=int((x < 4).sum()))


CASES = ["(i) solubility", "(ii) solubility + respiration", "(iii) empirical coupling"]
out = [dict(ssp="observed", year=2025, dT=0.0, case="observed", **stats_(B.DO.values))]
daily = {"observed": B.DO.values}
for ssp in ["SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]:
    for yr in (2050, 2100):
        dT = warming(yr, ssp)
        for c in CASES:
            x = project(dT, c)
            row = dict(ssp=ssp, year=yr, dT=dT, case=c, **stats_(x))
            if c.startswith("(iii)"):     # carry the bootstrap interval of the empirical slope
                row["mean_range"] = f"{stats_(project(dT, c, bB.hi))['mean']:.2f}-{stats_(project(dT, c, bB.lo))['mean']:.2f}"
                row["p5_range"] = f"{stats_(project(dT, c, bB.hi))['p5']:.2f}-{stats_(project(dT, c, bB.lo))['p5']:.2f}"
                row["frac_lt4_range"] = f"{stats_(project(dT, c, bB.hi))['frac_lt4']:.3f}-{stats_(project(dT, c, bB.lo))['frac_lt4']:.3f}"
            if c.startswith("(ii)"):
                row["ER_factor_mean"] = float(np.mean(rR(B['T'].values + dT) / rR(B['T'].values)))
            out.append(row)
            if yr == 2100:
                daily[f"{ssp} {c}"] = x
PRJ = pd.DataFrame(out)
PRJ.to_csv(OUT + "oxygen_projection.csv", index=False)
pd.DataFrame(daily, index=B.day.dt.date).to_csv(OUT + "oxygen_daily_2100.csv")
print(f"bottom-layer days: {len(B)} ({B.day.min().date()} to {B.day.max().date()})")
print(PRJ.round(3).to_string(index=False))

# sensitivity: one full annual cycle only (3 Feb 2025 to 2 Feb 2026) so that summer is not counted twice
B_all = B
B = B_all[(B_all.day >= "2025-02-03") & (B_all.day < "2026-02-03")]
sens = [dict(ssp="observed", case="observed", **stats_(B.DO.values))]
for ssp in ["SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]:
    for c in CASES:
        sens.append(dict(ssp=ssp, case=c, **stats_(project(warming(2100, ssp), c))))
S = pd.DataFrame(sens); S.to_csv(OUT + "oxygen_projection_oneyear.csv", index=False)
print("one annual cycle (2100):"); print(S.round(3).to_string(index=False))
