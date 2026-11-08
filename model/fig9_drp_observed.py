"""Fig. 9. Observed phosphorus and pH trends between Lake Karapiro and Hamilton compared with
the electrostatic expectation. Waikato Regional Council monthly monitoring (ECP compilation,
data/ECP_ChemData_JL_norm_to_flow.xlsx): Hamilton Narrows (1131-328, 2007-2022) and Hamilton
Traffic Bridge (1131-64, 1989-2015). Karapiro tailrace (1131-79) has DRP only from 2018 and is
not trended.
(a) Flow-adjusted dissolved reactive phosphorus (DRP), annual medians, as % departure from each
site's record mean, with the seasonal Sen slope. Black dashed: the change expected from the
electrostatic term alone (Fig. 8 reconstruction), drawn through the Traffic Bridge record mean.
(b) Seasonal Sen slopes of flow-adjusted DRP, total phosphorus (TP) and DRP:TP, % per decade,
with 90% bootstrap intervals (resampling years). Vertical band: model expectation for DRP and
DRP:TP. (c) Flow-normalised pH (WRTDS/EGRET, from the 2022 proposal analysis) as change from 1990
at Waipapa tailrace (1131-143, above Karapiro) and Hamilton Traffic Bridge (national network site HM3,
held as 1131-64; sampling ended January 2016, so the record stops at 2015),
against the modelled change at Lake Karapiro (Fig. 8). Ticks on (a): WRC DRP method changes.
Flow adjustment: residuals of log concentration (or pH) regressed on log discharge, per site."""
import numpy as np
import pandas as pd
from fig_common import *

rng = np.random.default_rng(1)
d = pd.read_excel("data/ECP_ChemData_JL_norm_to_flow.xlsx", "ECP_ChemData_incl_TE xlsx")
SITES = {"1131_328": "Hamilton Narrows", "1131_64": "Hamilton Traffic Bridge"}
SITE_COL = {"1131_328": COL_TEAL_400, "1131_64": COL_TEAL_800}
VARS = [("_DRP", "DRP", True), ("_TP", "TP", True), ("ratio", "DRP:TP", True), ("_pH", "pH", False)]

rec = pd.read_csv("results/reconstruction_karapiro.csv")
rec = rec[rec.year <= 2025].drop_duplicates("year").set_index("year")
D_OF = lambda y: np.interp(y, rec.index, rec.D)
PH_OF = lambda y: np.interp(y, rec.index, rec.pH)
Z_PO4 = -1.5


def seasonal_sen(t, y, m):
    sl = []
    for mm in np.unique(m):
        i = m == mm
        tt, yy = t[i], y[i]
        for a in range(len(tt) - 1):
            dt = tt[a + 1:] - tt[a]
            ok = dt != 0
            sl.append(((yy[a + 1:] - yy[a])[ok] / dt[ok]))
    sl = np.concatenate(sl) if sl else np.array([np.nan])
    return np.median(sl)


def prep(site, var, log):
    x = d[d.DeviceID == site].dropna(subset=["_Discharge"]).copy()
    x["ratio"] = x._DRP / x._TP
    x = x.dropna(subset=[var])
    x = x[x[var] > 0]
    x["t"] = x.Date.dt.year + (x.Date.dt.dayofyear - 0.5) / 365.25
    x["m"] = x.Date.dt.month
    x["yr"] = x.Date.dt.year
    y = np.log(x[var]) if log else x[var]
    lq = np.log(x._Discharge)
    x["r"] = y - np.polyval(np.polyfit(lq, y, 1), lq)
    return x


def trend(x, log, nboot=400):
    f = lambda z: seasonal_sen(z.t.values, z.r.values, z.m.values) * 10
    est = f(x)
    yrs = x.yr.unique()
    bs = []
    for _ in range(nboot):
        pick = rng.choice(yrs, len(yrs))
        z = pd.concat([x[x.yr == y].assign(t=x[x.yr == y].t + 1000 * k) if False else x[x.yr == y]
                       for k, y in enumerate(pick)])
        bs.append(f(z))
    lo, hi = np.nanpercentile(bs, [5, 95])
    conv = (lambda s: 100 * (np.exp(s) - 1)) if log else (lambda s: s)
    return conv(est), conv(lo), conv(hi)


rows = []
for site, name in SITES.items():
    for var, lab, log in VARS:
        x = prep(site, var, log)
        e, lo, hi = trend(x, log)
        rows.append(dict(site=name, var=lab, est=e, lo=lo, hi=hi, y0=x.yr.min(), y1=x.yr.max(), n=len(x)))
tr = pd.DataFrame(rows)
tr.to_csv("results/observed_trends_karapiro_hamilton.csv", index=False)

# model expectation per decade over 1990-2025
dD = (D_OF(2025) - D_OF(1990)) / 3.5
exp_drp = 100 * (10 ** (Z_PO4 * dD) - 1)
exp_pH = (PH_OF(2025) - PH_OF(1990)) / 3.5

fig, axs = plt.subplots(1, 3, figsize=(DOUBLE_COL, SINGLE_COL * 0.62),
                        gridspec_kw=dict(width_ratios=[1.15, 0.95, 1.15]))

# a: time series
ax = axs[0]
for site, name in SITES.items():
    x = prep(site, "_DRP", True)
    ann = x.groupby("yr").r.median()
    ann = ann[x.groupby("yr").size() >= 6]
    ax.plot(ann.index, 100 * (np.exp(ann - x.r.mean()) - 1), "o", ms=2.2, color=SITE_COL[site],
            mec="none", label=name)
    e = tr[(tr.site == name) & (tr["var"] == "DRP")].iloc[0]
    yy = np.array([e.y0, e.y1]); mid = (e.y0 + e.y1) / 2
    ax.plot(yy, 100 * ((1 + e.est / 100) ** ((yy - mid) / 10) - 1), color=SITE_COL[site], lw=0.9)
yy = np.arange(1989, 2026)
ax.plot(yy, 100 * (10 ** (Z_PO4 * (D_OF(yy) - D_OF(2002.5))) - 1), color=COL_BG_900, lw=0.8,
        ls=(0, (3, 1.5)), label="Electrostatic expectation")
ax.axhline(0, color=COL_BG_300, lw=0.4)
ax.set_xlim(1988, 2024)
ax.set_xlabel("Year")
ax.set_ylabel("Flow-adjusted DRP,\ndeparture from record mean (%)")
ax.legend(loc="upper left", frameon=False)
panel_label(ax, "a", x=0.92)
style_ax(ax)

# b: trend dots
ax = axs[1]
labs = ["DRP", "TP", "DRP:TP"]
ax.axvspan(exp_drp - 0.4, exp_drp + 0.4, color=COL_BG_300, lw=0, zorder=0)
ax.axvline(0, color=COL_BG_300, lw=0.4)
for j, (site, name) in enumerate(SITES.items()):
    for i, v in enumerate(labs):
        e = tr[(tr.site == name) & (tr["var"] == v)].iloc[0]
        y = len(labs) - 1 - i + (0.14 if j == 0 else -0.14)
        ax.plot([e.lo, e.hi], [y, y], color=SITE_COL[site], lw=0.8)
        ax.plot(e.est, y, "o", ms=3, color=SITE_COL[site], mec="none")
ax.set_yticks(range(len(labs)), labs[::-1])
ax.set_ylim(-0.6, len(labs) - 0.4)
ax.set_xlabel("Trend (% per decade)")
ax.text(exp_drp - 0.8, -0.5, "Model", ha="right", va="bottom", fontsize=5.5, color=COL_BG_700)
panel_label(ax, "b", x=0.92, y=0.12)
style_ax(ax)

# c: flow-normalised pH from the proposal analysis (EGRET), change from 1990
ax = axs[2]
fnx = pd.ExcelFile("data/Flow_Normalized_pH_Selected.xlsx")
FN = {"1131-143": ("Waipapa (above Kar\u0101piro)", COL_BROWN_400), "1131-64 and HM3": ("Hamilton", COL_TEAL_800)}
fn_tr = {}
for sh, (lab, col) in FN.items():
    f = fnx.parse(sh)
    f = f[pd.to_numeric(f.Year, errors="coerce").notna()].astype({"Year": float, "FN pH": float})
    if lab == "Hamilton":   # HM3 sampling ended January 2016; later flow-normalised values rest on flow alone
        f = f[f.Year <= 2015]
    base = np.interp(1990, f.Year, f["FN pH"])
    ax.plot(f.Year, f["FN pH"] - base, "o-", ms=2, lw=0.6, color=col, mec="none", label=lab)
    fn_tr[lab] = np.polyfit(f.Year, f["FN pH"], 1)[0] * 10
yy = np.arange(1988, 2023)
ax.plot(yy, PH_OF(yy) - PH_OF(1990), color=COL_BG_900, lw=0.8, ls=(0, (3, 1.5)), label="Model (CO$_2$ and warming)")
ax.axhline(0, color=COL_BG_300, lw=0.4)
ax.set_xlim(1987, 2023)
ax.set_xlabel("Year")
ax.set_ylabel("Flow-normalised pH,\nchange from 1990")
ax.text(2004, 0.035, "Waipapa (above Kar\u0101piro)", color=COL_BROWN_400, fontsize=6.5, ha="center", va="center")
ax.text(2015.7, -0.208, "Hamilton", color=COL_TEAL_800, fontsize=6.5, ha="left", va="center")
ax.text(2022.5, -0.037, "Model (CO$_2$ and warming)", color=COL_BG_900, fontsize=6.5, ha="right", va="top")
panel_label(ax, "c", x=0.03)
style_ax(ax)
for yv in (2002.7, 2005.4, 2017.5):
    axs[0].plot([yv, yv], [0.97, 1.0], transform=axs[0].get_xaxis_transform(), color=COL_BG_600, lw=0.6, clip_on=False)

fig.tight_layout(w_pad=1.4)
save(fig, "fig9_drp_observed")

pd.set_option("display.width", 200)
print(tr.round(3).to_string(index=False))
print(f"model 1990-2025: dD/dec {dD:+.4f}, DRP {exp_drp:+.1f}%/dec, pH {exp_pH:+.4f}/dec")
for name in SITES.values():
    e = tr[(tr.site == name) & (tr["var"] == "pH")].iloc[0].est
    # DRP change implied by the observed pH trend (PZC change from warming ~ -0.003/dec, ignored)
    print(f"{name}: DRP implied by observed pH trend {100*(10**(Z_PO4*(-e))-1):+.1f}%/dec")
for k, v in fn_tr.items():
    print(f"FN pH {k}: {v:+.3f}/dec; DRP implied {100*(10**(Z_PO4*(-v))-1):+.1f}%/dec")
