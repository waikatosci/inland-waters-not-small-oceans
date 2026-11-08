"""Fig. 1 of the GCB manuscript: drivers, controls in the water and the ecosystem response, with magnitudes
from this study (after the layout of Battin et al. 2023, Nature, Fig. 3)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
})

# palette (after Battin et al. 2023, Fig. 3)
BAND = "#E6E9EE"
DRV = "#F6B77B"      # global-change drivers (light orange)
DRV_D = "#EE7A2B"    # darker orange
CTRL = "#1A9E9E"     # in-water controls (teal)
ECO = "#9ED9D6"      # ecosystem compartments (light teal)
LOOP = "#E0601E"     # reinforcing loop arrows
GREY = "#7A8591"
TXT = "#1E252B"

fig = plt.figure(figsize=(7.2, 6.9), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 94)
ax.axis("off")

# bands
bands = [(78, 93, "Global change"), (54, 76, "Controls in\nthe water"), (1, 52, "Inland-water\necosystem")]
for y0, y1, lab in bands:
    ax.add_patch(Rectangle((7, y0), 92, y1 - y0, fc=BAND, ec="none", zorder=0))
    ax.text(3.6, (y0 + y1) / 2, lab, rotation=90, ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=TXT, linespacing=1.1)


def box(x, y, w, h, title, body, fc, tc="white", fs=7.0, ls="-", ec=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                fc=fc, ec=ec or fc, lw=1.2, ls=ls, zorder=2))
    ax.text(x + w / 2, y + h - 1.4, title, ha="center", va="top", fontsize=fs + 0.8,
            fontweight="bold", color=tc, zorder=3, linespacing=1.05)
    ax.text(x + w / 2, y + h - (7.0 if "\n" in title else 4.4), body, ha="center", va="top", fontsize=fs, color=tc,
            linespacing=1.25, zorder=3)


def circ(x, y, r, title, body, fs=6.9):
    ax.add_patch(Circle((x, y), r, fc=ECO, ec="none", zorder=2))
    ax.text(x, y + r * 0.45, title, ha="center", va="center", fontsize=fs + 0.7, fontweight="bold",
            color=TXT, zorder=3, linespacing=1.05)
    ax.text(x, y - r * 0.17, body, ha="center", va="center", fontsize=fs, color=TXT, zorder=3,
            linespacing=1.25)


def arrow(p0, p1, color=GREY, lw=1.3, ls="-", rad=0.0, z=1, ms=9):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, color=color, lw=lw, ls=ls,
                        connectionstyle=f"arc3,rad={rad}", zorder=z, shrinkA=0, shrinkB=0)
    ax.add_patch(a)


def label(x, y, s, color=GREY, fs=6.4, style="italic", ha="center", rot=0):
    ax.text(x, y, s, ha=ha, va="center", fontsize=fs, color=color, style=style, zorder=4,
            rotation=rot, bbox=dict(fc=BAND, ec="none", pad=0.5))


# ---- Global change ----
box(10, 80, 25, 10.5, "Rising atmospheric CO$_2$",
    "Sets only ~1/S of the pH response\nwhile the CO$_2$ excess is fixed", DRV, tc=TXT)
box(39.5, 80, 25, 10.5, "Warming",
    "Speeds respiration more\nthan photosynthesis", DRV, tc=TXT)
box(69.5, 80, 26.5, 10.5, "Biological invasions",
    "Gold clam, hornwort: take up Ca,\nmove nutrients, lower GPP/ER", DRV_D, tc="white", fs=6.6)

# ---- Controls in the water ----
box(10, 56.5, 19.5, 16, "Gas exchange", "Median river (S = 4.6)\nundergoes 34–54% of\nthe ocean-like pH\nresponse by 2100", CTRL, fs=6.6)
box(32, 56.5, 19.5, 16, "Surface charge", "Lower pH loosens metal\nbinding: divalent metals\n+1–12% (median data),\nup to 2.1× at the limit", "white",
    tc=TXT, fs=6.6, ls="--", ec=GREY)
box(54, 56.5, 19.5, 16, "Net metabolic source", "ER − GPP grows\n10–13% per °C; adds\n42–75% to forcing at\nLake Karāpiro", CTRL, fs=6.6)
box(76, 56.5, 20, 16, "Solubility and\nstratification", "Bottom-water O$_2$\n−0.42 mg L$^{-1}$ per °C\n(seasonal; upper bound)", CTRL, fs=6.6)

# ---- Ecosystem ----
circ(20, 38, 9.4, "Carbon supply", "Free CO$_2$ +14–48%;\nair supply in\nblooms +35–142%")
circ(50, 34, 11.0, "Metabolism", "ER outpaces GPP:\nmore CO$_2$, less O$_2$")
circ(84, 38, 9.4, "Deoxygenation", "Mean bottom O$_2$\n−9–20%; days below\n5 mg L$^{-1}$ up to 7×")
circ(84, 12.5, 9.4, "Redox release\nof nutrients", "Mn, Fe, Co, P\n×10–100 in anoxic\nbottom water")
circ(50, 10.2, 8.8, "Primary\nproducers", "potentially harmful\ncyanobacteria vs\nalgae, macrophytes", fs=6.5)
circ(20, 12.5, 9.4, "Carbonate\nsaturation", "Calcite Ω ≈ 0.03;\nfalls a further tenth\n(SSP2-4.5) to over a\nquarter (SSP5-8.5)")

# ---- drivers -> controls ----
arrow((19, 80), (19, 73.2), color=DRV_D, lw=1.6)
arrow((26, 80), (39, 73.2), color=DRV_D, lw=1.2, ls="--")
arrow((52, 80), (62, 73.2), color=DRV_D, lw=1.6)
arrow((58, 80), (83, 73.2), color=DRV_D, lw=1.6)
arrow((76, 80), (69, 73.2), color=DRV_D, lw=1.6)

# ---- controls -> ecosystem ----
arrow((19.5, 56.5), (20, 47.6), color=CTRL, lw=1.5)
arrow((12.5, 56.5), (12.6, 18.6), color=CTRL, lw=1.2, rad=0.22)
label(10.6, 30, "pH, CO$_3^{2-}$", color=CTRL, rot=90)
arrow((62, 56.5), (55.5, 43.7), color=CTRL, lw=1.5)
arrow((86, 56.5), (85, 47.6), color=CTRL, lw=1.5)
# surface charge (weak) -> producers, routed through the gap
arrow((38, 56.5), (42.6, 15.2), color=GREY, lw=1.1, ls="--", rad=0.18)
label(35.6, 47.5, "metals\n(uncertain)", color=GREY)

# ---- reinforcing loop ----
arrow((60.8, 36), (74.8, 38.8), color=LOOP, lw=2.2, ms=11)
label(67.6, 40.6, "O$_2$ loss", color=LOOP, fs=6.8)
arrow((84, 28.6), (84, 22.2), color=LOOP, lw=2.2, ms=11)
arrow((74.7, 11.6), (59.0, 10.6), color=LOOP, lw=2.2, ms=11, ls=(0, (4, 2)))
label(66.8, 7.6, "metals, P\n(hypothesised)", color=LOOP, fs=6.3)
arrow((50, 19.0), (50, 23.0), color=LOOP, lw=2.2, ms=11)
label(56.4, 21.0, "organic matter", color=LOOP, fs=6.3)
ax.text(68.5, 25.5, "reinforcing\nloop", ha="center", va="center", fontsize=8, fontweight="bold",
        color=LOOP, zorder=4)

# carbon supply -> producers; negative feedback to metabolism
arrow((26.8, 31.5), (42.2, 14.6), color=CTRL, lw=1.3)
label(31.0, 23.0, "carbon", color=CTRL)
arrow((29.4, 38.5), (39.0, 36.2), color=GREY, lw=1.1, ls="--")
label(34.2, 40.4, "raises GPP/ER (−)", color=GREY, fs=6.0)

fig.savefig("figures/fig_concept.png", dpi=300)
fig.savefig("figures/fig_concept.pdf")
print("ok")
