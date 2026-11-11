"""Probability of each risk level by age, for a driver who starts at 18.

Usage (from the project root, after `make run`):
    python3 plots/plot_results.py           # opens the figure in a window
    python3 plots/plot_results.py --save    # writes results/figures/levels_by_age.png instead

Reads results/levels_by_age.csv, written by the C program (Monte Carlo simulation of
100,000 drivers followed from age 18 to 80). Only drivers aged 18 or more are counted:
the C program drops younger records from the data, and this script plots ages >= 18 only.
This script only plots.
"""

import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"

# Bottom to top: the worst level sits on the axis, so the band of new drivers starts at 100%
LEVELS = ["High Risk", "Moderate Risk", "Standard", "Preferential"]
NAMES = ["Alto riesgo", "Riesgo moderado", "Estándar", "Preferencial"]
COLORS = ["#0d366b", "#1c5cab", "#3987e5", "#86b6ef"]  # ordinal blue ramp, dark = more risk
INK, INK_2, MUTED, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7", "#fcfcfb"
MIN_AGE = 18  # only adult drivers are counted


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--save", action="store_true", help="write a PNG instead of opening a window")
    args = parser.parse_args()
    if args.save:
        matplotlib.use("Agg")

    path = RES / "levels_by_age.csv"
    if not path.exists():
        raise SystemExit(f"{path.relative_to(ROOT)} not found; run `make run` first.")
    import pandas as pd
    d = pd.read_csv(path)
    d = d[d["age"] >= MIN_AGE]
    age = d["age"]
    pct = [100 * d[level] for level in LEVELS]

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.labelcolor": INK_2,
                         "axes.grid": True, "axes.axisbelow": True,
                         "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(11, 6), facecolor=SURFACE)
    fig.canvas.manager.set_window_title("Nivel de riesgo por edad")
    ax.set_facecolor(SURFACE)
    ax.stackplot(age, pct, colors=COLORS, labels=NAMES, edgecolor=SURFACE, linewidth=0.6)

    ax.set_xlim(age.iloc[0], age.iloc[-1])
    ax.set_ylim(0, 100)
    ax.set_xticks([MIN_AGE, *range(20, int(age.iloc[-1]) + 1, 10)])
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=0))
    ax.set_xlabel("Edad del conductor (años)")
    ax.set_ylabel("Conductores en cada nivel")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", visible=False)
    ax.set_title("Cómo cambia el nivel de riesgo a lo largo de la vida",
                 loc="left", fontweight="bold", color=INK, pad=30)
    ax.text(0, 1.02, "Simulación de 100.000 conductores que comienzan a los 18 años; modelo basado en 678.013 pólizas de Francia.",
            transform=ax.transAxes, color=INK_2, fontsize=9.5, va="bottom")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=False, title="Nivel de riesgo", alignment="left")

    if args.save:
        out = RES / "figures" / "levels_by_age.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=180, bbox_inches="tight")
        print(f"  {out.relative_to(ROOT)}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
