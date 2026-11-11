"""The 8 detailed figures of the progress reports (docs/AVANCE_*.md), labelled in Spanish.

Usage (from the project root, after `make run`):
    python3 plots/report_figures.py

Reads results/*.csv and results/metrics.txt and writes PNG files to results/figures/.
All computation happens in C; this script only plots.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"

LEVELS = ["Preferential", "Standard", "Moderate Risk", "High Risk"]
LEVEL_NAMES = dict(zip(LEVELS, ["Preferencial", "Estándar", "Riesgo moderado", "Alto riesgo"]))
SCENARIO_NAMES = {"Current (-5% / +25%)": "Actual (−5% / +25%)",
                  "Stricter malus (+40%)": "Malus más estricto (+40%)",
                  "Faster bonus (-10%)": "Bonus más rápido (−10%)"}

# Palette: validated default categorical slots + an ordinal blue ramp for the 4 levels
INK, INK_2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
LEVEL_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]

# Age enters as a spline and is drawn on its own in figure 02, so it is left out here.
FEATURE_LABELS = {
    "veh_age_new": "Vehículo nuevo (0 años)", "veh_age_11plus": "Vehículo de 11+ años",
    "power_5": "Potencia clase 5", "power_6": "Potencia clase 6", "power_7": "Potencia clase 7",
    "power_8": "Potencia clase 8", "power_9": "Potencia clase 9+",
    "bonus_malus": "Bonus-malus (+100 puntos)", "log_density": "Densidad (× e)", "diesel": "Diésel",
    "exposure_under_3m": "Control: póliza < 3 meses", "exposure_3_to_9m": "Control: póliza 3–9 meses",
}


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
        "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "legend.frameon": False, "font.size": 10, "lines.linewidth": 2,
        "savefig.dpi": 180, "savefig.bbox": "tight",
    })


def read_metrics():
    metrics = {}
    for line in (RES / "metrics.txt").read_text().splitlines():
        key, *values = line.split()
        metrics[key] = [float(v) for v in values]
    return metrics


def save(fig, name):
    fig.savefig(FIG / name)
    plt.close(fig)
    print(f"  results/figures/{name}")


def plot_coefficients():
    c = pd.read_csv(RES / "model_coefficients.csv", index_col="feature")
    c = c.loc[[f for f in FEATURE_LABELS if f in c.index]].iloc[::-1]
    lo, hi = np.exp(c["coef"] - 1.96 * c["std_err"]), np.exp(c["coef"] + 1.96 * c["std_err"])
    y = np.arange(len(c))
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.axvline(1, color=INK_2, linewidth=1.2, linestyle="--", zorder=1)
    ax.hlines(y, lo, hi, color=SERIES[0], linewidth=2.2, zorder=2)
    ax.plot(c["relative_effect"], y, "o", color=SERIES[0], markersize=7,
            markeredgecolor=SURFACE, markeredgewidth=1, zorder=3)
    ax.set_yticks(y, [FEATURE_LABELS[f] for f in c.index])
    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 2, 5, 10], ["0,5", "1", "2", "5", "10"])
    ax.set_xlim(0.45, 11)
    ax.set_xlabel("Cambio multiplicativo en la frecuencia esperada (IC 95%, escala log)")
    ax.set_title("Efecto de las variables sobre el riesgo estimado", pad=30)
    ax.text(0.01, 1.01, "Valores > 1 aumentan λ; valores < 1 la reducen. La línea discontinua marca efecto neutro.",
            transform=ax.transAxes, color=INK_2, fontsize=9, va="bottom")
    ax.grid(axis="y", visible=False)
    fig.subplots_adjust(left=0.36, top=0.88)
    save(fig, "01_model_effects.png")


def plot_age_fit():
    a = pd.read_csv(RES / "age_fit.csv")
    e = pd.read_csv(RES / "age_effect.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4))
    # Dot area grows with the exposure behind each observed frequency
    ax.scatter(a["age"], a["observed_frequency"], s=4 + 60 * a["exposure"] / a["exposure"].max(),
               color=SERIES[0], alpha=0.7, edgecolor=SURFACE, linewidth=0.5, label="Observada")
    ax.plot(a["age"], a["model_frequency"], color=SERIES[1], label="Modelo (GLM)")
    ax.set_xlabel("Edad del conductor (años)")
    ax.set_ylabel("Siniestros por año de exposición")
    ax.set_title("Frecuencia de siniestros por edad")
    ax.legend(loc="upper right")
    bx.axhline(1, color=INK_2, linewidth=1)
    bx.plot(e["age"], e["full_model"], color=SERIES[0], label="Modelo completo (con bonus-malus)")
    bx.plot(e["age"], e["true_risk_model"], color=SERIES[2], label="Riesgo real (sin bonus-malus)")
    bx.set_xlabel("Edad del conductor (años)")
    bx.set_ylabel("Frecuencia relativa (45 años = 1)")
    bx.set_title("Efecto de la edad (curva suave)")
    bx.legend(loc="upper right")
    for x in (ax, bx):
        x.set_xlim(18, a["age"].max())
    save(fig, "02_frequency_by_age.png")


def plot_risk_distribution(thresholds):
    h = pd.read_csv(RES / "lambda_histogram.csv")
    fig, ax = plt.subplots(figsize=(7.5, 4))
    bottom = np.zeros(len(h))
    for name, color in zip(LEVELS, LEVEL_COLORS):
        ax.bar(h["bin_low"], h[name], width=h["bin_high"] - h["bin_low"], bottom=bottom,
               align="edge", color=color, label=LEVEL_NAMES[name], edgecolor=SURFACE, linewidth=0.3)
        bottom += h[name].to_numpy()
    for t in thresholds:
        ax.axvline(t, color=INK_2, linewidth=1, linestyle="--")
    ax.set_xscale("log")
    ax.set_xticks([0.03, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 2],
                  ["0,03", "0,05", "0,07", "0,1", "0,2", "0,5", "1", "2"])
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(0.025, 2.2)
    ax.set_ylim(0, bottom.max() * 1.08)
    ax.set_xlabel("Siniestros esperados por año, λ (escala log)")
    ax.set_ylabel("Número de conductores")
    ax.set_title("Distribución del riesgo en la población", pad=30)
    threshold_text = " · ".join(f"Q{i + 1} = {t:.3f}" for i, t in enumerate(thresholds))
    ax.text(0.01, 1.01, f"Umbrales de los cuartiles: {threshold_text}. Líneas discontinuas.",
            transform=ax.transAxes, color=INK_2, fontsize=9, va="bottom")
    ax.legend(loc="upper right", title="Nivel de riesgo")
    save(fig, "03_risk_distribution.png")


def plot_mc_convergence():
    mc = pd.read_csv(RES / "mc_convergence.csv")
    levels = pd.read_csv(RES / "risk_levels.csv", index_col="level")
    fig, ax = plt.subplots(figsize=(9, 5))
    for name, color in zip(LEVELS, LEVEL_COLORS):
        lam = levels.loc[name, "example_lambda"]
        exact = levels.loc[name, "exact_prob_claim"]
        years = mc["years"].to_numpy()
        se = np.sqrt(exact * (1 - exact) / years)
        show_band = years >= 20
        ax.fill_between(years, np.maximum(0, exact - 1.96 * se),
                        np.minimum(1, exact + 1.96 * se), where=show_band,
                        color=color, alpha=0.09, linewidth=0)
        ax.plot(years, mc[name], color=color, linewidth=1.5, alpha=0.78,
                label=f"{LEVEL_NAMES[name]} (λ = {lam:.3f})")
        ax.axhline(exact, color=color, linewidth=1, linestyle=":")
    ax.set_xscale("log")
    ax.set_xticks([1, 10, 100, 1000, 5000], ["1", "10", "100", "1.000", "5.000"])
    ax.set_xlim(1, 5000)
    ax.set_ylim(0, 0.24)
    ax.set_xlabel("Años simulados (escala logarítmica)")
    ax.set_ylabel("Probabilidad estimada de ≥1 siniestro")
    ax.set_title("La simulación converge a la probabilidad exacta", pad=30)
    ax.text(0.01, 1.01, "Línea punteada: valor exacto 1 − e⁻λ. Banda: intervalo aproximado del 95%.",
            transform=ax.transAxes, color=INK_2, fontsize=9, va="bottom")
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    ax.legend(loc="upper right", ncols=2)
    save(fig, "04_montecarlo_convergence.png")


def plot_transition_matrix():
    P = pd.read_csv(RES / "transition_matrix.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(6.2, 5))
    ax.imshow(P.to_numpy(), cmap="Blues", vmin=0, vmax=1)
    for i in range(len(P)):
        for j in range(len(P)):
            v = P.iat[i, j]
            ax.text(j, i, f"{v:.1%}", ha="center", va="center", fontsize=10,
                    color="white" if v > 0.5 else INK)
    ax.set_xticks(range(len(P)), [LEVEL_NAMES[n] for n in P.columns], rotation=20)
    ax.set_yticks(range(len(P)), [LEVEL_NAMES[n] for n in P.index])
    ax.set_xlabel("Nivel el próximo año")
    ax.set_ylabel("Nivel este año")
    ax.set_title("Probabilidad de transición en un año", pad=30)
    ax.text(0.5, 1.02, "Cada fila suma 100%; se lee desde el nivel actual hacia el siguiente.",
            transform=ax.transAxes, ha="center", color=INK_2, fontsize=9)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    save(fig, "05_transition_matrix.png")


def plot_projection(stationary):
    p = pd.read_csv(RES / "portfolio_projection.csv")
    fig, ax = plt.subplots(figsize=(9, 5))
    end = p["year"].iloc[-1]
    for i, (name, color) in enumerate(zip(LEVELS, LEVEL_COLORS)):
        ax.plot(p["year"], p[f"sim_{name}"], color=color, label=LEVEL_NAMES[name],
                marker="o", markersize=4)
        ax.plot(p["year"], p[f"markov_{name}"], color=color, linestyle="--", linewidth=1.5)
    ax.axhline(0.25, color=AXIS, linewidth=0.8, linestyle=":")
    ax.set_xlim(0, end)
    ax.set_xlabel("Año")
    ax.set_ylabel("Proporción del portafolio")
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_title("Proyección del portafolio: simulación y cadena de Markov", pad=30)
    ax.text(0.01, 1.01, "Sólida: Monte Carlo. Discontinua: π₀Pᵗ. Referencia: 25% inicial; eje vertical ampliado.",
            transform=ax.transAxes, color=INK_2, fontsize=9, va="bottom")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), title="Nivel de riesgo")
    save(fig, "06_portfolio_projection.png")


def plot_scenarios():
    s = pd.read_csv(RES / "policy_scenarios.csv", index_col="level")
    x = np.arange(len(s))
    w = 0.8 / len(s.columns)
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for k, name in enumerate(s.columns):
        ax.bar(x + (k - (len(s.columns) - 1) / 2) * w, s[name], w * 0.95, color=SERIES[k],
               label=SCENARIO_NAMES.get(name, name))
    ax.set_xticks(x, [LEVEL_NAMES[n] for n in s.index])
    ax.set_ylim(0, max(0.4, s.to_numpy().max() * 1.12))
    ax.set_ylabel("Proporción del portafolio en el año 10")
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_title("Impacto de la regla bonus-malus de la aseguradora (año 10)")
    ax.grid(axis="x", visible=False)
    for container in ax.containers:
        ax.bar_label(container, fmt="%.0f%%", labels=[f"{v.get_height():.0%}" for v in container],
                     padding=2, fontsize=8, color=INK_2)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncols=1)
    fig.subplots_adjust(bottom=0.28)
    save(fig, "07_policy_scenarios.png")


def plot_upgrade_by_age(min_drivers=300):
    u = pd.read_csv(RES / "upgrade_by_age.csv")
    u = u[u["drivers"] >= min_drivers]  # fewer real drivers than this: too noisy to show
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    # Each panel is titled with its level, so the same two colors serve all three panels
    for ax, name in zip(axes, LEVELS[:-1]):
        g = u[u["level"] == name]
        ages = g["age"].to_numpy()
        total = 100 * g["prob_up"].to_numpy()
        aging = 100 * g["prob_up_without_claims"].to_numpy()
        breaks = np.flatnonzero(np.diff(ages) > 1) + 1
        for segment in np.split(np.arange(len(g)), breaks):
            ax.fill_between(ages[segment], aging[segment], color=LEVEL_COLORS[0],
                            alpha=0.38, linewidth=0)
            ax.plot(ages[segment], total[segment], color=SERIES[0], linewidth=2)
        ax.set_title(f"Nivel actual: {LEVEL_NAMES[name]}")
        ax.set_xlabel("Edad actual del conductor (años)")
        ax.set_xlim(18, 90)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=0))
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("Probabilidad de subir el próximo año")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([0], [0], color=SERIES[0], linewidth=2,
                      label="Probabilidad total de subir"),
               Patch(facecolor=LEVEL_COLORS[0], alpha=0.38,
                     label="Sube solo por cumplir un año, sin siniestros")]
    fig.legend(handles=handles, loc="lower center", ncols=2, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("¿Qué tan probable es subir de nivel de riesgo en un año? Por edad y nivel actual",
                 x=0.01, ha="left", fontweight="bold", color=INK)
    fig.text(0.01, 0.9, "Conductores reales de 18 años o más (freMTPL2), cada uno simulado un año "
                        f"10 veces. Se omiten edades con menos de {min_drivers} conductores en el nivel.",
             color=INK_2, fontsize=9)
    fig.subplots_adjust(top=0.8, bottom=0.2, wspace=0.08)
    save(fig, "08_upgrade_by_age.png")


def main():
    if not (RES / "metrics.txt").exists():
        raise SystemExit("results/metrics.txt not found; run `make run` first.")
    FIG.mkdir(parents=True, exist_ok=True)
    style()
    metrics = read_metrics()
    plot_coefficients()
    plot_age_fit()
    plot_risk_distribution(metrics["thresholds"])
    plot_mc_convergence()
    plot_transition_matrix()
    plot_projection(metrics["stationary_distribution"])
    plot_scenarios()
    plot_upgrade_by_age()


if __name__ == "__main__":
    main()
