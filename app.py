import csv
import io
import os
import subprocess
from pathlib import Path

import plotly.graph_objects as go
import streamlit as st


ROOT = Path(__file__).resolve().parent
DATA_FILE = ROOT / "data" / "raw" / "freMTPL2freq.arff"
C_PROGRAM = ROOT / ("age_risk.exe" if os.name == "nt" else "age_risk")


def run_analysis(age: int, years: int) -> list[dict[str, str]]:
    result = subprocess.run(
        [
            str(C_PROGRAM),
            "--data", str(DATA_FILE),
            "--age", str(age),
            "--years", str(years),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return list(csv.DictReader(io.StringIO(result.stdout)))


def draw_chart(rows: list[dict[str, str]], age: int) -> go.Figure:
    ages = [float(row["age_center"]) for row in rows]
    observed = [float(row["exact_probability"]) for row in rows]
    simulated = [float(row["mc_probability"]) for row in rows]
    selected = next(row for row in rows if row["selected"] == "1")

    figure = go.Figure()
    figure.add_trace(go.Scatter(
        x=ages,
        y=observed,
        mode="lines+markers",
        name="Estimación según los datos",
        line={"color": "#245b8f", "width": 3},
        marker={"size": 7},
        customdata=[[row["age_group"], row["rate"], row["policies"]] for row in rows],
        hovertemplate=(
            "Edad: %{customdata[0]} años<br>"
            "Probabilidad: %{y:.1%}<br>"
            "Tasa: %{customdata[1]} siniestros/año<br>"
            "Pólizas: %{customdata[2]}<extra></extra>"
        ),
    ))
    figure.add_trace(go.Scatter(
        x=ages,
        y=simulated,
        mode="markers",
        name="Simulación Monte Carlo",
        marker={"color": "#43a6a0", "size": 9},
        error_y={
            "type": "data",
            "symmetric": False,
            "array": [float(row["error_plus"]) for row in rows],
            "arrayminus": [float(row["error_minus"]) for row in rows],
            "color": "#43a6a0",
            "thickness": 1.5,
            "width": 4,
        },
        customdata=[[row["age_group"], row["rate"]] for row in rows],
        hovertemplate=(
            "Grupo de edad: %{customdata[0]} años<br>"
            "Probabilidad simulada: %{y:.1%}<br>"
            "Tasa: %{customdata[1]} siniestros/año<extra></extra>"
        ),
    ))
    figure.add_trace(go.Scatter(
        x=[age],
        y=[float(selected["mc_probability"])],
        mode="markers+text",
        name=f"Edad seleccionada: {age}",
        text=[f"{age} años"],
        textposition="top center",
        marker={"color": "#ed7048", "size": 16, "line": {"color": "white", "width": 2}},
        error_y={
            "type": "data",
            "symmetric": False,
            "array": [float(selected["error_plus"])],
            "arrayminus": [float(selected["error_minus"])],
            "color": "#ed7048",
            "thickness": 2,
            "width": 5,
        },
        hovertemplate=(
            "Edad seleccionada: %{x} años<br>"
            "Grupo: %{customdata[0]}<br>"
            "Probabilidad Monte Carlo: %{y:.1%}<extra></extra>"
        ),
        customdata=[[selected["age_group"]]],
    ))
    figure.update_layout(
        template="plotly_white",
        height=540,
        margin={"l": 70, "r": 30, "t": 35, "b": 65},
        xaxis={
            "title": "Edad del conductor (años)",
            "range": [18, 100],
            "dtick": 10,
            "gridcolor": "#e9edf2",
        },
        yaxis={
            "title": "Probabilidad de al menos un siniestro en un año",
            "tickformat": ".0%",
            "rangemode": "tozero",
            "gridcolor": "#e9edf2",
        },
        legend={"orientation": "h", "y": 1.12, "x": 0},
        hovermode="closest",
    )
    return figure


def main() -> None:
    st.set_page_config(page_title="Siniestros por edad", page_icon="🚗", layout="wide")
    st.title("¿Cómo cambia el riesgo de siniestro con la edad?")
    st.write(
        "El programa en C calcula las tasas por edad con las pólizas reales y ejecuta "
        "la simulación Monte Carlo. Esta página solo permite elegir los datos y mostrar la gráfica."
    )

    with st.form("conductor"):
        age_col, years_col, submit_col = st.columns([1, 1, 1])
        age = age_col.number_input("Edad del conductor", 18, 100, 30)
        years = years_col.number_input(
            "Años a simular", min_value=100, max_value=100_000, value=10_000, step=1_000
        )
        submitted = submit_col.form_submit_button("Actualizar gráfica", type="primary")

    if not DATA_FILE.exists():
        st.error(f"No se encontró el dataset: {DATA_FILE}")
        st.stop()
    if not C_PROGRAM.exists():
        st.error("No se encontró el programa C. En la terminal ejecuta: make")
        st.stop()

    if submitted or "results" not in st.session_state:
        try:
            st.session_state["results"] = run_analysis(int(age), int(years))
        except subprocess.CalledProcessError as error:
            st.error(f"El análisis en C falló: {error.stderr.strip()}")
            st.stop()
        except OSError as error:
            st.error(f"No se pudo iniciar el programa C: {error}")
            st.stop()

    rows = st.session_state["results"]
    selected = next(row for row in rows if row["selected"] == "1")
    st.caption(
        f"Edad {int(age)} años, grupo {selected['age_group']}: "
        f"{float(selected['rate']):.3f} siniestros por año según los datos. "
        f"Monte Carlo simula {int(years):,} años."
    )
    st.plotly_chart(
        draw_chart(rows, int(age)),
        use_container_width=True,
    )
    st.caption(
        "La línea muestra la probabilidad calculada en C a partir de las pólizas; "
        "los puntos muestran la proporción Monte Carlo y su intervalo aproximado del 95%."
    )


if __name__ == "__main__":
    main()
