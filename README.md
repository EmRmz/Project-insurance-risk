# Riesgo de siniestro por edad

Aplicación sencilla que usa las pólizas reales de `data/raw/freMTPL2freq.arff` para estimar la frecuencia de siniestros por grupo de edad y mostrar una simulación Monte Carlo.

## Cómo funciona

1. El programa en C lee las pólizas y agrupa a los conductores adultos por edad.
2. Para cada grupo calcula la tasa observada: total de siniestros dividido entre el total de años asegurados (exposición).
3. Con esa tasa, C simula siniestros durante la cantidad de años indicada y calcula la proporción de años con al menos un siniestro.
4. Python muestra los resultados en una gráfica interactiva y resalta el grupo de edad elegido.

La línea de la gráfica representa la probabilidad calculada en C a partir de la tasa observada. Los puntos representan la proporción obtenida por Monte Carlo; las barras muestran un intervalo aproximado del 95%.

La edad selecciona un grupo de datos; no es una predicción individual. El análisis describe pólizas de Francia y no considera el costo de los siniestros.

## Archivos

- `src/age_risk.c`: lectura de datos, tasas por edad y simulación Monte Carlo.
- `app.py`: interfaz y gráfica; no calcula tasas ni simula siniestros.
- `tests/test_age_risk.c`: pruebas del programa en C.
- `data/raw/freMTPL2freq.arff`: datos utilizados.

## Ejecutar en WSL

Instala las dependencias de Python y compila el programa en C:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
make
```

Inicia la aplicación:

```bash
python -m streamlit run app.py
```

Para ejecutar las pruebas del código C:

```bash
make test
```
