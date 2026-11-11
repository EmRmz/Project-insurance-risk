# Project-insurance-risk

Proyecto académico que estima la frecuencia de siniestros de autos y simula cómo puede cambiar el nivel de riesgo de una cartera de conductores.

## En pocas palabras

1. Usa datos reales de seguros de auto de Francia.
2. Aprende qué características se relacionan con más o menos siniestros.
3. Clasifica a los conductores en cuatro niveles: Preferencial, Estándar, Riesgo moderado y Alto riesgo.
4. Simula posibles resultados futuros y compara la simulación con una cadena de Markov.

El programa de análisis y el modelo están escritos en C. Python se usa únicamente para crear las gráficas.

## Datos y límites del modelo

Se usa el conjunto [freMTPL2freq de OpenML](https://www.openml.org/d/41214), con pólizas reales de Francia. El modelo estima la **frecuencia** de siniestros, no su costo.

Las categorías se definen usando cuartiles de riesgo de esa población: sirven para comparar conductores dentro de los datos, pero no son una medida universal de peligro ni determinan una prima real. Los resultados son simulaciones, no predicciones garantizadas sobre una persona.

## Cómo ejecutarlo

Necesitas un compilador C (GCC o Clang), `make`, `curl` y Python 3. La primera ejecución descarga los datos. Para crear las tablas y el modelo:

```bash
make run
```

Para compilar y correr las pruebas:

```bash
make test
```

Para instalar los paquetes de las gráficas y generarlas:

```bash
python -m pip install -r plots/requirements.txt
make figures
```

Las imágenes se guardan en `results/figures/`. `make figures` debe ejecutarse después de `make run`, porque usa las tablas generadas por el análisis.

## Evaluar un conductor

Primero genera el modelo y compila el programa:

```bash
make run
make assess_driver
```

Ejemplo en Windows PowerShell:

```powershell
.\assess_driver.exe --age 23 --veh-age 3 --bonus-malus 90 --density 1200 --horizon 10
```

En Linux o macOS, usa `./assess_driver` en lugar de `.\assess_driver.exe`.

El programa muestra la frecuencia esperada de siniestros por año, su categoría actual y una proyección. Puedes cambiar los datos con estas opciones:

| Opción | Significado | Valor predeterminado |
|---|---|---:|
| `--age` | Edad del conductor (requerida) | — |
| `--veh-age` | Antigüedad del vehículo en años (requerida) | — |
| `--veh-power` | Clase de potencia del vehículo | 6 |
| `--bonus-malus` | Coeficiente de historial; 50 es bajo y 100 equivale a conductor nuevo | 100 |
| `--density` | Densidad de población en habitantes/km² | 1000 |
| `--gas` | Tipo de combustible: `Regular` o `Diesel` | `Regular` |
| `--horizon` | Años de proyección, entre 0 y 50 | 10 |
| `--years` | Simulaciones Monte Carlo | 10000 |
| `--seed` | Semilla para repetir la simulación | Hora actual |

## Archivos principales

- `src/`: modelo, datos, simulación y cadena de Markov en C.
- `plots/`: scripts Python para generar las gráficas.
- `results/`: tablas CSV, modelo ajustado y gráficas generadas.
- `tests/`: pruebas del modelo y la simulación.
- `docs/AVANCE_1.md` y `docs/AVANCE_2.md`: informes detallados del proyecto y sus resultados.
