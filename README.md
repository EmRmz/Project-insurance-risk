# Riesgo de siniestro por edad

## Descripción

Proyecto final de Modelado y Simulación. Usa pólizas reales de seguros de auto de Francia (dataset freMTPL2freq) para estimar, en cada grupo de edad, la probabilidad de que un conductor tenga **al menos un siniestro en un año**. El programa en C lee los datos, calcula las tasas y hace la simulación; Python solo muestra los resultados en una gráfica.

**Modelo:** el número de siniestros de un conductor en un año se modela con una **distribución de Poisson** de tasa λ. La tasa se obtiene de los datos y la probabilidad se estima con el **método de Monte Carlo**, simulando N años con el **método de Knuth**. El resultado se compara con la probabilidad exacta para validar la simulación.

---

## Condiciones iniciales

| Condición | Valor |
| --- | --- |
| Edad válida | 18 a 100 años |
| Exposición válida | $E > 0$ |
| Tope por póliza | 4 siniestros y 1 año de exposición |
| Grupos de edad | 18-19, 20-24, …, 85-89, 90+ |
| Años simulados | $N = 10\,000$ |

---

## Tasa de siniestros

```math
\lambda = \frac{S}{E}
```

$S$ = siniestros del grupo, $E$ = años asegurados del grupo.
