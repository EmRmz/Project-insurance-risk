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

---

## Probabilidad de tener al menos un siniestro

```math
P = 1 - e^{-\lambda}
```

$e^{-\lambda}$ es la probabilidad de **no** tener ningún siniestro en el año; al restarla de 1 queda la probabilidad de tener uno o más.

---

## Estimación por simulación

```math
\hat{p} = \frac{A}{N}
```

$A$ = años simulados con al menos un siniestro, $N$ = años simulados en total.

---

## Cómo funciona la simulación

1. Con los datos se calcula la tasa $\lambda$ de cada grupo de edad.
2. Se simula un año: la computadora genera números al azar entre 0 y 1 y los multiplica uno tras otro. Cada multiplicación que deja el producto por encima de $e^{-\lambda}$ cuenta como un siniestro (método de Knuth).
3. Si en ese año hubo al menos un siniestro, se suma 1 a $A$.
4. Se repiten los pasos 2 y 3 durante $N$ años (por ejemplo, 10 000).
5. Se divide $A$ entre $N$ para obtener $\hat{p}$.
6. Se compara $\hat{p}$ con el valor exacto $P$: si son casi iguales, la simulación es correcta. Mientras más años se simulan, más se acercan.
