# Riesgo de siniestro por edad

Probabilidad de que un conductor tenga al menos un siniestro en un año, por grupo de edad, estimada con Monte Carlo y Poisson.

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

## Modelo de Poisson

```math
X \sim \text{Poisson}(\lambda), \qquad P(X = k) = \frac{\lambda^{k}\, e^{-\lambda}}{k!}
```

---

## Probabilidad exacta

```math
P(X \ge 1) = 1 - e^{-\lambda}
```

---

## Simulación de un año (método de Knuth)

```math
U \sim \text{Uniforme}(0, 1)
```

```math
X = \min \left\{ k \ge 0 \;:\; \prod_{j=0}^{k} U_j \le e^{-\lambda} \right\}
```

---

## Estimación Monte Carlo

```math
\hat{p} = \frac{A}{N}
```

$A$ = años simulados con $X \ge 1$.
