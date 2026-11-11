# Project-insurance-risk

Proyecto académico que modela la frecuencia de siniestros de autos y simula cómo cambia el nivel de riesgo de los conductores con el tiempo.

## ¿Qué modela?

El proyecto estima cuántos siniestros se esperan por año para un perfil de conductor. Usa datos de pólizas francesas del conjunto [freMTPL2freq de OpenML](https://www.openml.org/d/41214) y considera características como:

- Edad del conductor.
- Antigüedad y potencia del vehículo.
- Coeficiente bonus-malus, que resume el historial de siniestros.
- Densidad de población y tipo de combustible.

El modelo estima la **frecuencia** de siniestros, no el costo de las reparaciones ni la severidad de los accidentes.

## Modelo de frecuencia

Se usa un modelo lineal generalizado Poisson (GLM). Para una póliza \(i\), el número de siniestros \(N_i\) se modela como:

\[
N_i \sim \operatorname{Poisson}(e_i\lambda_i)
\]

donde \(e_i\) es el tiempo de exposición de la póliza en años y \(\lambda_i\) es la frecuencia anual esperada para ese perfil.

La frecuencia se calcula con:

\[
\log(\lambda_i)=\beta_0+\sum_{j=1}^{p}\beta_jx_{ij},
\qquad
\lambda_i=\exp\left(\beta_0+\sum_{j=1}^{p}\beta_jx_{ij}\right)
\]

Las variables \(x_{ij}\) representan las características del conductor y del vehículo; los coeficientes \(\beta_j\) se estiman a partir de los datos. La edad se representa con una curva suave por tramos. El ajuste se realiza en C mediante IRLS, un método iterativo equivalente a aplicar Newton para maximizar la verosimilitud del modelo.

La frecuencia también se puede expresar como la probabilidad de tener al menos un siniestro durante un año:

\[
P(N_i\geq 1)=1-\exp(-\lambda_i)
\]

## Variación entre conductores

Conductores con características similares pueden tener resultados distintos. Para representar esa variación, las simulaciones usan un multiplicador aleatorio \(\theta_i\):

\[
\theta_i\sim\operatorname{Gamma}\left(\frac{1}{\phi},\ \text{escala}=\phi\right),
\qquad
N_i\mid\theta_i\sim\operatorname{Poisson}(\theta_i\lambda_i)
\]

Así, \(\theta_i\) tiene media 1 y varianza \(\phi\), que se estima con los datos. Esta mezcla de distribuciones Poisson y Gamma produce más variabilidad entre conductores que una Poisson simple.

## Niveles de riesgo

El proyecto ordena las frecuencias estimadas \(\lambda\) de la población y usa sus cuartiles para formar cuatro grupos:

1. Preferencial.
2. Estándar.
3. Riesgo moderado.
4. Alto riesgo.

Por eso, los niveles son **relativos a la población analizada**: no son umbrales universales de peligrosidad ni una tarifa de seguro.

## Simulación de cambios en el tiempo

Para simular una cartera, el programa repite estos pasos cada año:

1. Genera el número de siniestros de cada conductor usando su frecuencia estimada.
2. Actualiza su bonus-malus: baja 5% si no tuvo siniestros y sube 25% por cada siniestro, respetando los límites definidos por la regla.
3. Aumenta la edad del conductor y vuelve a calcular su nivel de riesgo.

La simulación separa el riesgo usado para generar siniestros del nivel que la aseguradora asigna. En las trayectorias de vida, los siniestros se generan con un modelo sin bonus-malus; luego el historial actualizado sí se considera para volver a clasificar al conductor. Esto evita que el bonus-malus por sí solo cree una retroalimentación artificial en el riesgo subyacente.

## Cadena de Markov

Los cuatro niveles se usan como estados de una cadena de Markov. La probabilidad de pasar del nivel \(i\) al nivel \(j\) en un año se representa por \(P_{ij}\). Se estima contando las transiciones observadas en la simulación y dividiendo cada conteo entre el total de conductores que comenzaron en el nivel \(i\).

La distribución de la cartera después de \(t\) años se proyecta mediante:

\[
\boldsymbol{\pi}_t=\boldsymbol{\pi}_0P^t
\]

donde \(\boldsymbol{\pi}_0\) es la distribución inicial de conductores entre niveles. Esta proyección se compara con la simulación Monte Carlo. La cadena es una aproximación: dos conductores del mismo nivel pueden tener distintas edades y características, por lo que sus probabilidades futuras no necesariamente son iguales.

## Alcance e interpretación

Los datos describen pólizas de Francia y el modelo estima frecuencia, no costo total de siniestros. Las simulaciones muestran resultados posibles bajo los supuestos del modelo; no son predicciones garantizadas sobre una persona ni deben interpretarse como una evaluación universal del riesgo.

Los informes detallados del proyecto están en [docs/AVANCE_1.md](docs/AVANCE_1.md) y [docs/AVANCE_2.md](docs/AVANCE_2.md).
