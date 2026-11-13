# Riesgo de siniestro por edad

Proyecto final de Modelado y Simulación. Con pólizas reales de seguros de auto (`data/raw/freMTPL2freq.arff`) se estima, para cada grupo de edad, la probabilidad de que un conductor tenga al menos un siniestro en un año, usando el método de Monte Carlo. El programa en C hace la lectura y los cálculos; Python solo muestra los resultados en una gráfica.

## Pasos de modelado

1. **Observar (datos).** Se usa el dataset freMTPL2freq de Francia. De cada póliza se toman tres columnas: la edad del conductor (`DrivAge`), el número de siniestros (`ClaimNb`) y el tiempo asegurado en años (`Exposure`). Solo se usan conductores de 18 a 100 años con exposición mayor a 0.

2. **Marco teórico.** Monte Carlo obtiene una respuesta aproximada repitiendo muchas veces un experimento con números aleatorios y promediando. Para contar eventos que ocurren al azar en un periodo de tiempo (como choques en un año) se usa la distribución de Poisson, que depende de un solo parámetro: la tasa λ.

3. **Modelo matemático.** Las edades se agrupan en intervalos (18-19, 20-24, …, 90+). Para cada grupo:

   - Tasa observada: λ = S / E
   - Siniestros en un año: X ~ Poisson(λ)
   - Probabilidad exacta de al menos un siniestro: P = 1 − e^(−λ)

   Este valor exacto sirve para comprobar que la simulación da un resultado correcto.

4. **Métodos numéricos.** Para cada grupo se simulan N años:

   - Se generan números aleatorios U entre 0 y 1 con `rand()`.
   - Cada año se simula con el método de Knuth: se multiplican números U hasta que el producto es menor que e^(−λ); el número de multiplicaciones extra es X, la cantidad de siniestros de ese año.
   - Se cuenta A, los años con X ≥ 1.
   - Estimación Monte Carlo: p̂ = A / N
   - Margen de error del 95%: M = z · √( p̂(1 − p̂) / N ). El error baja como 1/√N, así que con más años simulados la estimación es más precisa.

## Interpretación de la gráfica

La línea muestra la probabilidad exacta P de cada grupo de edad. Los puntos son la estimación p̂ de Monte Carlo y sus barras el margen M. El punto naranja marca el grupo de la edad elegida.

El resultado es un promedio del grupo, no una predicción para una persona. Los datos son de Francia y solo se modela la frecuencia de siniestros, no su costo.
