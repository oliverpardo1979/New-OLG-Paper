# Eleccion endogena de regimen pensional

Este modulo amplia la implementacion validada sin eliminar los benchmarks
anteriores. Cada tipo compara simultaneamente tres alternativas:

1. trabajo informal;
2. trabajo formal con cuenta de capitalizacion;
3. trabajo formal con PAYG.

La capitalizacion es portable. El PAYG tiene una probabilidad de elegibilidad
creciente con el tipo y devuelve la cuenta nocional sin interes real cuando el
trabajador no cumple los requisitos. La calibracion usa un aporte total de 16%,
11.5% acreditado a capitalizacion y 13% a la cuenta nocional PAYG. Dos momentos
comparables del CEDE Pension Model disciplinan formalidad y elegibilidad; la
afiliacion PAYG de entrantes y el retorno real quedan como validaciones externas.
No se imponen cuotas, umbrales de eleccion ni shocks logit.

## Ejecutar las pruebas

```powershell
Rscript tests/run_endogenous_tests.R
```

Las pruebas exigen:

- convergencia y residuos inferiores a la tolerancia;
- consumo positivo y ahorro no negativo;
- participaciones positivas en las tres alternativas;
- dos umbrales interiores;
- el orden informal, capitalizacion y PAYG;
- consistencia de la identidad del ingreso y del presupuesto publico.

## Reproducir los resultados

```powershell
Rscript scripts/run_endogenous_analysis.R
```

El script resuelve primero una malla de 101 tipos y usa esa solucion para
inicializar la malla de 501 tipos. Tambien ejecuta una sensibilidad local de la
elegibilidad PAYG. Los resultados se guardan en `../results/`. Los archivos
`becerra2026_moment_comparison.csv` y
`becerra2026_parameter_provenance.csv` documentan, respectivamente, los ajustes
contra los datos y la condicion de cada parametro (directo, calibrado, validado
o retenido).

## Alcance

El modulo implementa el nuevo equilibrio estacionario. Las funciones de
transicion conservadas en `R/transition.R` corresponden a los benchmarks
exogenos anteriores y no deben usarse para una reforma de eleccion endogena.
