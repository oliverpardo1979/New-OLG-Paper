# Eleccion endogena de regimen pensional

Este modulo amplia la implementacion validada sin eliminar los benchmarks
anteriores. Cada tipo compara simultaneamente tres alternativas:

1. trabajo informal;
2. trabajo formal con cuenta de capitalizacion;
3. trabajo formal con PAYG.

La capitalizacion es portable. El PAYG tiene una probabilidad de elegibilidad
creciente con el tipo, devuelve una fraccion del saldo cuando el trabajador no
cumple los requisitos y aplica un beneficio con piso y techo. No se imponen
cuotas, umbrales de eleccion ni shocks logit.

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
elegibilidad PAYG. Los resultados se guardan en `outputs/`.

## Alcance

El modulo implementa el nuevo equilibrio estacionario. Las funciones de
transicion conservadas en `R/transition.R` corresponden a los benchmarks
exogenos anteriores y no deben usarse para una reforma de eleccion endogena.
