# Modelo OLG con reforma pensional endogena

El modulo implementa dos arquitecturas institucionales:

1. `regime = "compete"`: Ley 100, con eleccion entre informalidad,
   capitalizacion formal y PAYG formal;
2. `regime = "pillars"`: Ley 2381, con eleccion entre informalidad y un
   paquete formal obligatorio. El paquete asigna PAYG hasta 2.3 SMLMV y ahorro
   individual sobre el excedente.

La cuenta individual recibe 13.2% del excedente, la formula PAYG usa
`0.655 - 0.005 * s`, el piso contributivo es un SMLMV y la renta solidaria se
normaliza con la linea de pobreza extrema de 2023. El multiplicador PAYG de 1
es la normalizacion central; 2.9 es un estres no identificado.

## Pruebas

```powershell
Rscript tests/run_endogenous_tests.R
Rscript tests/run_policy_reform_tests.R
```

Las pruebas exigen consumo positivo, ahorro no negativo, participaciones
normalizadas, masas positivas en las tres categorias, cierre presupuestal y
residuos bajo tolerancia. La transicion se valida si converge al estado
estacionario o a un ciclo de dos generaciones no degenerado. En la calibracion
central, la brecha de repeticion a dos periodos es menor a 1%, aunque el estado
estacionario de pilares no es el limite dinamico.

## Analisis reproducible

```powershell
Rscript scripts/run_endogenous_analysis.R
Rscript scripts/run_policy_reform_analysis.R
```

El segundo script crea las trayectorias de capital, salarios, impuestos,
formalidad y componentes pensionales; microdatos por tipo; Gini y cuantiles de
salarios; equivalentes de consumo y participaciones de ganadores; una
comparacion estacionaria de estres; y las figuras del paper. Los ultimos tres
periodos se usan como buffer terminal y no se reportan como trayectoria.

## Alcance

Un periodo representa aproximadamente 40 anos. El modelo no contiene edades,
sexo, historias semanales de cotizacion ni un pilar semicontributivo separado.
La simulacion es contrafactual: no debe leerse como trayectoria anual ni como
evaluacion de una reforma ya implementada.
