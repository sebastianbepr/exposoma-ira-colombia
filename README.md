# Piloto exposómico — determinantes ambientales de la transmisión de virus respiratorios (Colombia, 2022)

Código y datos del análisis de la estructura de la matriz de exposición ambiental
que acompaña al manuscrito *Determinantes ambientales de la transmisión de virus
respiratorios* (Universidad Militar Nueva Granada).

**Autor:** Sebastián A. Beltrán Prieto, MD, MSc (Bioinformática)

## Contenido

```
exposoma-ira-colombia/
├── Script.R                              # análisis completo (R)
├── data/
│   └── BASE_DE_DATOS_ORGANIZADA.xlsx     # matriz climática diaria, hoja "CLIMA"
├── .gitignore
├── LICENSE
└── README.md
```

## Datos

`data/BASE_DE_DATOS_ORGANIZADA.xlsx`, hoja `CLIMA`. Formato ancho con encabezado
en dos filas: la primera identifica el departamento (un bloque de columnas por
departamento, 23 en total) y la segunda las variables dentro de cada bloque.

Cobertura: 482 días de 2022, 23 unidades territoriales (22 departamentos + Bogotá).

Variables por bloque: día, temperatura seca (promedio, mínima, máxima),
temperatura húmeda, humedad relativa (máxima, mínima, promedio), radiación UV,
precipitación diaria y altura.

Variables derivadas en el script:

| Variable | Definición |
|---|---|
| `ATD` | amplitud térmica diurna (`T_max − T_min`) |
| `es` | presión de vapor de saturación (hPa), fórmula de Bolton |
| `e` | presión de vapor (hPa) = `es · HR/100` |
| `HA` | humedad absoluta (g/m³) = `216.7 · e / (T_med + 273.15)` |

## Análisis (`Script.R`)

1. Lectura del formato ancho con encabezado de dos filas y paso a formato largo.
2. Cálculo de variables derivadas (ATD, presión de vapor, humedad absoluta).
3. Diagnóstico de completitud por departamento y revisión de la escala de precipitación.
4. Conjunto analítico: departamentos con ≥ 120 observaciones completas.
5. Estructura de correlación entre `T_med`, `ATD`, `HR` y `HA` (global y por departamento).
6. ACP sobre la matriz de exposición diaria.
7. ACP sobre las medias departamentales + k-medias (k = 3) para una tipología territorial.
8. Estabilidad de las cargas del ACP entre temporada de lluvias y temporada seca.

## Requisitos

R ≥ 4.1 y los paquetes:

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "factoextra", "corrplot"))
```

`factoextra` y `corrplot` solo hacen falta para las visualizaciones (comentadas en el script).

## Cómo ejecutarlo

El script usa rutas relativas, así que debe correrse con el directorio de trabajo
en la raíz del repositorio:

```r
setwd("ruta/a/exposoma-ira-colombia")
source("Script.R")
```

O bien abrir el proyecto en RStudio desde esta carpeta.

## Cita

Si usa este código o los datos, por favor cite el manuscrito asociado.

## Licencia

Código bajo licencia MIT (ver `LICENSE`). Los datos climáticos derivan de fuentes
del IDEAM; verifique sus condiciones de uso antes de redistribuirlos.
