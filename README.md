# Exposición ambiental colombiana como base para el desarrollo de un modelo exposómico de alerta temprana

Código, datos y salidas del análisis de la estructura de la matriz de exposición ambiental
(Colombia, 2022) que acompaña al manuscrito del mismo título, en evaluación en la
*Revista Panamericana de Salud Pública*.

**Autor:** Sebastián A. Beltrán Prieto, MD, MSc (Bioinformática)

## Contenido

```
exposoma-ira-colombia/
├── Script_04.R                           # análisis completo (R), versión 04 del 21/10/2026
├── data/
│   └── BASE_DE_DATOS_ORGANIZADA.xlsx     # matriz climática diaria, hoja "CLIMA"
├── outputs/                              # tablas y figuras generadas por Script_04.R
├── .here                                 # marca la raíz del proyecto para here::here()
├── LICENSE
└── README.md
```

## Datos

`data/BASE_DE_DATOS_ORGANIZADA.xlsx`, hoja `CLIMA`, con registros del IDEAM. Formato ancho con
encabezado en dos filas: la primera identifica el departamento (un bloque de columnas por
departamento) y la segunda las variables dentro de cada bloque.

- Cobertura: 348 días de 2022 (2 de enero a 15 de diciembre), 23 departamentos.
- Panel analítico: departamentos con al menos 120 días con registro completo (temperatura
  mínima, media y máxima y humedad relativa) → 14 departamentos, 3 154 días-departamento
  (`outputs/completitud.csv`).
- Algunos encabezados de departamento conservan notas de trabajo del archivo original
  (p. ej., "CHOCO ORGANIZAR TEMPERATURAS"). Se mantienen sin modificar para no alterar los
  datos; ninguno de esos departamentos entra en el panel analítico.

Variables de exposición analizadas: temperatura media (`T_med`), amplitud térmica diurna
(`ATD`), humedad relativa (`HR`) y humedad absoluta (`HA`, derivada en el script).

## Análisis (`Script_04.R`)

El script está organizado en bloques marcados `__INDEX_xx__`:

| Bloque | Contenido | Salida relacionada en `outputs/` |
|---|---|---|
| 00 | Paquetes, rutas y parámetros (semilla `20220101`) | — |
| 01 | Lectura del libro ancho con encabezado de dos filas | — |
| 02 | Control de calidad: rangos de plausibilidad y coherencia interna | `calidad.csv` |
| 03 | Variables derivadas: ATD y humedad absoluta | — |
| 04 | Panel analítico: completitud y racha máxima | `completitud.csv` |
| 05 | Correlaciones global, intra- y entre departamentos | `correlacion_global.csv`, `correlacion_intra.csv`, `correlacion_entre.csv`, `humedad_por_departamento.csv` |
| 06 | ACP diaria: Kaiser, sedimentación y análisis paralelo de Horn | `retencion.csv`, `cargas.csv`, `figura_S1_sedimentacion.png` |
| 07 | Sensibilidad: sin humedad absoluta e intradepartamental | `sensibilidad.csv`, `congruencia_intra.csv` |
| 08 | Estabilidad de cargas: bootstrap por conglomerados (congruencia de Tucker) | `bootstrap.csv` |
| 09 | Estabilidad temporal: ACP mensual y por temporada frente a la anual | `estabilidad.csv` |
| 10 | ACP de medias departamentales y tipología (Ward, k-medias, silueta) | `tipologia.csv`, `silueta.csv` |
| 11 | Figuras | `figura_1.png`, `figura_2.png` (y versiones `.eps`) |
| 12 | Exportación de tablas e información de sesión | `sessionInfo.txt` |
| 13 | Resumen final en consola | — |

Las figuras `figura_1.png` y `figura_2.png` corresponden a las Figuras 1 y 2 del manuscrito, y
`figura_S1_sedimentacion.png` al gráfico de sedimentación suplementario. La carpeta
`outputs/` contiene las salidas de la ejecución usada para el manuscrito; al correr el
script se sobrescriben.

## Requisitos

Ejecutado con R 4.6.0 (ver `outputs/sessionInfo.txt`). Paquetes y versiones usadas:

| Paquete | Versión |
|---|---|
| here | 1.0.2 |
| readxl | 1.5.0 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| purrr | 1.2.2 |
| ggplot2 | 4.0.3 |
| cluster | 2.1.8.2 |
| ggrepel | 0.9.8 |
| scales | 1.4.0 |

El script también llama funciones de `stringi` y `tibble` con `::`.

```r
install.packages(c("here", "readxl", "dplyr", "tidyr", "purrr", "ggplot2",
                   "cluster", "ggrepel", "scales", "stringi", "tibble"))
```

## Cómo ejecutarlo

Las rutas se resuelven desde la raíz del proyecto con `here::here()` (el archivo `.here`
marca esa raíz), sin `setwd()`. Basta con abrir R en la carpeta del repositorio, o en
cualquier subcarpeta, y ejecutar:

```r
source("Script_04.R")
```

## Uso de inteligencia artificial

Como se declara en el manuscrito, Claude (Anthropic) se usó para la revisión y la
reorganización del código de análisis en R. El código fue ejecutado por el autor sobre sus
datos y todas las salidas fueron verificadas por él.

## Cita

Si usa este código o los datos, por favor cite el manuscrito asociado.

## Licencia

Código bajo licencia MIT (ver `LICENSE`). Los datos climáticos derivan de fuentes del IDEAM;
verifique sus condiciones de uso antes de redistribuirlos.
