# =============================================================================
# Piloto exposómico — estructura de la matriz ambiental (Colombia, 2022)
# Manuscrito: determinantes ambientales de la transmisión de virus respiratorios
# Autor: Sebastian A. Beltran Prieto MD, MSc (Bioinformatica)
# Fecha: 10/09/2026
# =============================================================================

# ---- 0. Configuración inicial y paquetes -------------------------------------
# Descomentar para instalar los paquetes necesarios si no están disponibles:
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "factoextra", "corrplot"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Uso de ruta relativa. 

ruta <- "data/BASE_DE_DATOS_ORGANIZADA.xlsx"  

# ---- 1. Lectura de formato ancho con encabezado en dos filas -----------------
crudo <- read_excel(ruta, sheet = "CLIMA", col_names = FALSE)
fila_dep <- as.character(unlist(crudo[1, ]))
fila_var <- as.character(unlist(crudo[2, ]))

# El nombre del departamento sólo aparece en la primera columna de cada bloque
dep <- fila_dep
for (j in seq_along(dep)) {
  if (is.na(dep[j]) || dep[j] == "NA" || trimws(dep[j]) == "") {
    dep[j] <- if (j > 1) dep[j - 1] else NA_character_
  }
}

cuerpo <- crudo[-c(1, 2), ]
fecha  <- as.Date(as.numeric(unlist(cuerpo[, 1])), origin = "1899-12-30")

norm_var <- function(v) {
  v <- toupper(trimws(v))
  v[v %in% c("PRECIPITACIONES DIARIA", "PRECIPITACIONES")] <- "PRECIPITACION"
  v
}
fila_var <- norm_var(fila_var)

largo <- lapply(unique(na.omit(dep)), function(d) {
  idx <- which(dep == d)
  sub <- cuerpo[, idx, drop = FALSE]
  names(sub) <- make.unique(fila_var[idx])
  # Limpieza de columnas duplicadas
  sub <- sub[, !duplicated(names(sub)), drop = FALSE]
  sub$DEP <- d
  sub$FECHA <- fecha
  sub
}) %>% bind_rows()

num <- c("TEMP SECA PROMEDIO", "TEMP SECA MINIMA", "TEMP SECA MAXIMA",
         "HUMEDAD RELATIVA PROMEDIO", "PRECIPITACION", "RADIACION UV", "ALTURA")
for (v in intersect(num, names(largo))) {
  largo[[v]] <- suppressWarnings(as.numeric(largo[[v]]))
}

# ---- 2. Variables derivadas --------------------------------------------------
d <- largo %>%
  rename(T_med = `TEMP SECA PROMEDIO`,
         T_min = `TEMP SECA MINIMA`,
         T_max = `TEMP SECA MAXIMA`,
         HR    = `HUMEDAD RELATIVA PROMEDIO`) %>%
  mutate(
    ATD = T_max - T_min,                                   # amplitud térmica diurna
    es  = 6.112 * exp(17.67 * T_med / (T_med + 243.5)),    # presión de vapor sat. (hPa)
    e   = es * HR / 100,                                   # presión de vapor (hPa)
    HA  = 216.7 * e / (T_med + 273.15)                     # humedad absoluta (g/m3)
  )

# ---- 3. Diagnóstico de completitud y de valores implausibles -----------------
completitud <- d %>%
  group_by(DEP) %>%
  summarise(across(c(T_med, ATD, HR, PRECIPITACION), ~ sum(!is.na(.x))), n = n())
print(as.data.frame(completitud))

# Precipitación: revisar la escala antes de usarla
d %>% 
  filter(PRECIPITACION > 100) %>%
  summarise(n = n(), min = min(PRECIPITACION), max = max(PRECIPITACION)) %>% 
  print()

# ---- 4. Conjunto analítico ---------------------------------------------------
v  <- c("T_med", "ATD", "HR", "HA")
ok <- d %>% filter(if_all(all_of(v), ~ !is.na(.x)))
deps_ok <- ok %>% count(DEP) %>% filter(n >= 120) %>% pull(DEP)
ok <- ok %>% filter(DEP %in% deps_ok)
cat("Departamentos:", length(deps_ok), " Observaciones:", nrow(ok), "\n")

# ---- 5. Estructura de correlación --------------------------------------------
M <- cor(ok[, v], use = "complete.obs")
print(round(M, 2))
# Descomentar para visualizar
# corrplot::corrplot(M, method = "number", type = "upper")

# HR frente a HA: global y dentro de cada departamento
cat("r(HR,HA) agrupando todos los departamentos:",
    round(cor(ok$HR, ok$HA), 2), "\n")
ok %>% group_by(DEP) %>%
  summarise(r_HR_HA = round(cor(HR, HA), 2)) %>%
  arrange(r_HR_HA) %>% as.data.frame() %>% print()

# ---- 6. ACP sobre la matriz de exposición diaria ------------------------------
acp <- prcomp(ok[, v], scale. = TRUE)
print(summary(acp))
print(round(acp$rotation[, 1:2], 2))
# Descomentar para visualizar
# factoextra::fviz_pca_var(acp, repel = TRUE)

# ---- 7. ACP entre departamentos (configuraciones territoriales) ---------------
medias <- ok %>% group_by(DEP) %>% summarise(across(all_of(v), mean))
acp_d  <- prcomp(medias[, v], scale. = TRUE)
print(summary(acp_d))
biplot(acp_d, xlabs = medias$DEP)

# Tipología territorial: k-medias sobre los componentes retenidos
set.seed(1)
km <- kmeans(scale(medias[, v]), centers = 3, nstart = 50)
print(data.frame(DEP = medias$DEP, grupo = km$cluster))

# ---- 8. Estabilidad de la estructura entre temporadas -------------------------
# Régimen bimodal: comparar las cargas del ACP entre trimestres húmedos y secos
ok <- ok %>% mutate(
  mes = as.integer(format(FECHA, "%m")),
  temporada = case_when(
    mes %in% c(4,5,10,11) ~ "lluvias",
    mes %in% c(1,2,7,8)   ~ "seca",
    TRUE                  ~ "transicion"
  )
)

for (s in c("lluvias", "seca")) {
  sub <- ok %>% filter(temporada == s)
  a <- prcomp(sub[, v], scale. = TRUE)
  cat("\n---", s, "— varianza PC1-PC2:",
      round(summary(a)$importance[2, 1:2], 3), "\n")
  print(round(a$rotation[, 1:2], 2))
}