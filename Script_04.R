# =================================================================================================
# Exposomic pilot: structure of the environmental exposure matrix (Colombia, 2022)
# Manuscript: Exposición ambiental colombiana como base para el desarrollo de un modelo exposómico de alerta temprana
# Author: Sebastian A. Beltran Prieto MD, MSc (Bioinformatics)
# Version: 04 (RPSP figures: no panel letters, accented names, decimal comma) - 21/10/2026
# =================================================================================================
#
# INDEX
#   __INDEX_00__  Setup: packages, paths, analysis parameters
#   __INDEX_01__  Reading the wide workbook with a two-row header
#   __INDEX_02__  Quality control: plausibility ranges and internal consistency
#   __INDEX_03__  Derived variables: diurnal temperature range and absolute humidity
#   __INDEX_04__  Analytic panel: completeness and longest uninterrupted run
#   __INDEX_05__  Correlation structure: pooled, within- and between-department
#   __INDEX_06__  PCA on daily data: retention criteria (Kaiser, scree, Horn parallel analysis)
#   __INDEX_07__  PCA sensitivity: without derived humidity and within-department (centred)
#   __INDEX_08__  Loading stability: cluster bootstrap by department (Tucker congruence)
#   __INDEX_09__  Temporal stability: monthly PCA versus the annual solution
#   __INDEX_10__  Between-department PCA and territorial typology (Ward, k-means, silhouette)
#   __INDEX_11__  Figures (Figure 2a, Figure 2b, supplementary scree plot)
#   __INDEX_12__  Export of tables and session information
#   __INDEX_13__  Final summary
# =================================================================================================


# __INDEX_00__ Setup ------------------------------------------------------------------------------
# All paths are resolved from the project root with here::here(); no setwd() is used anywhere.

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(cluster)
  library(ggrepel)
})

clima_ruta_xlsx   <- here::here("data", "BASE_DE_DATOS_ORGANIZADA.xlsx")
resultados_dir    <- here::here("outputs")
dir.create(resultados_dir, showWarnings = FALSE, recursive = TRUE)

exposicion_vars_chr    <- c("T_med", "ATD", "HR", "HA")
dias_minimos_int       <- 120L     # minimum complete days per department
bootstrap_reps_int     <- 1000L    # cluster-bootstrap replicates
paralelo_reps_int      <- 500L     # Horn parallel-analysis replicates
componentes_k_int      <- 2L       # components compared in stability analyses

set.seed(20220101)

# Tucker's congruence coefficient between two loading vectors. Values >= 0.95 are conventionally
# read as factor equality and 0.85-0.94 as fair similarity (Lorenzo-Seva & ten Berge, 2006).
congruencia_tucker <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

# Orthogonal Procrustes alignment of a loading matrix to a reference. It removes the sign and
# order indeterminacy of PCA solutions before loadings from different samples are compared.
alinear_procrustes <- function(cargas_mat, referencia_mat) {
  svd_lst <- svd(t(cargas_mat) %*% referencia_mat)
  cargas_mat %*% (svd_lst$u %*% t(svd_lst$v))
}

# Loadings expressed as variable-component correlations (eigenvector x sqrt(eigenvalue)), which
# are the quantities that should be reported in the manuscript.
cargas_correlacion <- function(acp_prcomp, k = componentes_k_int) {
  acp_prcomp$rotation[, seq_len(k), drop = FALSE] %*% diag(acp_prcomp$sdev[seq_len(k)])
}


# __INDEX_01__ Reading the wide workbook --------------------------------------------------------
# The sheet stores one block of columns per department. Row 1 carries the department name only in
# the first column of each block; row 2 carries the variable names; column 1 carries the date.

message("[01] Reading workbook: ", clima_ruta_xlsx)

clima_crudo_tbl <- read_excel(clima_ruta_xlsx, sheet = "CLIMA", col_names = FALSE,
                              .name_repair = "minimal")

clima_departamento_chr <- as.character(unlist(clima_crudo_tbl[1, ]))
clima_variable_chr     <- as.character(unlist(clima_crudo_tbl[2, ]))

# Forward-fill department names across each block; column 1 (date) is excluded explicitly so
# that a label such as "FECHA" is never treated as a department.
clima_departamento_chr[clima_departamento_chr %in% c("NA", "") ] <- NA_character_
clima_departamento_chr[1] <- NA_character_
clima_departamento_chr <- tidyr::fill(tibble(d = clima_departamento_chr), d)$d
clima_departamento_chr[1] <- NA_character_

clima_variable_chr <- toupper(trimws(clima_variable_chr))
clima_variable_chr[clima_variable_chr %in% c("PRECIPITACIONES DIARIA", "PRECIPITACIONES")] <-
  "PRECIPITACION"

clima_cuerpo_tbl <- clima_crudo_tbl[-c(1, 2), ]
clima_fecha_vec  <- as.Date(as.numeric(unlist(clima_cuerpo_tbl[, 1])), origin = "1899-12-30")

clima_largo_df <- unique(na.omit(clima_departamento_chr)) %>%
  map(function(departamento) {
    idx_int <- which(clima_departamento_chr == departamento)
    bloque_tbl <- clima_cuerpo_tbl[, idx_int, drop = FALSE]
    nombres_chr <- clima_variable_chr[idx_int]
    if (any(duplicated(nombres_chr))) {
      warning("Duplicated variable names in block ", departamento, ": ",
              paste(unique(nombres_chr[duplicated(nombres_chr)]), collapse = ", "),
              ". The first occurrence is kept.")
    }
    bloque_tbl <- bloque_tbl[, !duplicated(nombres_chr), drop = FALSE]
    names(bloque_tbl) <- nombres_chr[!duplicated(nombres_chr)]
    bloque_tbl %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))),
             DEP = departamento, FECHA = clima_fecha_vec)
  }) %>%
  bind_rows() %>%
  filter(!is.na(FECHA))

message("[01] Departments read: ", n_distinct(clima_largo_df$DEP),
        " | Date range: ", paste(range(clima_largo_df$FECHA), collapse = " to "))


# __INDEX_02__ Quality control ------------------------------------------------------------------
# Implausible values are set to NA rather than silently entering the correlation matrix. Ranges are
# deliberately wide so that genuine paramo and Caribbean extremes are retained.

message("[02] Quality control")

clima_curado_df <- clima_largo_df %>%
  rename(T_med = `TEMP SECA PROMEDIO`, T_min = `TEMP SECA MINIMA`,
         T_max = `TEMP SECA MAXIMA`,   HR    = `HUMEDAD RELATIVA PROMEDIO`) %>%
  mutate(
    qc_temperatura = !is.na(T_med) & (T_med < -5 | T_med > 40),
    qc_humedad     = !is.na(HR) & (HR <= 0 | HR > 100),
    qc_orden       = !is.na(T_max) & !is.na(T_min) &
                     (T_max < T_min | T_med < T_min - 0.5 | T_med > T_max + 0.5),
    T_med = if_else(qc_temperatura, NA_real_, T_med),
    HR    = if_else(qc_humedad,     NA_real_, HR),
    T_min = if_else(qc_orden,       NA_real_, T_min),
    T_max = if_else(qc_orden,       NA_real_, T_max)
  )

calidad_resumen_df <- clima_curado_df %>%
  group_by(DEP) %>%
  summarise(dias = n(),
            temperatura_fuera_rango = sum(qc_temperatura),
            humedad_fuera_rango     = sum(qc_humedad),
            temperatura_incoherente = sum(qc_orden),
            .groups = "drop")
print(as.data.frame(calidad_resumen_df))


# __INDEX_03__ Derived variables ------------------------------------------------------------------
# Absolute humidity is computed from mean temperature and mean relative humidity (Magnus-Bolton
# saturation vapour pressure). Because HA is a deterministic function of T_med and HR, the four
# variables cannot carry four independent dimensions; INDEX_07 quantifies this explicitly.

clima_curado_df <- clima_curado_df %>%
  mutate(
    ATD = T_max - T_min,
    es  = 6.112 * exp(17.67 * T_med / (T_med + 243.5)),
    e   = es * HR / 100,
    HA  = 216.7 * e / (T_med + 273.15)
  )


# __INDEX_04__ Analytic panel ------------------------------------------------------------------
# The manuscript previously stated "120 dias de registro ininterrumpido". The inclusion rule is in
# fact >= 120 complete days; the longest uninterrupted run is reported so the text can be exact.

message("[04] Building the analytic panel")

racha_maxima <- function(fechas_vec) {
  if (length(fechas_vec) == 0) return(0L)
  r <- rle(c(1, diff(sort(fechas_vec))) == 1)
  max(c(0L, r$lengths[r$values]), na.rm = TRUE) + 1L
}

completitud_df <- clima_curado_df %>%
  filter(if_all(all_of(exposicion_vars_chr), ~ !is.na(.x))) %>%
  group_by(DEP) %>%
  summarise(dias_completos = n(), racha_maxima = racha_maxima(FECHA), .groups = "drop") %>%
  arrange(desc(dias_completos))
print(as.data.frame(completitud_df))

departamentos_incluidos_chr <- completitud_df %>%
  filter(dias_completos >= dias_minimos_int) %>% pull(DEP)

exposicion_diaria_df <- clima_curado_df %>%
  filter(DEP %in% departamentos_incluidos_chr,
         if_all(all_of(exposicion_vars_chr), ~ !is.na(.x))) %>%
  select(DEP, FECHA, all_of(exposicion_vars_chr))

message("[04] Departments: ", length(departamentos_incluidos_chr),
        " | Department-days: ", nrow(exposicion_diaria_df))


# __INDEX_05__ Correlation structure --------------------------------------------------------------
# A pooled correlation mixes two different quantities: day-to-day covariation inside departments
# (within) and differences between departmental climates (between). Only the decomposition shows
# which of the two drives the low pooled HR-HA correlation reported in the manuscript.

message("[05] Correlation decomposition")

exposicion_centrada_df <- exposicion_diaria_df %>%
  group_by(DEP) %>%
  mutate(across(all_of(exposicion_vars_chr), ~ .x - mean(.x))) %>%
  ungroup()

exposicion_medias_df <- exposicion_diaria_df %>%
  group_by(DEP) %>%
  summarise(across(all_of(exposicion_vars_chr), mean), .groups = "drop")

correlacion_global_mat <- cor(exposicion_diaria_df[, exposicion_vars_chr])
correlacion_intra_mat  <- cor(exposicion_centrada_df[, exposicion_vars_chr])
correlacion_entre_mat  <- cor(exposicion_medias_df[, exposicion_vars_chr])

message("Pooled correlation matrix");         print(round(correlacion_global_mat, 2))
message("Within-department correlation");      print(round(correlacion_intra_mat, 2))
message("Between-department correlation (n = ", nrow(exposicion_medias_df), ")")
print(round(correlacion_entre_mat, 2))

humedad_correlacion_dep_df <- exposicion_diaria_df %>%
  group_by(DEP) %>%
  summarise(n = n(),
            r_HR_HA   = cor(HR, HA),
            r_Tmed_HA = cor(T_med, HA),
            .groups = "drop") %>%
  arrange(r_HR_HA)
print(as.data.frame(humedad_correlacion_dep_df), digits = 2)

humedad_resumen_lst <- list(
  r_global = correlacion_global_mat["HR", "HA"],
  r_intra_mediana = median(humedad_correlacion_dep_df$r_HR_HA),
  r_intra_rango = range(humedad_correlacion_dep_df$r_HR_HA),
  r_intra_centrado = correlacion_intra_mat["HR", "HA"],
  r_entre = correlacion_entre_mat["HR", "HA"]
)


# __INDEX_06__ PCA on daily data and retention criteria ------------------------------------------
# Three convergent criteria justify the number of retained components: Kaiser (eigenvalue > 1),
# the scree elbow, and Horn's parallel analysis (observed eigenvalue above the 95th percentile of
# eigenvalues from uncorrelated data of identical dimensions). This answers reviewer point 2.

message("[06] Daily PCA and retention criteria")

exposicion_acp_prcomp <- prcomp(exposicion_diaria_df[, exposicion_vars_chr], scale. = TRUE)

paralelo_autovalores_mat <- replicate(paralelo_reps_int, {
  aleatorio_mat <- apply(as.matrix(exposicion_diaria_df[, exposicion_vars_chr]), 2, sample)
  eigen(cor(aleatorio_mat), only.values = TRUE)$values
})

retencion_df <- tibble(
  componente         = paste0("PC", seq_along(exposicion_acp_prcomp$sdev)),
  autovalor          = exposicion_acp_prcomp$sdev^2,
  varianza_pct       = 100 * autovalor / sum(autovalor),
  varianza_acum_pct  = cumsum(varianza_pct),
  paralelo_p95       = apply(paralelo_autovalores_mat, 1, quantile, probs = 0.95),
  kaiser_retener     = autovalor > 1,
  paralelo_retener   = autovalor > paralelo_p95
)
print(as.data.frame(retencion_df), digits = 3)

exposicion_cargas_mat <- cargas_correlacion(exposicion_acp_prcomp)
colnames(exposicion_cargas_mat) <- paste0("PC", seq_len(componentes_k_int))
message("Loadings (eigenvectors)");                print(round(exposicion_acp_prcomp$rotation, 2))
message("Loadings (variable-component r)");       print(round(exposicion_cargas_mat, 2))


# __INDEX_07__ PCA sensitivity analyses -----------------------------------------------------------
# (a) Without HA: if two components still dominate, the reduced dimensionality is not an artefact
#     of including a variable derived from T_med and HR.
# (b) Within-department PCA: variables centred and scaled inside each department, so the solution
#     describes day-to-day covariation only. This is the structure a time-series exposure model
#     actually faces; the pooled PCA is dominated by altitude-driven between-department contrasts.

message("[07] Sensitivity PCAs")

resumen_acp <- function(acp_prcomp, etiqueta) {
  tibble(analisis = etiqueta,
         componente = paste0("PC", seq_along(acp_prcomp$sdev)),
         autovalor = acp_prcomp$sdev^2,
         varianza_pct = 100 * acp_prcomp$sdev^2 / sum(acp_prcomp$sdev^2))
}

sin_ha_vars_chr <- setdiff(exposicion_vars_chr, "HA")
sin_ha_acp_prcomp <- prcomp(exposicion_diaria_df[, sin_ha_vars_chr], scale. = TRUE)

exposicion_intra_df <- exposicion_diaria_df %>%
  group_by(DEP) %>%
  mutate(across(all_of(exposicion_vars_chr), ~ as.numeric(scale(.x)))) %>%
  ungroup()
intra_acp_prcomp <- prcomp(exposicion_intra_df[, exposicion_vars_chr], scale. = FALSE)

intra_cargas_mat <- alinear_procrustes(cargas_correlacion(intra_acp_prcomp),
                                       exposicion_cargas_mat)

sensibilidad_df <- bind_rows(
  resumen_acp(exposicion_acp_prcomp, "principal_4_variables"),
  resumen_acp(sin_ha_acp_prcomp,     "sin_humedad_absoluta"),
  resumen_acp(intra_acp_prcomp,      "intra_departamento")
)
print(as.data.frame(sensibilidad_df), digits = 3)

message("Within-department loadings (Procrustes-aligned to the pooled solution)")
print(round(intra_cargas_mat, 2))
intra_congruencia_vec <- map_dbl(seq_len(componentes_k_int), ~
  congruencia_tucker(intra_cargas_mat[, .x], exposicion_cargas_mat[, .x]))
message("Tucker congruence, within vs pooled: ",
        paste(round(intra_congruencia_vec, 3), collapse = " / "))


# __INDEX_08__ Cluster bootstrap by department --------------------------------------------------
# Days are autocorrelated and nested in departments, so resampling days would overstate
# precision. Whole departments are resampled with replacement, which respects both dependencies.

message("[08] Cluster bootstrap (", bootstrap_reps_int, " replicates)")

bootstrap_df <- map_dfr(seq_len(bootstrap_reps_int), function(b) {
  muestra_chr <- sample(departamentos_incluidos_chr, replace = TRUE)
  muestra_df  <- map_dfr(muestra_chr, ~ filter(exposicion_diaria_df, DEP == .x))
  acp_b <- prcomp(muestra_df[, exposicion_vars_chr], scale. = TRUE)
  cargas_b <- alinear_procrustes(cargas_correlacion(acp_b), exposicion_cargas_mat)
  tibble(replica = b,
         var_pc1 = 100 * acp_b$sdev[1]^2 / sum(acp_b$sdev^2),
         var_pc12 = 100 * sum(acp_b$sdev[1:2]^2) / sum(acp_b$sdev^2),
         autovalor_pc2 = acp_b$sdev[2]^2,
         autovalor_pc3 = acp_b$sdev[3]^2,
         cong_pc1 = congruencia_tucker(cargas_b[, 1], exposicion_cargas_mat[, 1]),
         cong_pc2 = congruencia_tucker(cargas_b[, 2], exposicion_cargas_mat[, 2]))
})

bootstrap_resumen_df <- bootstrap_df %>%
  select(-replica) %>%
  pivot_longer(everything(), names_to = "indicador") %>%
  group_by(indicador) %>%
  summarise(mediana = median(value),
            ic_inf = quantile(value, 0.025),
            ic_sup = quantile(value, 0.975),
            .groups = "drop")
print(as.data.frame(bootstrap_resumen_df), digits = 3)


# __INDEX_09__ Temporal stability ---------------------------------------------------------------
# The earlier rainy/dry partition used the Andean bimodal calendar for every department, although
# Orinoquia and Amazonia follow a monomodal regime. A regime-agnostic alternative is to fit one PCA
# per calendar month and compare each against the annual solution. The original partition is kept
# only for continuity with the previous draft.

message("[09] Temporal stability")

exposicion_diaria_df <- exposicion_diaria_df %>%
  mutate(mes = as.integer(format(FECHA, "%m")),
         temporada = case_when(mes %in% c(4, 5, 10, 11) ~ "lluvias_andina",
                               mes %in% c(1, 2, 7, 8)   ~ "seca_andina",
                               TRUE                     ~ "transicion"))

estabilidad_subconjunto <- function(sub_df, etiqueta) {
  acp_s <- prcomp(sub_df[, exposicion_vars_chr], scale. = TRUE)
  cargas_s <- alinear_procrustes(cargas_correlacion(acp_s), exposicion_cargas_mat)
  tibble(periodo = etiqueta,
         n = nrow(sub_df),
         var_pc1 = 100 * acp_s$sdev[1]^2 / sum(acp_s$sdev^2),
         var_pc2 = 100 * acp_s$sdev[2]^2 / sum(acp_s$sdev^2),
         cong_pc1 = congruencia_tucker(cargas_s[, 1], exposicion_cargas_mat[, 1]),
         cong_pc2 = congruencia_tucker(cargas_s[, 2], exposicion_cargas_mat[, 2]))
}

estabilidad_mensual_df <- exposicion_diaria_df %>%
  group_split(mes) %>%
  map_dfr(~ estabilidad_subconjunto(.x, sprintf("mes_%02d", unique(.x$mes))))

estabilidad_temporada_df <- exposicion_diaria_df %>%
  filter(temporada != "transicion") %>%
  group_split(temporada) %>%
  map_dfr(~ estabilidad_subconjunto(.x, unique(.x$temporada)))

print(as.data.frame(bind_rows(estabilidad_temporada_df, estabilidad_mensual_df)), digits = 3)


# __INDEX_10__ Between-department PCA and territorial typology ------------------------------------
# With 14 departments the typology is exploratory. The number of groups is chosen by average
# silhouette width, and agreement between Ward hierarchical clustering and k-means is reported
# instead of reading groups visually from the biplot.

message("[10] Between-department PCA and typology")

medias_acp_prcomp <- prcomp(exposicion_medias_df[, exposicion_vars_chr], scale. = TRUE)
medias_escaladas_mat <- scale(exposicion_medias_df[, exposicion_vars_chr])
rownames(medias_escaladas_mat) <- exposicion_medias_df$DEP
medias_distancia_dist <- dist(medias_escaladas_mat)
medias_ward_hclust <- hclust(medias_distancia_dist, method = "ward.D2")

silueta_df <- map_dfr(2:5, function(k) {
  grupos_ward <- cutree(medias_ward_hclust, k)
  grupos_km   <- kmeans(medias_escaladas_mat, centers = k, nstart = 100)$cluster
  tibble(k = k,
         silueta_ward   = mean(silhouette(grupos_ward, medias_distancia_dist)[, 3]),
         silueta_kmeans = mean(silhouette(grupos_km,   medias_distancia_dist)[, 3]))
})
print(as.data.frame(silueta_df), digits = 3)

grupos_k_int <- silueta_df$k[which.max(silueta_df$silueta_ward)]
tipologia_df <- tibble(
  DEP = exposicion_medias_df$DEP,
  grupo_ward = cutree(medias_ward_hclust, grupos_k_int),
  grupo_kmeans = kmeans(medias_escaladas_mat, centers = grupos_k_int, nstart = 100)$cluster
) %>%
  bind_cols(as_tibble(medias_acp_prcomp$x[, 1:2]))

tipologia_acuerdo_tbl <- table(tipologia_df$grupo_ward, tipologia_df$grupo_kmeans)
message("Selected k = ", grupos_k_int, " (maximum Ward silhouette)")
print(as.data.frame(tipologia_df), digits = 2)
print(tipologia_acuerdo_tbl)


# __INDEX_11__ Figures --------------------------------------------------------------------------
message("[11] Figures")

tema_articulo <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot")

# Spanish typography for the journal: decimal comma and accented department names
coma_decimal_fn <- scales::label_number(decimal.mark = ",", big.mark = " ")
nombres_departamento_chr <- c(
  ANTIOQUIA = "Antioquia", ATLANTICO = "Atl\u00e1ntico", BOGOTA = "Bogot\u00e1 D.C.",
  BOLIVAR = "Bol\u00edvar", BOYACA = "Boyac\u00e1", CALDAS = "Caldas", CAQUETA = "Caquet\u00e1",
  CASANARE = "Casanare", CAUCA = "Cauca", CESAR = "Cesar", CHOCO = "Choc\u00f3",
  CORDOBA = "C\u00f3rdoba", CUNDINAMARCA = "Cundinamarca", GUAINIA = "Guain\u00eda",
  GUAVIARE = "Guaviare", HUILA = "Huila", LA_GUAJIRA = "La Guajira", MAGDALENA = "Magdalena",
  META = "Meta", NARINO = "Nari\u00f1o", NORTE_DE_SANTANDER = "Norte de Santander",
  PUTUMAYO = "Putumayo", QUINDIO = "Quind\u00edo", RISARALDA = "Risaralda",
  SANTANDER = "Santander", SUCRE = "Sucre", TOLIMA = "Tolima", VALLE_DEL_CAUCA = "Valle del Cauca",
  VAUPES = "Vaup\u00e9s", VICHADA = "Vichada", AMAZONAS = "Amazonas", ARAUCA = "Arauca",
  SAN_ANDRES = "San Andr\u00e9s")
nombre_departamento_fn <- function(x) {
  clave_chr <- toupper(stringi::stri_trans_general(trimws(x), "Latin-ASCII"))
  clave_chr <- gsub("[^A-Z]+", "_", clave_chr)
  clave_chr <- sub("^BOGOTA.*", "BOGOTA", clave_chr)
  clave_chr <- sub("^GUAJIRA$", "LA_GUAJIRA", clave_chr)
  clave_chr <- sub("^VALLE$", "VALLE_DEL_CAUCA", clave_chr)
  salida_chr <- unname(nombres_departamento_chr[clave_chr])
  ifelse(is.na(salida_chr), x, salida_chr)
}

figura_2a_plot <- humedad_correlacion_dep_df %>%
  mutate(DEP = reorder(nombre_departamento_fn(DEP), r_HR_HA)) %>%
  ggplot(aes(x = r_HR_HA, y = DEP)) +
  geom_vline(xintercept = humedad_resumen_lst$r_global, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = 0, xend = r_HR_HA, yend = DEP), colour = "grey75") +
  geom_point(size = 2.2, colour = "#1f5a85") +
  scale_x_continuous(limits = c(min(0, min(humedad_correlacion_dep_df$r_HR_HA)), 1),
                     labels = coma_decimal_fn) +
  annotate("text", x = humedad_resumen_lst$r_global, y = nrow(humedad_correlacion_dep_df) + 0.8,
           hjust = -0.05, size = 2.8,
           label = paste0("r agrupado = ", formatC(humedad_resumen_lst$r_global, format = "f",
                                                   digits = 2, decimal.mark = ","))) +
  coord_cartesian(clip = "off") +
  labs(x = "r intradepartamental (HR, HA)", y = NULL) +
  tema_articulo

medias_cargas_df <- as.data.frame(cargas_correlacion(medias_acp_prcomp)) %>%
  setNames(c("PC1", "PC2")) %>%
  mutate(variable = c("Temperatura media", "Amplitud t\u00e9rmica", "HR", "HA"))
escala_flecha_num <- 0.9 * max(abs(tipologia_df[, c("PC1", "PC2")])) /
  max(abs(medias_cargas_df[, c("PC1", "PC2")]))
medias_var_pct <- round(100 * medias_acp_prcomp$sdev^2 / sum(medias_acp_prcomp$sdev^2), 1)
pct_coma_chr <- formatC(medias_var_pct, format = "f", digits = 1, decimal.mark = ",")

# Group labels from the data: the group with the lower mean temperature is the Andean one
t_media_grupo_num <- tapply(exposicion_medias_df$T_med, tipologia_df$grupo_ward, mean)
grupo_andino_int <- as.integer(names(which.min(t_media_grupo_num)))
tipologia_figura_df <- tipologia_df %>%
  mutate(DEP = nombre_departamento_fn(DEP),
         grupo_chr = if_else(grupo_ward == grupo_andino_int, "Andino", "Tierras c\u00e1lidas"))

figura_2b_plot <- ggplot(tipologia_figura_df, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_segment(data = medias_cargas_df,
               aes(x = 0, y = 0, xend = PC1 * escala_flecha_num, yend = PC2 * escala_flecha_num),
               arrow = arrow(length = unit(2, "mm")), colour = "grey35") +
  geom_text_repel(data = medias_cargas_df,
                  aes(x = PC1 * escala_flecha_num * 1.1, y = PC2 * escala_flecha_num * 1.1,
                      label = variable), size = 3, colour = "grey20", fontface = "italic",
                  seed = 2, box.padding = 0.4, min.segment.length = Inf) +
  geom_point(aes(shape = grupo_chr), size = 2.4, colour = "#1f5a85") +
  geom_text_repel(aes(label = DEP), size = 2.6, colour = "grey10", seed = 1,
                  min.segment.length = 0.3, box.padding = 0.35, max.overlaps = Inf) +
  scale_x_continuous(expand = expansion(mult = 0.15), labels = coma_decimal_fn) +
  scale_y_continuous(labels = coma_decimal_fn) +
  coord_cartesian(clip = "off") +
  labs(x = paste0("CP1 (", pct_coma_chr[1], " %)"), y = paste0("CP2 (", pct_coma_chr[2], " %)"),
       shape = "Grupo") +
  tema_articulo

figura_scree_plot <- retencion_df %>%
  select(componente, observado = autovalor, paralelo_p95) %>%
  pivot_longer(-componente, names_to = "serie", values_to = "autovalor") %>%
  ggplot(aes(componente, autovalor, group = serie, linetype = serie)) +
  geom_hline(yintercept = 1, colour = "grey70") +
  geom_line() + geom_point() +
  scale_linetype_manual(values = c(observado = "solid", paralelo_p95 = "dashed"),
                        labels = c(observado = "Observado",
                                   paralelo_p95 = "An\u00e1lisis paralelo (percentil 95)")) +
  scale_x_discrete(labels = function(x) sub("PC", "CP", x)) +
  scale_y_continuous(labels = coma_decimal_fn) +
  labs(x = NULL, y = "Autovalor", linetype = NULL,
       title = "Sedimentaci\u00f3n y an\u00e1lisis paralelo de Horn") +
  tema_articulo

ggsave(file.path(resultados_dir, "figura_1.png"), figura_2a_plot, width = 4, height = 4, dpi = 600)
ggsave(file.path(resultados_dir, "figura_2.png"), figura_2b_plot, width = 5, height = 4, dpi = 600)
ggsave(file.path(resultados_dir, "figura_S1_sedimentacion.png"), figura_scree_plot,
       width = 5, height = 3.5, dpi = 600)

# RPSP requests editable figures (.eps) exported from the statistical software
ggsave(file.path(resultados_dir, "figura_1.eps"), figura_2a_plot, width = 4, height = 4,
       device = cairo_ps)
ggsave(file.path(resultados_dir, "figura_2.eps"), figura_2b_plot, width = 5, height = 4,
       device = cairo_ps)


# __INDEX_12__ Export ---------------------------------------------------------------------------
message("[12] Exporting tables")

tablas_final_lst <- list(
  calidad = calidad_resumen_df, completitud = completitud_df,
  humedad_por_departamento = humedad_correlacion_dep_df, retencion = retencion_df,
  sensibilidad = sensibilidad_df, bootstrap = bootstrap_resumen_df,
  estabilidad = bind_rows(estabilidad_temporada_df, estabilidad_mensual_df),
  silueta = silueta_df, tipologia = tipologia_df,
  correlacion_global = as.data.frame(correlacion_global_mat) %>% tibble::rownames_to_column("var"),
  correlacion_intra = as.data.frame(correlacion_intra_mat) %>% tibble::rownames_to_column("var"),
  correlacion_entre = as.data.frame(correlacion_entre_mat) %>% tibble::rownames_to_column("var"),
  cargas = tibble(variable = exposicion_vars_chr,
                  global_pc1 = exposicion_cargas_mat[, 1], global_pc2 = exposicion_cargas_mat[, 2],
                  intra_pc1 = intra_cargas_mat[, 1], intra_pc2 = intra_cargas_mat[, 2],
                  vector_pc1 = exposicion_acp_prcomp$rotation[, 1],
                  vector_pc2 = exposicion_acp_prcomp$rotation[, 2]),
  congruencia_intra = tibble(componente = paste0("PC", seq_len(componentes_k_int)),
                             congruencia = intra_congruencia_vec)
)
iwalk(tablas_final_lst, ~ write.csv(.x, file.path(resultados_dir, paste0(.y, ".csv")),
                                    row.names = FALSE))
writeLines(capture.output(sessionInfo()), file.path(resultados_dir, "sessionInfo.txt"))


# __INDEX_13__ Final summary --------------------------------------------------------------------
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

message("\n==================== SUMMARY ====================")
message("Panel: ", length(departamentos_incluidos_chr), " departments, ",
        nrow(exposicion_diaria_df), " department-days")
message("HR-HA r: pooled ", fmt(humedad_resumen_lst$r_global),
        " | within-department median ", fmt(humedad_resumen_lst$r_intra_mediana),
        " (range ", fmt(humedad_resumen_lst$r_intra_rango[1]), " to ",
        fmt(humedad_resumen_lst$r_intra_rango[2]), ")",
        " | between-department ", fmt(humedad_resumen_lst$r_entre))
message("Eigenvalues: ", paste(fmt(retencion_df$autovalor), collapse = ", "),
        " | Horn p95: ", paste(fmt(retencion_df$paralelo_p95), collapse = ", "))
message("Retained: Kaiser ", sum(retencion_df$kaiser_retener),
        " | parallel analysis ", sum(retencion_df$paralelo_retener))
message("Variance PC1+PC2: pooled ", fmt(retencion_df$varianza_acum_pct[2], 1),
        "% | without HA ", fmt(sum(head(filter(sensibilidad_df,
          analisis == "sin_humedad_absoluta")$varianza_pct, 2)), 1),
        "% | within-department ", fmt(sum(head(filter(sensibilidad_df,
          analisis == "intra_departamento")$varianza_pct, 2)), 1), "%")
message("Bootstrap congruence PC1 / PC2 (median): ",
        fmt(median(bootstrap_df$cong_pc1), 3), " / ", fmt(median(bootstrap_df$cong_pc2), 3))
message("Monthly congruence PC1 range: ",
        paste(fmt(range(estabilidad_mensual_df$cong_pc1), 3), collapse = " to "),
        " | PC2 range: ",
        paste(fmt(range(estabilidad_mensual_df$cong_pc2), 3), collapse = " to "))
message("Typology: k = ", grupos_k_int, " (Ward silhouette ",
        fmt(max(silueta_df$silueta_ward)), ")")
message("Outputs written to: ", resultados_dir)
message("=================================================")
