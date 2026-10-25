# Aday García Merino,  25/10/2026
# =============================================================================
# Estudio exploratorio del riesgo de huracanes tardíos (bootstrap)
# Pasos: (1) población HURDAT2  (2) variables físicas + ACE
#        (3) pérdidas NOAA Storm Events  (4) bootstrap  (5) resultados
# Solo R base. Probado en lógica contra el formato HURDAT2 real.
# ============================================================================= 
options(timeout = 1200, stringsAsFactors = FALSE)
set.seed(2026)

# ----------------------------- 0. PARÁMETROS --------------------------------
HURDAT_URL   <- "https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2025-092326.txt"
HURDAT_LOCAL <- "hurdat2.txt"          # si lo descargas a mano, ponlo aquí
SE_DIR_URL   <- "https://www.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles/"
SE_LOCAL_DIR <- "StormEvents"          # carpeta con los StormEvents_details-*.csv.gz (la del zip)

MES_CORTE    <- 10     # "tardío": primer registro como huracán desde el 1 de octubre
ANIO_MIN     <- 1996   # Storm Events recoge todos los tipos de evento desde 1996
ANIO_MAX     <- 2025   # último año completo en HURDAT2 actual
MODO_PERDIDA <- "nombre"  # "nombre" (recomendado) | "core" | "core+narrativa"
# "nombre": solo suma eventos cuya narrativa menciona el nombre del ciclón (evita
# atribuir a un huracán los daños de otro con ventana de fechas solapada).
B            <- 10000  # réplicas bootstrap
UMBRALES_USD <- c(1e8, 5e8, 1e9, 5e9)  # umbrales de pérdida para P(L > u)

dir.create(SE_LOCAL_DIR, showWarnings = FALSE)

# ------------------- 1. LECTURA Y PARSEO DE HURDAT2 -------------------------
if (!file.exists(HURDAT_LOCAL)) download.file(HURDAT_URL, HURDAT_LOCAL, mode = "wb")
lin <- readLines(HURDAT_LOCAL, warn = FALSE)

es_cab <- grepl("^[A-Z]{2}[0-9]{6},", lin)       # líneas de cabecera de cada ciclón
cab    <- do.call(rbind, strsplit(lin[es_cab], ","))
cab    <- data.frame(id = trimws(cab[, 1]), nombre = trimws(cab[, 2]))

obs <- read.table(text = lin[!es_cab], sep = ",", strip.white = TRUE, fill = TRUE,
                  header = FALSE, quote = "", comment.char = "",
                  col.names = c("fecha", "hora", "ident", "estado", "lat", "lon",
                                "viento", "presion", paste0("r", 1:13)))
obs$id <- rep(cab$id, times = rle(cumsum(es_cab)[!es_cab])$lengths)
obs$nombre <- cab$nombre[match(obs$id, cab$id)]

obs$fecha <- as.Date(as.character(obs$fecha), "%Y%m%d")
obs$lat <- ifelse(grepl("S", obs$lat), -1, 1) * as.numeric(sub("[NS]", "", obs$lat))
obs$lon <- ifelse(grepl("W", obs$lon), -1, 1) * as.numeric(sub("[EW]", "", obs$lon))
obs$viento[obs$viento < 0]   <- NA
obs$presion[obs$presion < 0] <- NA
obs <- obs[grepl("^AL", obs$id), ]                # cuenca atlántica

# ------------- 2. POBLACIÓN DE HURACANES TARDÍOS + VARIABLES FÍSICAS --------
por_ciclon <- split(obs, obs$id)

resumen <- do.call(rbind, lapply(por_ciclon, function(d) {
  hu <- d[d$estado == "HU" | (!is.na(d$viento) & d$viento >= 64), ]
  if (nrow(hu) == 0) return(NULL)                  # nunca fue huracán
  tc <- d[d$estado %in% c("TS", "HU", "SS"), ]
  # ACE: 1e-4 * sum(v^2) sobre registros sinópticos (00,06,12,18Z) con v >= 34 kt
  sin <- d[d$hora %in% c(0, 600, 1200, 1800) & !is.na(d$viento) &
             d$viento >= 34 & d$estado %in% c("TS", "HU", "SS"), ]
  data.frame(
    id          = d$id[1],
    nombre      = d$nombre[1],
    anio        = as.integer(format(min(hu$fecha), "%Y")),
    fecha_HU    = min(hu$fecha),
    fecha_ini   = min(tc$fecha),
    fecha_fin   = max(d$fecha),
    viento_max  = max(d$viento, na.rm = TRUE),
    pres_min    = if (all(is.na(d$presion))) NA else min(d$presion, na.rm = TRUE),
    ACE         = 1e-4 * sum(sin$viento^2, na.rm = TRUE),
    landfall_L  = any(d$ident == "L")              # registro de tocada de tierra
  )
}))
rownames(resumen) <- NULL

tardios <- subset(resumen,
                  as.integer(format(fecha_HU, "%m")) >= MES_CORTE &
                    anio >= ANIO_MIN & anio <= ANIO_MAX)

cat("\n== Población de referencia: huracanes con primer registro HU desde el 1-oct ==\n")
cat("Total en 1851-2025 :", sum(as.integer(format(resumen$fecha_HU, "%m")) >= MES_CORTE), "\n")
cat("En", ANIO_MIN, "-", ANIO_MAX, ":", nrow(tardios), "\n")
print(tardios[order(tardios$fecha_HU), c("id", "nombre", "fecha_HU", "viento_max",
                                         "pres_min", "ACE", "landfall_L")], row.names = FALSE)

# ----------------- 3. PÉRDIDAS: NOAA STORM EVENTS (por evento) --------------
# Convierte "10.00K", "1.50M", "2.00B" en dólares
parse_dano <- function(x) {
  x <- toupper(trimws(x)); x[x == "" | is.na(x)] <- "0"
  mult <- c(K = 1e3, M = 1e6, B = 1e9)
  u <- substr(x, nchar(x), nchar(x))
  m <- ifelse(u %in% names(mult), mult[u], 1)
  num <- suppressWarnings(as.numeric(ifelse(u %in% names(mult), substr(x, 1, nchar(x) - 1), x)))
  num[is.na(num)] <- 0
  num * m
}

# Localiza/descarga el CSV.gz de detalles de un año
fichero_se <- function(anio) {
  loc <- list.files(SE_LOCAL_DIR, pattern = sprintf("details.*_d%d_.*csv\\.gz$", anio), full.names = TRUE)
  if (length(loc) > 0) return(loc[1])
  idx <- readLines(SE_DIR_URL, warn = FALSE)             # listado del directorio
  m <- regmatches(idx, regexpr(sprintf("StormEvents_details-ftp_v1\\.0_d%d_c[0-9]+\\.csv\\.gz", anio), idx))
  if (length(m) == 0) stop("No encuentro el fichero de ", anio, ". Descárgalo a mano en ", SE_LOCAL_DIR)
  dest <- file.path(SE_LOCAL_DIR, m[1])
  download.file(paste0(SE_DIR_URL, m[1]), dest, mode = "wb", quiet = TRUE)
  dest
}

leer_se <- function(anio) {
  f <- fichero_se(anio)
  hdr <- names(read.csv(gzfile(f), nrows = 1))
  keep <- c("BEGIN_YEARMONTH", "BEGIN_DAY", "END_YEARMONTH", "END_DAY", "STATE",
            "EVENT_TYPE", "DAMAGE_PROPERTY", "DAMAGE_CROPS",
            "EPISODE_NARRATIVE", "EVENT_NARRATIVE")
  cc <- ifelse(hdr %in% keep, "character", "NULL")
  d <- read.csv(gzfile(f), colClasses = cc)
  d$f_ini <- as.Date(sprintf("%s%02d", d$BEGIN_YEARMONTH, as.integer(d$BEGIN_DAY)), "%Y%m%d")
  d$f_fin <- as.Date(sprintf("%s%02d", d$END_YEARMONTH,   as.integer(d$END_DAY)),   "%Y%m%d")
  d$dano_prop  <- parse_dano(d$DAMAGE_PROPERTY)
  d$dano_cult  <- parse_dano(d$DAMAGE_CROPS)
  d$dano_total <- d$dano_prop + d$dano_cult
  # solo estados/territorios (fuera zonas marinas y lagos)
  d[!grepl("ATLANTIC|GULF|LAKE|PACIFIC|HAWAII WATERS|AMERICAN SAMOA|GUAM", d$STATE), ]
}

TIPOS_CORE <- c("Hurricane (Typhoon)", "Hurricane", "Tropical Storm",
                "Tropical Depression", "Storm Surge/Tide")
TIPOS_EXTRA <- c("Heavy Rain", "Flash Flood", "Flood", "Coastal Flood", "Tornado",
                 "High Wind", "Strong Wind", "Thunderstorm Wind", "Lightning",
                 "Rip Current", "Debris Flow")

anios <- sort(unique(tardios$anio))
se <- setNames(lapply(anios, leer_se), anios)

asignar_perdida <- function(s) {
  d <- se[[as.character(s$anio)]]
  # ventana temporal: vida del ciclón (+3 días para inundaciones interiores)
  en_vent <- d$f_ini >= s$fecha_ini & d$f_ini <= s$fecha_fin + 3
  core <- en_vent & d$EVENT_TYPE %in% TIPOS_CORE & d$f_ini <= s$fecha_fin + 1
  L_core <- sum(d$dano_total[core])
  L_nar <- 0; L_nombre <- 0
  if (s$nombre != "UNNAMED") {
    # Solo miramos las narrativas de los pocos eventos candidatos (ventana + tipo),
    # no de las ~60.000 filas del año: es lo que hacía el código eterno.
    cand <- which(en_vent & d$EVENT_TYPE %in% c(TIPOS_CORE, TIPOS_EXTRA))
    if (length(cand) > 0) {
      nom <- toupper(s$nombre)
      txt <- toupper(paste(d$EPISODE_NARRATIVE[cand], d$EVENT_NARRATIVE[cand]))
      pre <- grepl(nom, txt, fixed = TRUE)           # filtro rápido (sin regex)
      menciona <- rep(FALSE, length(cand))
      if (any(pre))                                  # confirmación con límite de palabra
        menciona[pre] <- grepl(paste0("\\b", nom, "\\b(?!\\s+(COUNTY|PARISH|RIVER|LAKE|CITY|TOWNSHIP|SPRINGS|BAY|ISLAND|CREEK))"),
                               txt[pre], perl = TRUE)
      L_nombre <- sum(d$dano_total[cand[menciona]])
      L_nar <- sum(d$dano_total[cand[menciona & d$EVENT_TYPE[cand] %in% TIPOS_EXTRA]])
    }
  }
  data.frame(id = s$id, L_core = L_core, L_core_narr = L_core + L_nar, L_nombre = L_nombre,
             n_ev_core = sum(core), estados = paste(sort(unique(d$STATE[core])), collapse = "; "))
}

perd <- do.call(rbind, lapply(seq_len(nrow(tardios)), function(i) asignar_perdida(tardios[i, ])))
base <- merge(tardios, perd, by = "id")
base$L <- switch(MODO_PERDIDA, nombre = base$L_nombre, core = base$L_core, base$L_core_narr)

cat("\n== Pérdidas asignadas (USD nominales, daños directos NOAA) ==\n")
print(base[order(-base$L), c("nombre", "anio", "viento_max", "ACE", "L_core", "L_core_narr", "L_nombre")],
      row.names = FALSE)

# ------------------ 3b. NORMALIZACIÓN POR INFLACIÓN (IPC) --------------------
ANIO_BASE     <- 2025   # dólares constantes de este año
DESCARGAR_CPI <- FALSE  # TRUE: intenta bajar CPIAUCNS de FRED (si falla usa la tabla)

# CPI-U, media anual, todos los artículos, 1982-84 = 100 (BLS).
# Comprueba la serie en FRED (CPIAUCNS) si la vas a citar en el informe.
CPI <- c("1996" = 156.9,  "1997" = 160.5,  "1998" = 163.0,  "1999" = 166.6,
         "2000" = 172.2,  "2001" = 177.1,  "2002" = 179.9,  "2003" = 184.0,
         "2004" = 188.9,  "2005" = 195.3,  "2006" = 201.6,  "2007" = 207.342,
         "2008" = 215.303, "2009" = 214.537, "2010" = 218.056, "2011" = 224.939,
         "2012" = 229.594, "2013" = 232.957, "2014" = 236.736, "2015" = 237.017,
         "2016" = 240.007, "2017" = 245.120, "2018" = 251.107, "2019" = 255.657,
         "2020" = 258.811, "2021" = 270.970, "2022" = 292.655, "2023" = 304.702,
         "2024" = 313.689, "2025" = 321.943)
if (DESCARGAR_CPI) {
  tryCatch({
    f  <- read.csv("https://fred.stlouisfed.org/graph/fredgraph.csv?id=CPIAUCNS")
    an <- as.integer(substr(f[[1]], 1, 4)); v <- suppressWarnings(as.numeric(f[[2]]))
    ok <- !is.na(v)
    n_mes <- tapply(v[ok], an[ok], length)
    med   <- tapply(v[ok], an[ok], mean)
    CPI   <- med[n_mes == 12]            # solo años completos
    message("CPI descargado de FRED")
  }, error = function(e) message("No se pudo bajar el CPI; uso la tabla incluida."))
}
stopifnot(as.character(ANIO_BASE) %in% names(CPI))

base$factor_ipc <- unname(CPI[as.character(ANIO_BASE)] / CPI[as.character(base$anio)])
base$L_real     <- base$L * base$factor_ipc

cat("\n== Pérdidas normalizadas por IPC (USD de", ANIO_BASE, ") ==\n")
top <- base[base$L > 0, c("nombre", "anio", "factor_ipc", "L", "L_real")]
top$L <- round(top$L / 1e6, 1); top$L_real <- round(top$L_real / 1e6, 1)
names(top)[4:5] <- c("MUSD_nominal", "MUSD_real")
print(top[order(-top$MUSD_real), ], row.names = FALSE)

# ------------------------------ 4. BOOTSTRAP --------------------------------
bootstrap_L <- function(L, B, umbrales) {
  n <- length(L)
  Lstar   <- matrix(L[sample.int(n, n * B, replace = TRUE)], nrow = B)   # B x n
  media_b <- rowMeans(Lstar)
  q95_b   <- apply(Lstar, 1, quantile, probs = 0.95)
  q99_b   <- apply(Lstar, 1, quantile, probs = 0.99)
  exc_b   <- sapply(umbrales, function(u) rowMeans(Lstar > u))
  ic <- function(x) quantile(x, c(0.025, 0.975), names = FALSE)
  
  tabla <- data.frame(
    estadistico = c("Media", "P95", "P99"),
    muestral    = c(mean(L), quantile(L, 0.95), quantile(L, 0.99)),
    boot_media  = c(mean(media_b), mean(q95_b), mean(q99_b)),
    IC95_inf    = c(ic(media_b)[1], ic(q95_b)[1], ic(q99_b)[1]),
    IC95_sup    = c(ic(media_b)[2], ic(q95_b)[2], ic(q99_b)[2])
  )
  tabla[, -1] <- round(tabla[, -1] / 1e6, 1)
  
  tabla_exc <- data.frame(
    umbral_MUSD = umbrales / 1e6,
    p_muestral  = sapply(umbrales, function(u) mean(L > u)),
    p_boot      = colMeans(exc_b),
    IC95_inf    = apply(exc_b, 2, function(x) ic(x)[1]),
    IC95_sup    = apply(exc_b, 2, function(x) ic(x)[2])
  )
  list(n = n, tabla = tabla, exc = round(tabla_exc, 3),
       media_b = media_b, q95_b = q95_b, q99_b = q99_b, ic = ic)
}

L_nom  <- base$L[base$L > 0]          # muestra: huracanes tardíos con impacto (L > 0)
L_real <- base$L_real[base$L > 0]
n <- length(L_nom)
cat("\nTamaño muestral n =", n, " (huracanes tardíos con pérdida > 0)\n")
if (n < 10) warning("n muy pequeño: los percentiles bootstrap serán muy poco fiables.")

set.seed(2026); res_nom  <- bootstrap_L(L_nom,  B, UMBRALES_USD)
set.seed(2026); res_real <- bootstrap_L(L_real, B, UMBRALES_USD)  # mismas remuestras

cat("\n== Bootstrap NOMINAL (MUSD corrientes) ==\n");  print(res_nom$tabla,  row.names = FALSE)
cat("\n== Bootstrap NORMALIZADO (MUSD de", ANIO_BASE, ") ==\n"); print(res_real$tabla, row.names = FALSE)
cat("\n== P(L > umbral) NOMINAL (condicionada a impacto en EE. UU.) ==\n"); print(res_nom$exc,  row.names = FALSE)
cat("\n== P(L > umbral) NORMALIZADO (USD de", ANIO_BASE, ") ==\n");         print(res_real$exc, row.names = FALSE)

# Referencia: Isaias (~3.500 MUSD de consenso, pérdida asegurada; solo orientativo)
cat("\nFrecuencia de pérdidas normalizadas > 3.500 MUSD:", mean(L_real > 3.5e9),
    " (", sum(L_real > 3.5e9), "de", n, ")\n")

# ------------------------------- 5. GRÁFICOS --------------------------------
par(mfrow = c(2, 2))
hist(log10(L_real), breaks = 10, col = "grey85",
     main = paste("Pérdidas observadas (USD", ANIO_BASE, ")"), xlab = "log10(USD)")
for (v in list(list(res_real$media_b, "Media"), list(res_real$q95_b, "P95"),
               list(res_real$q99_b, "P99"))) {
  hist(v[[1]] / 1e6, breaks = 50, col = "steelblue", border = "white",
       main = paste("Bootstrap normalizado:", v[[2]]), xlab = "MUSD")
  abline(v = res_real$ic(v[[1]]) / 1e6, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))

# =============================================================================
# 6. EVT: PEAKS OVER THRESHOLD (GPD) sobre las pérdidas normalizadas
#    Sigue los apuntes de Teoría de Valores Extremos (UC3M, sec. 2.2-2.3 y 4.2)
# =============================================================================
K_EXC         <- 7      # nº de excedencias: u = (k+1)-ésima mayor pérdida
K_GRID        <- 4:12   # rejilla de k para el gráfico de estabilidad
B_EVT         <- 2000   # réplicas bootstrap (incertidumbre de parámetros y riesgos)
NIVELES       <- c(0.95, 0.99, 0.995)
ANIOS         <- ANIO_MAX - ANIO_MIN + 1
UMBRAL_ISAIAS <- 3500   # MUSD (consenso, pérdida asegurada: solo orientativo)

X    <- sort(L_real) / 1e6                      # MUSD de ANIO_BASE, solo L > 0
nX   <- length(X)
Xtot <- c(rep(0, nrow(base) - nX), X)           # todos los tardíos (incluye ceros)
stopifnot(K_EXC %in% K_GRID, nX > max(K_GRID))

# ---- 6.1 GPD por máxima verosimilitud (sec. 2.3) ----------------------------
# Restringimos xi > -0.5 (régimen regular; con xi < -1 la verosimilitud no está acotada)
gpd_nll <- function(par, y) {
  xi <- par[1]; b <- exp(par[2])
  if (xi <= -0.5) return(1e10)
  if (abs(xi) < 1e-8) return(length(y) * log(b) + sum(y) / b)
  z <- 1 + xi * y / b
  if (any(z <= 0)) return(1e10)
  length(y) * log(b) + (1 / xi + 1) * sum(log(z))
}
gpd_fit <- function(y) {
  l0  <- log(mean(y))
  ini <- list(c(0.3, l0), c(0.8, l0), c(-0.2, l0))   # varios arranques
  fits <- lapply(ini, function(p0)
    optim(p0, gpd_nll, y = y, method = "Nelder-Mead",
          control = list(maxit = 5000, reltol = 1e-12)))
  best <- fits[[which.min(sapply(fits, function(f) f$value))]]
  c(xi = best$par[[1]], beta = exp(best$par[[2]]))
}

# ---- 6.2 VaR y TVaR por POT (sec. 4.2) y probabilidad de excedencia ----------
var_pot <- function(p, u, xi, b, k, n) {
  if (k / n <= 1 - p) return(NA_real_)   # cuantil por debajo de u: la GPD no aplica
  if (abs(xi) < 1e-8) u + b * log(k / (n * (1 - p)))
  else u + b / xi * ((n / k * (1 - p))^(-xi) - 1)
}
tvar_pot <- function(p, u, xi, b, k, n) {
  if (xi >= 1) return(Inf)                           # sin media: TVaR infinito
  v <- var_pot(p, u, xi, b, k, n)
  if (abs(xi) < 1e-8) v + b else v / (1 - xi) + (b - xi * u) / (1 - xi)
}
p_exc <- function(x, u, xi, b, k, n) {               # P(L > x), x >= u
  if (abs(xi) < 1e-8) return(k / n * exp(-(x - u) / b))
  z <- 1 + xi * (x - u) / b
  if (z <= 0) 0 else k / n * z^(-1 / xi)
}
nivel_retorno <- function(T, u, xi, b, lambda) {     # perdida superada 1 vez cada T años
  if (abs(xi) < 1e-8) u + b * log(lambda * T) else u + b / xi * ((lambda * T)^xi - 1)
}

# ---- 6.3 Elección del umbral: excesos medios y estabilidad (sec. 2.3.1) ------
diag_k <- do.call(rbind, lapply(K_GRID, function(kk) {
  uu <- X[nX - kk]; yy <- X[(nX - kk + 1):nX] - uu
  f  <- gpd_fit(yy)
  data.frame(k = kk, u = round(uu, 1), xi = round(f[[1]], 2), beta = round(f[[2]], 0),
             VaR995_cond  = round(var_pot(0.995, uu, f[[1]], f[[2]], kk, nX)),
             P_Isaias_cond = round(p_exc(UMBRAL_ISAIAS, uu, f[[1]], f[[2]], kk, nX), 3))
}))
cat("\n== EVT: estabilidad de los parámetros al variar el umbral (MUSD) ==\n")
print(diag_k, row.names = FALSE)

# ---- 6.4 Ajuste final con k = K_EXC -----------------------------------------
u  <- X[nX - K_EXC]
y  <- X[X > u] - u
k  <- length(y)
f  <- gpd_fit(y); xi <- f[[1]]; b <- f[[2]]
cat("\n== Ajuste GPD final ==\n")
cat(sprintf("u = %.1f MUSD | k = %d excedencias | xi = %.3f | beta = %.0f\n", u, k, xi, b))
if (xi >= 1) cat("Aviso: xi >= 1 -> la GPD ajustada no tiene media (TVaR infinito).\n")

nbase <- c("condicionada a impacto (n=24)" = nX, "todos los tardios (n=61)" = length(Xtot))
res_pot <- do.call(rbind, lapply(names(nbase), function(nm) {
  data.frame(base = nm, p = NIVELES,
             VaR  = round(sapply(NIVELES, function(p) var_pot(p, u, xi, b, k, nbase[[nm]]))),
             TVaR = round(sapply(NIVELES, function(p) tvar_pot(p, u, xi, b, k, nbase[[nm]]))))
}))
cat("\n== VaR y TVaR por POT (MUSD de", ANIO_BASE, ") ==\n"); print(res_pot, row.names = FALSE)
cat("\nComparación con el bootstrap empírico (P99 muestral, MUSD):",
    round(quantile(X, 0.99)), " | máximo muestral:", round(max(X)), "\n")

cat("\n== P(L >", UMBRAL_ISAIAS, "MUSD) ==\n")
cat(sprintf("GPD condicionada: %.3f (empírica %.3f) | GPD sobre todos los tardíos: %.3f (empírica %.3f)\n",
            p_exc(UMBRAL_ISAIAS, u, xi, b, k, nX), mean(X > UMBRAL_ISAIAS),
            p_exc(UMBRAL_ISAIAS, u, xi, b, k, length(Xtot)), mean(Xtot > UMBRAL_ISAIAS)))

lambda <- k / ANIOS                                  # excedencias de u por año
rl <- data.frame(T_anios = c(10, 50, 100, 200),
                 nivel_MUSD = round(sapply(c(10, 50, 100, 200),
                                           function(T) nivel_retorno(T, u, xi, b, lambda))))
cat("\n== Niveles de retorno (pérdida superada en media 1 vez cada T años) ==\n")
cat(sprintf("(tasa de excedencias de u: %.3f por año)\n", lambda)); print(rl, row.names = FALSE)

# ---- 6.5 Validación por bootstrap (sec. 2.3.1): remuestreo de las X_i --------
boot_pot <- function(x, u, B, niveles) {
  out <- matrix(NA_real_, B, 2 + 2 * length(niveles))
  colnames(out) <- c("xi", "beta", paste0("VaR", niveles), paste0("TVaR", niveles))
  for (i in seq_len(B)) {
    xb <- sample(x, length(x), replace = TRUE)
    yb <- xb[xb > u] - u
    if (length(yb) < 5) next
    fb <- tryCatch(gpd_fit(yb), error = function(e) NULL)
    if (is.null(fb)) next
    kb <- length(yb)
    out[i, ] <- c(fb[[1]], fb[[2]],
                  sapply(niveles, function(p) var_pot(p, u, fb[[1]], fb[[2]], kb, length(x))),
                  sapply(niveles, function(p) tvar_pot(p, u, fb[[1]], fb[[2]], kb, length(x))))
  }
  out[complete.cases(out), , drop = FALSE]
}
resumen_boot <- function(bt) {                       # type = 1: sin aritmética con Inf
  q <- function(p) apply(bt, 2, function(x) quantile(x, probs = p, type = 1, na.rm = TRUE))
  d <- data.frame(medida = colnames(bt), mediana = q(0.5),
                  IC95_inf = q(0.025), IC95_sup = q(0.975), row.names = NULL)
  es_xi <- d$medida == "xi"                          # xi con 2 decimales; el resto, enteros
  d[es_xi, -1]  <- round(d[es_xi, -1], 2)
  d[!es_xi, -1] <- round(d[!es_xi, -1])
  d
}
set.seed(2026); bt_c <- boot_pot(X,    u, B_EVT, NIVELES)
set.seed(2026); bt_t <- boot_pot(Xtot, u, B_EVT, NIVELES)
cat("\n== Bootstrap POT, base condicionada (n=24): réplicas válidas", nrow(bt_c), "de", B_EVT,
    "| P(xi >= 1) =", round(mean(bt_c[, "xi"] >= 1), 2), "==\n")
print(resumen_boot(bt_c), row.names = FALSE)
cat("\n== Bootstrap POT, todos los tardíos (n=61): réplicas válidas", nrow(bt_t), "de", B_EVT,
    "| P(xi >= 1) =", round(mean(bt_t[, "xi"] >= 1), 2), "==\n")
print(resumen_boot(bt_t), row.names = FALSE)

# ---- 6.6 Gráficos de diagnóstico --------------------------------------------
par(mfrow = c(2, 2))
us <- X[1:(nX - 3)]                                   # al menos 3 excedencias
me <- sapply(us, function(v) mean(X[X > v] - v))
plot(us, me, type = "b", pch = 19, xlab = "u (MUSD)", ylab = "e(u)",
     main = "Excesos medios"); abline(v = u, lty = 2, col = "red")
plot(diag_k$u, diag_k$xi, type = "b", pch = 19, xlab = "u (MUSD)", ylab = "xi(u)",
     main = "Estabilidad de xi"); abline(h = 0, lty = 3); abline(v = u, lty = 2, col = "red")
plot(diag_k$u, diag_k$beta, type = "b", pch = 19, xlab = "u (MUSD)", ylab = "beta(u)",
     main = "Estabilidad de beta"); abline(v = u, lty = 2, col = "red")
xs <- sort(X, decreasing = TRUE); emp <- seq_len(nX) / (nX + 1)
xg <- exp(seq(log(u), log(max(X) * 3), length.out = 100))
plot(xs, emp, log = "xy", pch = 19, xlab = "Pérdida (MUSD)", ylab = "P(L > x | impacto)",
     main = "Cola: datos vs GPD", xlim = range(c(xs, xg)))
lines(xg, sapply(xg, function(x) p_exc(x, u, xi, b, k, nX)), col = "red", lwd = 2)
abline(v = u, lty = 2)
par(mfrow = c(1, 1))

# =============================================================================
# 7. EVT CON TODOS LOS HURACANES DE LA TEMPORADA
#    Idea: xi y beta (forma y escala de la cola) se estiman con TODOS los huracanes
#    1996-2025 de cualquier mes -> muchas más excedencias; la FRECUENCIA de excedencias
#    (zeta) se toma solo de los huracanes TARDÍOS.  Es un modelo híbrido:
#       VaR_tardío(p) = u + beta/xi * [ ( N_tard/k_tard * (1-p) )^(-xi) - 1 ]
#    Supuesto clave: la forma de la cola es la misma en toda la temporada (se contrasta abajo).
# =============================================================================
Q_GRID <- seq(0.50, 0.92, by = 0.02)   # cuantiles candidatos para el umbral (sobre L > 0)
Q_U    <- 0.80                         # cuantil elegido (ajústalo tras ver el diagnóstico)
B_POOL <- 2000

# ---- 7.1 Población: todos los huracanes del Atlántico ANIO_MIN-ANIO_MAX --------
todos <- subset(resumen, anio >= ANIO_MIN & anio <= ANIO_MAX)
todos$tardio <- as.integer(format(todos$fecha_HU, "%m")) >= MES_CORTE
cat("\n== Sección 7: huracanes", ANIO_MIN, "-", ANIO_MAX, "(cualquier mes):", nrow(todos),
    "| tardíos:", sum(todos$tardio), "==\n")

# ---- 7.2 Cargar los años de Storm Events que aún faltan -----------------------
falta <- setdiff(unique(todos$anio), as.integer(names(se)))
for (a in falta) se[[as.character(a)]] <- leer_se(a)

# ---- 7.3 Asignar pérdidas (mismo método) y normalizar por IPC ------------------
perd_todos <- do.call(rbind, lapply(seq_len(nrow(todos)),
                                    function(i) asignar_perdida(todos[i, ])))
todos <- merge(todos, perd_todos[, c("id", "L_nombre")], by = "id")
todos$L      <- todos$L_nombre
todos$L_real <- todos$L * unname(CPI[as.character(ANIO_BASE)] / CPI[as.character(todos$anio)])

cat("Control (tardíos): suma de pérdidas normalizadas sección 7 =",
    round(sum(todos$L_real[todos$tardio]) / 1e6), "MUSD | sección 4 =",
    round(sum(L_real) / 1e6), "MUSD (deben coincidir si MODO_PERDIDA = 'nombre')\n")
cat("\n-- Mayores pérdidas de la temporada completa (revisa que sean huracanes reales) --\n")
topP <- todos[order(-todos$L_real), c("nombre", "anio", "tardio", "viento_max", "L_real")][1:15, ]
topP$L_real <- round(topP$L_real / 1e6); names(topP)[5] <- "MUSD_real"
print(topP, row.names = FALSE)

# ---- 7.4 Muestra agrupada ------------------------------------------------------
Pm <- todos[todos$L > 0, ]                       # huracanes (cualquier mes) con impacto
XP <- sort(Pm$L_real) / 1e6
nP <- length(XP); N_all <- nrow(todos); N_tard <- sum(todos$tardio)
cat(sprintf("\nCon impacto (L>0): %d de %d huracanes | de ellos tardíos: %d\n",
            nP, N_all, sum(Pm$tardio)))
stopifnot(nP > 20)

# ---- 7.5 Elección del umbral: estabilidad por cuantiles + excesos medios -------
diag_q <- do.call(rbind, lapply(Q_GRID, function(q) {
  uu <- unname(quantile(XP, q, type = 1)); yy <- XP[XP > uu] - uu
  if (length(yy) < 8) return(NULL)
  ff <- gpd_fit(yy)
  kt <- sum(Pm$tardio & Pm$L_real / 1e6 > uu)
  data.frame(q = q, u = round(uu, 1), k = length(yy), k_tard = kt,
             xi = round(ff[[1]], 2), beta = round(ff[[2]]),
             VaR995_tard = round(var_pot(0.995, uu, ff[[1]], ff[[2]], kt, N_tard)))
}))
cat("\n== Estabilidad de xi, beta y VaR99.5 híbrido (tardíos) al variar el umbral ==\n"); print(diag_q, row.names = FALSE)

# ---- 7.6 Ajuste final con la muestra agrupada ---------------------------------
uP <- unname(quantile(XP, Q_U, type = 1)); yP <- XP[XP > uP] - uP; kP <- length(yP)
fP <- gpd_fit(yP); xiP <- fP[[1]]; bP <- fP[[2]]
cat(sprintf("\n== GPD agrupada: u = %.1f MUSD | k = %d | xi = %.3f | beta = %.0f ==\n",
            uP, kP, xiP, bP))
cat(sprintf("(solo tardíos, sección 6: xi = %.3f con k = %d)\n", xi, k))

# ---- 7.7 ¿La cola de los tardíos es distinta? (contrastes, potencia baja) -------
k_tard <- sum(Pm$tardio & Pm$L_real / 1e6 > uP)
tab2 <- matrix(c(k_tard, N_tard - k_tard, kP - k_tard, (N_all - N_tard) - (kP - k_tard)), 2)
cat(sprintf("\nFrecuencia de excedencias de u: tardíos %d/%d = %.3f | resto %d/%d = %.3f | Fisher p = %.3f\n",
            k_tard, N_tard, k_tard / N_tard, kP - k_tard, N_all - N_tard,
            (kP - k_tard) / (N_all - N_tard), fisher.test(tab2)$p.value))
yT <- (Pm$L_real / 1e6)[Pm$tardio & Pm$L_real / 1e6 > uP] - uP
yR <- (Pm$L_real / 1e6)[!Pm$tardio & Pm$L_real / 1e6 > uP] - uP
if (length(yT) >= 3 && length(yR) >= 3)
  cat(sprintf("Excesos tardíos vs resto (KS, ties ignorados): p = %.3f  (con k_tard=%d la potencia es baja)\n",
              suppressWarnings(ks.test(yT, yR)$p.value), length(yT)))

# ---- 7.8 Medidas de riesgo híbridas para los tardíos --------------------------
stopifnot(k_tard >= 1)
res_hib <- data.frame(
  p = NIVELES,
  VaR_tardio  = round(sapply(NIVELES, function(p) var_pot(p, uP, xiP, bP, k_tard, N_tard))),
  TVaR_tardio = round(sapply(NIVELES, function(p) tvar_pot(p, uP, xiP, bP, k_tard, N_tard))),
  VaR_toda_temporada = round(sapply(NIVELES, function(p) var_pot(p, uP, xiP, bP, kP, N_all))))
cat("\n== VaR/TVaR por huracán tardío (híbrido: xi,beta agrupados; zeta de tardíos), MUSD ==\n")
print(res_hib, row.names = FALSE)

lam_t <- k_tard / ANIOS
cat(sprintf("\nNiveles de retorno (tasa de excedencias de u, tardíos: %.3f/año):\n", lam_t))
print(data.frame(T_anios = c(10, 50, 100, 200),
                 MUSD = round(sapply(c(10, 50, 100, 200),
                                     function(T) nivel_retorno(T, uP, xiP, bP, lam_t)))), row.names = FALSE)
if (UMBRAL_ISAIAS >= uP)
  cat(sprintf("P(L > %d MUSD | tardío con impacto) = %.3f (empírica %.3f)\n", UMBRAL_ISAIAS,
              p_exc(UMBRAL_ISAIAS, uP, xiP, bP, k_tard, sum(Pm$tardio)), mean(X > UMBRAL_ISAIAS)))

# ---- 7.9 Bootstrap por huracán (la unidad es el huracán, no la pérdida) --------
boot_pool <- function(dat, uP, B, niveles) {
  out <- matrix(NA_real_, B, 2 + length(niveles) + 2)
  colnames(out) <- c("xi", "beta", paste0("VaRtard", niveles), "TVaRtard0.995", "k_tard")
  for (i in seq_len(B)) {
    db <- dat[sample.int(nrow(dat), replace = TRUE), ]
    xb <- db$L_real / 1e6; yb <- xb[xb > uP] - uP
    kt <- sum(xb > uP & db$tardio); Nt <- sum(db$tardio)
    if (length(yb) < 8 || kt < 1) next
    fb <- tryCatch(gpd_fit(yb), error = function(e) NULL)
    if (is.null(fb)) next
    out[i, ] <- c(fb[[1]], fb[[2]],
                  sapply(niveles, function(p) var_pot(p, uP, fb[[1]], fb[[2]], kt, Nt)),
                  tvar_pot(0.995, uP, fb[[1]], fb[[2]], kt, Nt), kt)
  }
  out[!is.na(out[, "xi"]), , drop = FALSE]     # NA en VaR = fórmula no aplicable (zeta <= 1-p)
}
set.seed(2026); bt_p <- boot_pool(todos, uP, B_POOL, NIVELES)
cat("\n== Bootstrap híbrido: réplicas válidas", nrow(bt_p), "de", B_POOL,
    "| P(xi >= 1) =", round(mean(bt_p[, "xi"] >= 1), 2), "==\n")
print(resumen_boot(bt_p), row.names = FALSE)
cat("Réplicas donde la fórmula no es aplicable (zeta <= 1-p), por medida:\n")
print(colSums(is.na(bt_p[, grep("VaR", colnames(bt_p))])))

q_xi <- function(x) round(quantile(x, c(0.025, 0.5, 0.975), type = 1), 2)
cat("\nIC 95% de xi (2.5%, mediana, 97.5%):\n")
cat("  solo tardíos (sección 6):", q_xi(bt_c[, "xi"]), "\n")
cat("  muestra agrupada        :", q_xi(bt_p[, "xi"]), "\n")

# ---- 7.10 Gráficos -------------------------------------------------------------
par(mfrow = c(2, 2))
usP <- XP[1:(nP - 5)]
plot(usP, sapply(usP, function(v) mean(XP[XP > v] - v)), type = "b", pch = 19,
     xlab = "u (MUSD)", ylab = "e(u)", main = "Excesos medios (agrupado)")
abline(v = uP, lty = 2, col = "red")
plot(diag_q$u, diag_q$xi, type = "b", pch = 19, xlab = "u (MUSD)", ylab = "xi(u)",
     main = "Estabilidad de xi (agrupado)"); abline(h = 0, lty = 3); abline(v = uP, lty = 2, col = "red")
ord <- order(Pm$L_real, decreasing = TRUE)
xsP <- Pm$L_real[ord] / 1e6; colP <- ifelse(Pm$tardio[ord], "red", "grey40")
xgP <- exp(seq(log(uP), log(max(xsP) * 3), length.out = 100))
plot(xsP, seq_along(xsP) / (nP + 1), log = "xy", pch = 19, col = colP,
     xlab = "Pérdida (MUSD)", ylab = "P(L > x | impacto)", xlim = range(c(xsP, xgP)),
     main = "Cola agrupada (rojo = tardíos)")
lines(xgP, sapply(xgP, function(x) p_exc(x, uP, xiP, bP, kP, nP)), col = "blue", lwd = 2)
abline(v = uP, lty = 2)
hist(bt_c[, "xi"], breaks = 40, col = rgb(1, 0, 0, 0.4), border = "white", xlim = c(-0.5, 3),
     main = "xi: tardíos (rojo) vs agrupado (azul)", xlab = "xi")
hist(bt_p[, "xi"], breaks = 40, col = rgb(0, 0, 1, 0.4), border = "white", add = TRUE)
par(mfrow = c(1, 1))

