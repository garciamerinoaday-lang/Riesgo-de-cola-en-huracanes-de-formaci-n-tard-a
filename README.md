# Riesgo de cola en huracanes de formación tardía

Estudio exploratorio con **bootstrap** y **teoría de valores extremos (POT/GPD)** sobre las pérdidas de los huracanes atlánticos que alcanzan la categoría de huracán el 1 de octubre o después, con el huracán **Isaias (2026)** como caso motivador.

<p align="center">
  <img src="figuras/fig3_cola.png" width="720" alt="La cola según cada modelo y los datos de toda la temporada">
</p>

## ¿De qué va?

Isaias fue el primer huracán de la temporada atlántica más tardío en más de un siglo, en un año muy tranquilo, y aun así tocó tierra en Florida con pérdidas relevantes. Eventos así son raros: casi no hay antecedentes, pero sus consecuencias pesan mucho en el reaseguro de catástrofes. La pregunta del proyecto es **qué se puede afirmar con honestidad sobre el riesgo de cola cuando la muestra es mínima**.

El análisis sigue estos pasos:

1. **Población de referencia** a partir de HURDAT2 (NOAA): huracanes del Atlántico con vientos ≥ 64 kt y primer registro como huracán a partir del 1 de octubre (1996–2025: 61 tardíos de 233 huracanes).
2. **Pérdidas por evento** a partir de NOAA Storm Events, emparejando los daños con cada huracán por la mención de su nombre en la narrativa del evento, dentro de una ventana temporal, y normalizando por inflación (IPC) a dólares de 2025.
3. **Bootstrap** de la media, los percentiles altos y la probabilidad de superar umbrales, con la muestra de tardíos con pérdida positiva (n = 24).
4. **EVT (Peaks Over Threshold)**: ajuste de una GPD, VaR y TVaR en forma cerrada, niveles de retorno y bootstrap de los parámetros, primero solo con tardíos.
5. **Modelo híbrido**: la forma de la cola (ξ, β) se estima con todos los huracanes de la temporada y la frecuencia de excedencias solo con los tardíos.

## Resultados principales

| | Solo tardíos | Modelo híbrido |
|---|---|---|
| Índice de cola ξ (estimado) | 0,72 | 0,09 |
| IC 95 % bootstrap de ξ | [−0,5; 2,9] | [−0,5; 0,5] |
| VaR 99,5 % por huracán tardío, mediana bootstrap (MUSD de 2025) | ≈ 39 200 | ≈ 39 200 |
| IC 95 % del VaR 99,5 % (MUSD de 2025) | 15 000 – 1 200 000 | 16 000 – 63 000 |

- Tres huracanes (Sandy, Michael y Wilma) concentran el 87 % de la pérdida normalizada: la media histórica es muy inestable.
- El bootstrap no puede extrapolar más allá del mayor evento observado.
- El POT con solo ~7 excedencias es inestable; compartir la forma de la cola con el resto de la temporada reduce mucho la incertidumbre, aunque no cambia la mediana del VaR.
- La incertidumbre restante procede sobre todo de la frecuencia: solo 4 tardíos superan el umbral del modelo híbrido.

Los detalles, los contrastes de igualdad entre tardíos y resto y las limitaciones están en el [informe](informe/informe_huracanes_revisado.pdf).

## Estructura del repositorio

```
.
├── README.md
├── R/
│   └── huracanes_tardios_bootstrap.R   # análisis completo (R base)
├── python/
│   ├── figs.py                         # figuras del informe (suavizado gaussiano)
│   └── linkedin.py                     # imagen para redes sociales
├── informe/
│   ├── informe_huracanes_revisado.tex
│   └── informe_huracanes_revisado.pdf
├── figuras/                            # PNG generados por figs.py
└── data/                               # datos descargados (no versionados)
```

## Datos

Los datos son abiertos y **no se incluyen** en el repositorio por tamaño:

- **HURDAT2** (NOAA/NHC): <https://www.nhc.noaa.gov/data/#hurdat>. El script lo descarga automáticamente; el nombre del fichero cambia con cada actualización, así que si falla, descárgalo a mano y guárdalo como `hurdat2.txt`.
- **Storm Events Database** (NOAA/NCEI): <https://www.ncei.noaa.gov/access/storm-events-database/>. Se usan los ficheros `StormEvents_details-ftp_v1.0_dAAAA_c*.csv.gz` de 1996 a 2025, en una carpeta `StormEvents/`. Antes de 1996 la base solo recoge unos pocos tipos de evento, por eso el periodo de estudio empieza ese año.
- **IPC (CPI-U)**, media anual: la serie está incluida en el script. Solo los valores de 1996 y 2025 están contrastados externamente; conviene verificar el resto con FRED (CPIAUCNS) si se va a citar.

## Cómo reproducir

### Análisis en R

Requiere R (solo paquetes base). Desde la raíz del proyecto:

```r
# 1. Coloca los CSV.gz de Storm Events en StormEvents/ (o ajusta SE_LOCAL_DIR)
# 2. Ejecuta el script completo
source("R/huracanes_tardios_bootstrap.R")
```

Los parámetros principales (periodo, mes de corte, número de réplicas, umbrales, modo de emparejamiento, cuantil del umbral híbrido) están al principio del script. La semilla es `set.seed(2026)`, con 10 000 réplicas en el bootstrap empírico y 2 000 en los de EVT. La lectura de ~30 ficheros anuales con narrativas es pesada en memoria y puede tardar varios minutos.

### Figuras

Requiere Python ≥ 3.9 con `numpy`, `scipy` y `matplotlib`:

```bash
python python/figs.py
```

Las figuras usan núcleos gaussianos (ancho de banda elegido por validación cruzada de la verosimilitud). Las distribuciones bootstrap de las figuras se regeneran en Python con los mismos datos, por lo que su forma coincide con la del script R aunque la semilla no sea idéntica. Los puntos de toda la temporada de la figura de la cola se reconstruyen con datos exactos de la salida de R. `linkedin.py` escribe en una ruta absoluta de salida: ajústala antes de ejecutarlo.

### Informe

Se compila con XeLaTeX o LuaLaTeX (usa `fontspec`); con pdfLaTeX hay que cambiar el preámbulo como indican los comentarios del `.tex`.

## Galería

| | |
|---|---|
| ![Distribución de las pérdidas](figuras/fig1_distribucion.png) | ![Bootstrap](figuras/fig2_bootstrap.png) |
| ![Estabilidad de ξ](figuras/fig4_estabilidad.png) | ![Comparación del VaR](figuras/fig5_var_comparacion.png) |

## Limitaciones

Es un ejercicio **exploratorio**, no un modelo de tarificación ni una estimación de capital regulatorio. Las pérdidas son daños directos (no asegurados), nominales y corregidas solo por IPC, sin ajuste por exposición; Storm Events infravalora los grandes eventos; el emparejamiento por narrativa es conservador; la muestra es muy pequeña (solo cuatro tardíos superan el umbral del modelo híbrido); y no se trata la dependencia entre eventos ni la falta de estacionariedad. La comparación con Isaias es orientativa: sus estimaciones son pérdidas aseguradas preliminares de distintos modelizadores y Isaias no forma parte de la muestra.

## Referencias

- Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.
- de Vicente Maldonado, J. (2026). *Teoría de los Valores Extremos*. Apuntes de Estadística Actuarial, Universidad Carlos III de Madrid.
- NOAA/NHC, HURDAT2; NOAA/NCEI, Storm Events Database.

## Autor

**Aday García Merino**

## Licencia

Pendiente de definir: añade un fichero `LICENSE` (por ejemplo, MIT para el código).
