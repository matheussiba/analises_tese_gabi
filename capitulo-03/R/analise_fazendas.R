# Inventário MaxiC, uma fazenda por vez (BH ou ER), com o mesmo código.
#
# Entradas em 0-raw/:
#   parcelas_<COD>.gpkg                                 coordenadas e tratamento (fonte única)
#   planilha_coleta_quali-quantitativa_fazenda_<COD>.xlsx  fustes e ficha da parcela
#
# Regras deste script:
#  - Coordenadas e tratamento vêm do GeoPackage. As colunas P1–P4 e TRATAMENTO
#    das planilhas são ignoradas.
#  - Cada parcela tem 30 x 30 m = 900 m² e 150 mudas plantadas.
#  - Biomassa: Ferez et al. (2015), por indivíduo, com a área seccional a 30 cm.

user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(user_lib) && dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(BIOMASS)
  library(vegan)
})

sf::sf_use_s2(FALSE)

# ---- Constantes --------------------------------------------------------------

AREA_PARCELA_M2 <- 900
MUDAS_PLANTADAS <- 150L
# ln(AGBW) = b0 + b_sa ln(SA) + b_ht ln(Ht) + b_wd ln(rho)
# SA em m², Ht em m, rho em g/cm³ (a forma do artigo; com SA em cm² e rho em
# kg/m³ a conta dá milhões de quilos por árvore).
FEREZ <- c(b0 = 6.039, b_sa = 0.945, b_ht = 0.961, b_wd = 1.022)
FRACAO_CARBONO <- 0.47
LIMIAR_DAP_CM <- 5
LIMIAR_MUDA_CM <- 1
# D30 acima disso é erro de digitação: nenhum fuste das duas fazendas tem DAP
# acima de 29 cm. O fuste fica sem área seccional, o indivíduo sem biomassa,
# e o caso vai para d30_suspeito.csv. A planilha não é alterada.
D30_MAX_PLAUSIVEL_CM <- 50
WD_MIN <- 0.08
WD_MAX <- 1.39
N_SIM_WD <- 1000L
SEMENTE <- 260319L

NOMES_FAZENDA <- c(BH = "Belo Horizonte", ER = "Entre Rios")
PAL_TRAT <- c(
  "T1" = "#0072B2", "T2" = "#E69F00", "T3" = "#009E73", "T4" = "#CC79A7",
  "T5" = "#D55E00", "T6" = "#56B4E9", "T7" = "#999933"
)

# ---- Utilidades --------------------------------------------------------------

ascii_lower <- function(x) {
  tolower(stringi::stri_trans_general(trimws(as.character(x)), "Latin-ASCII"))
}

numero_br <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "N/A", "na", "Na", "null")] <- NA_character_
  y <- gsub("[[:space:]]", "", y)
  y <- gsub(",", ".", y)
  suppressWarnings(as.numeric(y))
}

data_excel <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "N/A")] <- NA_character_
  num <- suppressWarnings(as.numeric(y))
  out <- as.Date(rep(NA_real_, length(y)), origin = "1970-01-01")
  serial <- is.finite(num) & num > 20000 & num < 80000
  out[serial] <- as.Date(num[serial], origin = "1899-12-30")
  resto <- !serial & !is.na(y)
  if (any(resto)) out[resto] <- suppressWarnings(as.Date(y[resto]))
  out
}

percentual <- function(x) {
  v <- numero_br(x)
  ifelse(!is.na(v) & v > 1, v / 100, v)
}

sim_nao <- function(x) {
  y <- ascii_lower(x)
  dplyr::case_when(
    y %in% c("sim", "s", "yes") ~ "sim",
    y %in% c("nao", "n", "no") ~ "nao",
    TRUE ~ NA_character_
  )
}

classe_fito <- function(x) {
  y <- ascii_lower(x)
  dplyr::case_when(
    str_detect(y, "avanc") ~ "avancado",
    str_detect(y, "inicial") ~ "inicial",
    str_detect(y, "sadio") ~ "sadio",
    is.na(y) | y == "na" ~ NA_character_,
    TRUE ~ "outro"
  )
}

nome_exibicao <- function(x) {
  y <- str_squish(as.character(x))
  dplyr::case_when(
    ascii_lower(y) %in% c("ni", "n/i", "nao identificada", "nao identificado") ~ "não identificada",
    ascii_lower(y) %in% c("outra_especie", "outra especie") ~ "outra espécie",
    TRUE ~ y
  )
}

id_parcela <- function(x) {
  y <- str_replace(ascii_lower(x), "^p", "")
  as.integer(round(numero_br(y)))
}

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits, decimal.mark = ",", big.mark = ""))
}

fmt_pct <- function(x, digits = 0) {
  ifelse(is.na(x), "—", paste0(fmt_num(100 * x, digits), "%"))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "—", ifelse(p < 0.001, "< 0,001", formatC(p, format = "f", digits = 3, decimal.mark = ",")))
}

ordem_tratamento <- function(x) {
  n <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  n[is.na(n)] <- 10^6
  order(n, x)
}

niveis_tratamento <- function(x) {
  u <- unique(as.character(x))
  u[ordem_tratamento(u)]
}

ajustar_holm <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "holm")
  out
}

shannon_nat <- function(n) {
  n <- n[n > 0]
  if (length(n) == 0) return(NA_real_)
  p <- n / sum(n)
  -sum(p * log(p))
}

predominante <- function(especie, peso) {
  ok <- !is.na(especie) & is.finite(peso) & peso > 0
  if (!any(ok)) return(NA_character_)
  tab <- tapply(peso[ok], especie[ok], sum)
  paste(sort(names(tab)[tab == max(tab)]), collapse = "; ")
}

epsg_utm_sul <- function(lon) {
  zona <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  31960 + min(max(as.integer(zona), 18L), 25L)
}

# ---- Leitura -----------------------------------------------------------------

pasta_entrada <- function(raiz) {
  candidatos <- c(file.path(raiz, "0-raw"), file.path(raiz, "dados"))
  hit <- candidatos[dir.exists(candidatos)]
  if (length(hit) == 0) stop("Não há pasta 0-raw nem dados em ", raiz)
  hit[1]
}

achar_coluna <- function(nomes, alvo, obrigatoria = TRUE) {
  hit <- nomes[ascii_lower(nomes) == ascii_lower(alvo)]
  if (length(hit) == 0) {
    if (obrigatoria) stop("Coluna ausente: ", alvo)
    return(NULL)
  }
  hit[1]
}

# Polígonos e tratamentos: só do GeoPackage.
ler_poligonos <- function(gpkg, cod) {
  pol <- sf::st_read(gpkg, layer = "parcelas", quiet = TRUE) %>% sf::st_transform(4326)
  if (!all(c("parcela", "tratamento") %in% names(pol))) {
    stop("O GeoPackage precisa das colunas parcela e tratamento: ", gpkg)
  }
  pol <- pol %>%
    transmute(
      codigo_fazenda = cod,
      parcela = as.integer(parcela),
      tratamento = str_squish(as.character(tratamento))
    )
  if (any(duplicated(pol$parcela))) stop("Parcela repetida no GeoPackage: ", gpkg)
  epsg <- epsg_utm_sul(sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(pol))))[, 1])
  pol$area_desenhada_m2 <- as.numeric(sf::st_area(sf::st_transform(pol, epsg)))
  cen <- sf::st_coordinates(sf::st_centroid(sf::st_transform(pol, epsg)))
  cen_ll <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(data.frame(x = cen[, 1], y = cen[, 2]), coords = c("x", "y"), crs = epsg), 4326
  ))
  pol$lon <- cen_ll[, 1]
  pol$lat <- cen_ll[, 2]
  pol
}

ler_qualitativa <- function(xlsx, cod) {
  abas <- excel_sheets(xlsx)
  aba <- abas[str_detect(tolower(abas), "qualit")][1]
  q <- read_excel(xlsx, sheet = aba, col_types = "text")
  nm <- names(q)
  col <- function(a) q[[achar_coluna(nm, a)]]
  tibble(
    parcela = id_parcela(col("ID_PARCELA")),
    fazenda = str_squish(as.character(col("FAZENDA"))),
    talhao = as.integer(numero_br(col("TALHAO"))),
    data_coleta = data_excel(col("DATA")),
    solo = str_squish(as.character(col("SOLO"))),
    declividade = str_squish(as.character(col("DECLIVIDADE"))),
    formigas = sim_nao(col("PRES_FORMIGAS")),
    cupim = sim_nao(col("PRES_CUPIM")),
    n_cupinzeiros = numero_br(col("NUMERO_CUPINZEIROS")),
    perc_folha_larga = percentual(col("PERC_FOLHA_LARGA")),
    perc_folha_estreita = percentual(col("PERC_FOLHA_ESTREITA")),
    obs = as.character(col("OBS"))
  ) %>% filter(!is.na(parcela))
}

ler_fustes <- function(xlsx, cod) {
  abas <- excel_sheets(xlsx)
  aba <- abas[str_detect(tolower(abas), "maxic|quant")][1]
  x <- read_excel(xlsx, sheet = aba, col_types = "text")
  nm <- names(x)
  col <- function(a) x[[achar_coluna(nm, a)]]
  regen <- achar_coluna(nm, "REGENERACAO", obrigatoria = FALSE)
  tibble(
    codigo_fazenda = cod,
    parcela = id_parcela(col("ID_PARCELA")),
    linha = as.integer(numero_br(col("N_LINHA"))),
    id_ind = {
      y <- str_squish(as.character(col("ID_IND")))
      y[y %in% c("", "NA", "N/A", "na")] <- NA_character_
      y
    },
    especie = str_squish(as.character(col("ESPECIE"))),
    d30_cm = numero_br(col("D30_CM")),
    dap_cm = numero_br(col("DAP_CM")),
    altura_m = numero_br(col("ALTURA_TOTAL_M")),
    fito = classe_fito(col("ESTADO_FITOSSANITARIO")),
    regeneracao = if (is.null(regen)) NA_character_ else sim_nao(x[[regen]])
  ) %>%
    filter(!is.na(parcela), !is.na(id_ind), !is.na(especie), especie != "")
}

ler_fazenda <- function(raiz, cod) {
  dados <- pasta_entrada(raiz)
  xlsx <- file.path(dados, sprintf("planilha_coleta_quali-quantitativa_fazenda_%s.xlsx", cod))
  gpkg <- file.path(dados, sprintf("parcelas_%s.gpkg", cod))
  if (!file.exists(xlsx)) stop("Falta a planilha: ", xlsx)
  if (!file.exists(gpkg)) stop("Falta o GeoPackage: ", gpkg)
  pol <- ler_poligonos(gpkg, cod)
  quali <- ler_qualitativa(xlsx, cod)
  fustes <- ler_fustes(xlsx, cod)
  so_planilha <- setdiff(quali$parcela, pol$parcela)
  if (length(so_planilha) > 0) {
    stop("Parcelas na planilha sem polígono no GeoPackage: ", paste(so_planilha, collapse = ", "))
  }
  sem_poligono <- setdiff(unique(fustes$parcela), pol$parcela)
  if (length(sem_poligono) > 0) {
    stop("Fustes em parcelas sem polígono no GeoPackage: ", paste(sem_poligono, collapse = ", "))
  }
  parcelas <- pol %>%
    left_join(quali, by = "parcela") %>%
    mutate(fazenda = ifelse(is.na(fazenda), NOMES_FAZENDA[[cod]], fazenda))
  list(
    cod = cod, parcelas = parcelas, fustes = fustes,
    arquivos = c(xlsx = basename(xlsx), gpkg = basename(gpkg))
  )
}

# ---- Indivíduos e densidade da madeira --------------------------------------

agregar_individuos <- function(fustes) {
  gravidade <- c(sadio = 1L, inicial = 2L, avancado = 3L, outro = 2L)
  fustes %>%
    group_by(codigo_fazenda, parcela, id_ind) %>%
    summarise(
      n_especies_id = n_distinct(especie),
      especie = especie[which.max(tabulate(match(especie, unique(especie))))],
      linha = if (all(is.na(linha))) NA_integer_ else min(linha, na.rm = TRUE),
      n_fustes = n(),
      # Área seccional a 30 cm do solo, somada nos fustes do indivíduo (m²).
      # Um fuste sem D30 válido não entra; se nenhum fuste serve, fica NA.
      sa_m2 = if (any(is.finite(d30_cm) & d30_cm > 0 & d30_cm <= D30_MAX_PLAUSIVEL_CM)) {
        v <- is.finite(d30_cm) & d30_cm > 0 & d30_cm <= D30_MAX_PLAUSIVEL_CM
        sum(pi * (d30_cm[v] / 200)^2)
      } else {
        NA_real_
      },
      n_fustes_d30 = sum(is.finite(d30_cm) & d30_cm > 0 & d30_cm <= D30_MAX_PLAUSIVEL_CM),
      dap_max_cm = if (any(is.finite(dap_cm))) max(dap_cm, na.rm = TRUE) else NA_real_,
      altura_m = if (any(is.finite(altura_m))) stats::median(altura_m[is.finite(altura_m)]) else NA_real_,
      regeneracao = if (all(is.na(regeneracao))) NA else any(regeneracao == "sim", na.rm = TRUE),
      fito = {
        g <- gravidade[fito]
        g <- g[!is.na(g)]
        if (length(g) == 0) NA_character_ else names(gravidade)[max(g)]
      },
      .groups = "drop"
    ) %>%
    mutate(
      altura_m = ifelse(is.finite(altura_m) & altura_m > 0, altura_m, NA_real_),
      # Protocolo: muda é o indivíduo com DAP inferior a 1 cm.
      muda = !(is.finite(dap_max_cm) & dap_max_cm >= LIMIAR_MUDA_CM),
      # A ficha de BH não tem coluna de regeneração: tudo conta como plantio.
      origem = ifelse(regeneracao %in% TRUE, "regeneracao", "plantio")
    )
}

densidade_madeira <- function(individuos, cod) {
  tab <- individuos %>% distinct(especie)
  partes <- str_split(str_squish(tab$especie), "\\s+")
  genus <- vapply(partes, function(p) p[1], character(1))
  species <- vapply(partes, function(p) {
    if (length(p) < 2) NA_character_ else paste(p[-1], collapse = " ")
  }, character(1))
  species[ascii_lower(species) %in% c("sp.", "sp", "spp.", "spp")] <- NA_character_
  sem_nome <- ascii_lower(genus) %in% c("ni", "n/i", "outra_especie", "outra")
  genus[sem_nome] <- NA_character_
  species[sem_nome] <- NA_character_
  # A média de fallback é a do próprio sítio.
  m <- getWoodDensity(
    genus = genus, species = species,
    stand = rep(cod, nrow(tab)), region = "SouthAmericaTrop", verbose = FALSE
  )
  tab %>%
    mutate(
      wd = m$meanWD, wd_sd = m$sdWD,
      wd_nivel = as.character(m$levelWD), wd_n = m$nInd
    )
}

# Ferez et al. (2015), por indivíduo. SA em m², Ht em m, rho em g/cm³.
agb_ferez <- function(sa_m2, ht_m, wd) {
  ok <- is.finite(sa_m2) & sa_m2 > 0 & is.finite(ht_m) & ht_m > 0 & is.finite(wd) & wd > 0
  out <- rep(NA_real_, length(sa_m2))
  out[ok] <- exp(
    FEREZ[["b0"]] + FEREZ[["b_sa"]] * log(sa_m2[ok]) +
      FEREZ[["b_ht"]] * log(ht_m[ok]) + FEREZ[["b_wd"]] * log(wd[ok])
  )
  out
}

calcular_biomassa <- function(individuos, wd) {
  individuos %>%
    left_join(wd, by = "especie") %>%
    mutate(
      agb_kg = agb_ferez(sa_m2, altura_m, wd),
      carbono_kg = agb_kg * FRACAO_CARBONO
    )
}

rtnorm_vec <- function(media, desvio, inferior = WD_MIN, superior = WD_MAX) {
  u <- stats::runif(length(media), stats::pnorm(inferior, media, desvio), stats::pnorm(superior, media, desvio))
  stats::qnorm(u, media, desvio)
}

# Só a incerteza da densidade da madeira. Não inclui o erro da equação nem o de medição.
propagar_wd <- function(ind, n_sim = N_SIM_WD, seed = SEMENTE) {
  alvo <- ind %>% filter(is.finite(agb_kg), is.finite(wd), wd > 0)
  vazio <- tibble(parcela = integer(), agb_kg_wd_inf = numeric(), agb_kg_wd_sup = numeric(), agb_kg_wd_sd = numeric())
  if (nrow(alvo) == 0) return(vazio)
  set.seed(seed)
  desvio <- alvo$wd_sd
  desvio[!is.finite(desvio) | desvio < 0] <- 0
  n <- nrow(alvo)
  mult <- matrix(1, nrow = n, ncol = n_sim)
  livre <- which(desvio > 0)
  if (length(livre) > 0) {
    for (j in seq_len(n_sim)) {
      sorteio <- rtnorm_vec(alvo$wd[livre], desvio[livre])
      mult[livre, j] <- (sorteio / alvo$wd[livre])^FEREZ[["b_wd"]]
    }
  }
  agb <- sweep(mult, 1, alvo$agb_kg, `*`)
  bind_rows(lapply(split(seq_len(n), alvo$parcela), function(idx) {
    total <- colSums(agb[idx, , drop = FALSE])
    q <- stats::quantile(total, c(0.025, 0.975), names = FALSE)
    tibble(parcela = alvo$parcela[idx[1]], agb_kg_wd_inf = q[1], agb_kg_wd_sup = q[2], agb_kg_wd_sd = stats::sd(total))
  }))
}

# ---- Parcelas, espécies ------------------------------------------------------

resumir_parcelas <- function(parcelas, individuos) {
  base <- parcelas %>% sf::st_drop_geometry()
  # Nomes próprios para o indivíduo: dentro do summarise, `agb_kg` passaria a ser
  # o total da parcela e quebraria as contas seguintes.
  sint <- individuos %>%
    rename(agb_i = agb_kg, sa_i = sa_m2) %>%
    group_by(parcela) %>%
    summarise(
      n_individuos = n(),
      n_mudas = sum(muda),
      n_plantio = sum(origem == "plantio"),
      n_regen = sum(origem == "regeneracao"),
      riqueza = n_distinct(especie),
      shannon = shannon_nat(as.numeric(table(especie))),
      especie_pred_n = predominante(nome_exibicao(especie), rep(1, n())),
      especie_pred_agb = predominante(nome_exibicao(especie), agb_i),
      agb_kg = sum(agb_i, na.rm = TRUE),
      agb_plantio_kg = sum(agb_i[origem == "plantio"], na.rm = TRUE),
      agb_regen_kg = sum(agb_i[origem == "regeneracao"], na.rm = TRUE),
      n_sem_biomassa = sum(!is.finite(agb_i)),
      sa_m2 = sum(sa_i, na.rm = TRUE),
      altura_media_m = mean(altura_m[is.finite(altura_m)]),
      prop_sadio = mean(fito == "sadio", na.rm = TRUE),
      prop_regen = mean(origem == "regeneracao"),
      .groups = "drop"
    )
  base %>%
    left_join(sint, by = "parcela") %>%
    mutate(
      tem_inventario = !is.na(n_individuos),
      sobrevivencia = n_plantio / MUDAS_PLANTADAS,
      agb_mg_ha = agb_kg / 1000 / (AREA_PARCELA_M2 / 10000),
      carbono_kg = agb_kg * FRACAO_CARBONO,
      carbono_mg_ha = agb_mg_ha * FRACAO_CARBONO,
      densidade_ha = n_individuos / (AREA_PARCELA_M2 / 10000)
    )
}

resumo_especies <- function(individuos) {
  individuos %>%
    group_by(especie) %>%
    summarise(
      n = n(),
      agb_kg = sum(agb_kg, na.rm = TRUE),
      wd = wd[1], wd_nivel = wd_nivel[1],
      .groups = "drop"
    ) %>%
    mutate(
      especie_nome = nome_exibicao(especie),
      pct_n = n / sum(n),
      pct_agb = if (sum(agb_kg) > 0) agb_kg / sum(agb_kg) else NA_real_
    ) %>%
    arrange(desc(n), especie)
}

# ---- Testes por tratamento ---------------------------------------------------

VARS_TRAT <- c(
  agb_plantio_kg = "Biomassa (kg/parcela)",
  sobrev_plantio = "Sobrevivência (% das 150 mudas)",
  riqueza_plantio = "Riqueza",
  altura_plantio_m = "Altura média (m)",
  sa_plantio_m2 = "Área seccional a 30 cm (m²/parcela)",
  prop_sadio_plantio = "Proporção sadia",
  prop_dap5_plantio = "Proporção com DAP ≥ 5 cm"
)
VARS_PROPORCAO <- c("sobrev_plantio", "prop_sadio_plantio", "prop_dap5_plantio")

metricas_plantio <- function(individuos, resumo) {
  trat <- resumo %>% filter(tem_inventario) %>% select(parcela, tratamento, talhao)
  individuos %>%
    filter(origem == "plantio") %>%
    group_by(parcela) %>%
    summarise(
      agb_plantio_kg = sum(agb_kg, na.rm = TRUE),
      n_plantio = n(),
      riqueza_plantio = n_distinct(especie),
      altura_plantio_m = mean(altura_m[is.finite(altura_m)]),
      sa_plantio_m2 = sum(sa_m2, na.rm = TRUE),
      n_sadio = sum(fito == "sadio", na.rm = TRUE),
      n_fito = sum(!is.na(fito)),
      prop_sadio_plantio = n_sadio / n_fito,
      prop_dap5_plantio = mean(is.finite(dap_max_cm) & dap_max_cm >= LIMIAR_DAP_CM),
      .groups = "drop"
    ) %>%
    right_join(trat, by = "parcela") %>%
    mutate(
      sobrev_plantio = n_plantio / MUDAS_PLANTADAS,
      tratamento = factor(tratamento, levels = niveis_tratamento(tratamento))
    )
}

poder_anova <- function(eps, n, k, alpha = 0.05) {
  if (!is.finite(eps) || !is.finite(n) || !is.finite(k) || eps < 0 || eps >= 1 || k < 2) return(NA_real_)
  df1 <- k - 1
  df2 <- n - k
  if (df2 < 1) return(NA_real_)
  fcrit <- stats::qf(1 - alpha, df1, df2)
  1 - stats::pf(fcrit, df1, df2, ncp = n * eps / (1 - eps))
}

eps_para_poder <- function(n, k, poder = 0.8, alpha = 0.05) {
  alvo <- function(e) poder_anova(e, n, k, alpha) - poder
  if (!is.finite(alvo(0.01)) || !is.finite(alvo(0.9))) return(NA_real_)
  if (alvo(0.01) > 0) return(0.01)
  if (alvo(0.9) < 0) return(NA_real_)
  stats::uniroot(alvo, c(0.01, 0.9))$root
}

cliff_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  mean(sign(outer(x, y, "-")))
}

descrever_tratamentos <- function(pp) {
  bind_rows(lapply(names(VARS_TRAT), function(v) {
    pp %>%
      group_by(tratamento) %>%
      summarise(
        mediana = median(.data[[v]], na.rm = TRUE),
        media = mean(.data[[v]], na.rm = TRUE),
        minimo = min(.data[[v]], na.rm = TRUE),
        maximo = max(.data[[v]], na.rm = TRUE),
        n = sum(is.finite(.data[[v]])),
        .groups = "drop"
      ) %>%
      mutate(variavel = v, variavel_nome = VARS_TRAT[[v]], posicao = rank(-mediana, ties.method = "min"))
  }))
}

prob_melhor <- function(pp, B = 2000L, seed = SEMENTE) {
  set.seed(seed)
  trats <- levels(pp$tratamento)
  idx <- split(seq_len(nrow(pp)), pp$tratamento)
  bind_rows(lapply(names(VARS_TRAT), function(v) {
    y <- pp[[v]]
    melhor <- setNames(numeric(length(trats)), trats)
    for (b in seq_len(B)) {
      medias <- vapply(idx[trats], function(i) mean(y[sample(i, length(i), replace = TRUE)], na.rm = TRUE), numeric(1))
      top <- which(medias == max(medias))
      melhor[top] <- melhor[top] + 1 / length(top)
    }
    tibble(variavel = v, variavel_nome = VARS_TRAT[[v]], tratamento = trats, prob_melhor = melhor / B)
  })) %>% mutate(tratamento = factor(tratamento, levels = trats))
}

um_contra_resto <- function(pp) {
  bind_rows(lapply(names(VARS_TRAT), function(v) {
    out <- bind_rows(lapply(levels(pp$tratamento), function(t) {
      x <- pp[[v]][pp$tratamento == t]
      y <- pp[[v]][pp$tratamento != t]
      tibble(
        tratamento = t, delta = cliff_delta(x, y),
        p = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
      )
    }))
    out$p_holm <- ajustar_holm(out$p)
    out %>% mutate(variavel = v, variavel_nome = VARS_TRAT[[v]], .before = 1)
  })) %>% mutate(tratamento = factor(tratamento, levels = levels(pp$tratamento)))
}

kw_tratamento <- function(pp) {
  out <- bind_rows(lapply(names(VARS_TRAT), function(v) {
    d <- pp %>% filter(is.finite(.data[[v]]))
    kt <- kruskal.test(d[[v]] ~ d$tratamento)
    n <- nrow(d)
    eps <- unname(kt$statistic) * (n + 1) / (n^2 - 1)
    tibble(
      variavel = v, variavel_nome = VARS_TRAT[[v]],
      H = unname(kt$statistic), gl = unname(kt$parameter), p = kt$p.value,
      efeito = eps, n = n,
      poder_aprox = poder_anova(eps, n, nlevels(d$tratamento))
    )
  }))
  out$p_holm <- ajustar_holm(out$p)
  k <- nlevels(pp$tratamento)
  n_familia <- sum(is.finite(out$p))
  out$eps2_80 <- mapply(eps_para_poder, out$n, k, MoreArgs = list(poder = 0.8, alpha = 0.05))
  out$eps2_80_holm <- mapply(eps_para_poder, out$n, k, MoreArgs = list(poder = 0.8, alpha = 0.05 / max(n_familia, 1)))
  out
}

modelos_tratamento <- function(pp, individuos, resumo) {
  linha <- function(resposta, modelo, estat, gl, p, nota) {
    tibble(resposta = resposta, modelo = modelo, estatistica = estat, gl = gl, p = p, nota = nota)
  }
  out <- list()
  q0 <- glm(cbind(n_plantio, MUDAS_PLANTADAS - n_plantio) ~ tratamento, family = quasibinomial, data = pp)
  a0 <- anova(q0, test = "F")
  out$surv <- linha(
    "Sobrevivência", "Quase-binomial sobre as 150 mudas, teste F",
    a0$F[2], a0$Df[2], a0$`Pr(>F)`[2],
    "Vivos de plantio sobre 150 plantadas. A sobredispersão entra no teste."
  )
  q1 <- glm(cbind(n_sadio, n_fito - n_sadio) ~ tratamento, family = quasibinomial, data = pp)
  aq <- anova(q1, test = "F")
  out$sad <- linha(
    "Proporção sadia", "Quase-binomial, teste F",
    aq$F[2], aq$Df[2], aq$`Pr(>F)`[2],
    "Usa o número de indivíduos de cada parcela. A sobredispersão entra no teste."
  )
  al <- anova(lm(log1p(agb_plantio_kg) ~ tratamento, data = pp))
  out$agb <- linha(
    "Biomassa (kg/parcela)", "Linear em log(1 + kg), teste F",
    al$`F value`[1], al$Df[1], al$`Pr(>F)`[1],
    "A escala log reduz o peso das parcelas mais pesadas."
  )
  trat <- resumo %>% filter(tem_inventario) %>% select(parcela, tratamento)
  hi <- individuos %>%
    filter(origem == "plantio", is.finite(altura_m), altura_m > 0) %>%
    left_join(trat, by = "parcela") %>%
    mutate(tratamento = factor(tratamento, levels = levels(pp$tratamento)))
  ajustar <- function(f1, f0) {
    b1 <- suppressMessages(lme4::lmer(f1, data = hi, REML = FALSE))
    b0 <- suppressMessages(lme4::lmer(f0, data = hi, REML = FALSE))
    list(teste = suppressMessages(anova(b0, b1)), fit = b1)
  }
  h1 <- ajustar(log(altura_m) ~ tratamento + (1 | parcela), log(altura_m) ~ 1 + (1 | parcela))
  out$h <- linha(
    "Altura do indivíduo", "Misto, parcela aleatória, razão de verossimilhança",
    h1$teste$Chisq[2], h1$teste$Df[2], h1$teste$`Pr(>Chisq)`[2],
    "Cada árvore entra, mas a parcela é a unidade do tratamento."
  )
  h2 <- ajustar(
    log(altura_m) ~ tratamento + (1 | especie) + (1 | parcela),
    log(altura_m) ~ 1 + (1 | especie) + (1 | parcela)
  )
  out$h_sp <- linha(
    "Altura do indivíduo, mesma espécie", "Misto, espécie e parcela aleatórias, razão de verossimilhança",
    h2$teste$Chisq[2], h2$teste$Df[2], h2$teste$`Pr(>Chisq)`[2],
    "Compara a altura dentro da espécie. Tira o efeito da mistura de espécies."
  )
  fe <- lme4::fixef(h2$fit)
  ef <- tibble(
    tratamento = levels(pp$tratamento),
    razao_altura = c(1, exp(fe[grep("^tratamento", names(fe))]))
  )
  tab <- bind_rows(out)
  tab$p_holm <- ajustar_holm(tab$p)
  list(tabela = tab, altura_mesma_especie = ef)
}

mde_tratamento <- function(pp) {
  n_rep <- max(2, round(stats::median(table(pp$tratamento))))
  bind_rows(lapply(c("agb_plantio_kg", "sobrev_plantio", "altura_plantio_m", "prop_sadio_plantio"), function(v) {
    s <- sigma(lm(pp[[v]] ~ pp$tratamento))
    k <- nlevels(pp$tratamento)
    tibble(
      variavel = v, variavel_nome = VARS_TRAT[[v]],
      media_geral = mean(pp[[v]], na.rm = TRUE), dp_dentro = s, n_rep = n_rep,
      dif_detectavel = power.t.test(n = n_rep, sd = s, power = 0.8, sig.level = 0.05)$delta,
      dif_detectavel_holm = power.t.test(n = n_rep, sd = s, power = 0.8, sig.level = 0.05 / choose(k, 2))$delta
    )
  })) %>% mutate(dif_pct_media = dif_detectavel / media_geral)
}

robustez_lider <- function(pp) {
  bind_rows(lapply(names(VARS_TRAT), function(v) {
    med <- tapply(pp[[v]], pp$tratamento, median, na.rm = TRUE)
    ordem <- sort(med, decreasing = TRUE)
    lider <- names(ordem)[1]
    sub <- pp[pp$tratamento == lider, ]
    retirar <- which.max(sub[[v]])
    resto <- pp[!(pp$tratamento == lider & pp$parcela == sub$parcela[retirar]), ]
    med2 <- tapply(resto[[v]], resto$tratamento, median, na.rm = TRUE)
    fica <- vapply(seq_len(nrow(pp)), function(i) {
      m <- tapply(pp[[v]][-i], pp$tratamento[-i], median, na.rm = TRUE)
      names(m)[which.max(m)] == lider
    }, logical(1))
    tibble(
      variavel = v, variavel_nome = VARS_TRAT[[v]], lider = lider,
      mediana_lider = ordem[[1]], segundo = names(ordem)[2], mediana_segundo = ordem[[2]],
      parcela_retirada = sub$parcela[retirar],
      posicao_sem_melhor_parcela = rank(-med2, ties.method = "min")[[lider]],
      prop_jackknife_lider = mean(fica)
    )
  }))
}

permanova_grupo <- function(individuos, resumo, coluna, min_n = 2) {
  mapa <- resumo %>% filter(tem_inventario, !is.na(.data[[coluna]])) %>% mutate(grupo_f = as.character(.data[[coluna]]))
  n_g <- mapa %>% count(grupo_f)
  mapa <- mapa %>% filter(grupo_f %in% n_g$grupo_f[n_g$n >= min_n])
  vazio <- list(r2 = NA_real_, p = NA_real_, n = 0L, n_grupos = 0L)
  comp <- individuos %>% filter(parcela %in% mapa$parcela) %>% count(parcela, especie)
  if (nrow(comp) == 0 || n_distinct(mapa$grupo_f) < 2) return(vazio)
  mat <- comp %>%
    pivot_wider(names_from = especie, values_from = n, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- as.character(mat$parcela)
  mat$parcela <- NULL
  mat <- as.matrix(mat)
  grupos <- mapa$grupo_f[match(as.integer(rownames(mat)), mapa$parcela)]
  set.seed(SEMENTE)
  fit <- as.data.frame(adonis2(mat ~ factor(grupos), method = "hellinger", permutations = 999))
  list(r2 = fit$R2[1], p = fit$`Pr(>F)`[1], n = nrow(mat), n_grupos = n_distinct(grupos))
}

# Talhão: só entra se a fazenda tiver dois ou mais talhões com repetição.
analise_talhoes <- function(resumo, individuos) {
  inv <- resumo %>% filter(tem_inventario, !is.na(talhao))
  n_t <- inv %>% count(talhao) %>% filter(n >= 2)
  if (nrow(n_t) < 2) return(NULL)
  vars <- c(agb_kg = "Biomassa (kg/parcela)", n_individuos = "Indivíduos vivos", riqueza = "Riqueza", altura_media_m = "Altura média (m)")
  kw <- bind_rows(lapply(names(vars), function(v) {
    d <- inv %>% filter(talhao %in% n_t$talhao, is.finite(.data[[v]]))
    kt <- kruskal.test(d[[v]] ~ factor(d$talhao))
    tibble(variavel = v, variavel_nome = vars[[v]], H = unname(kt$statistic), gl = unname(kt$parameter), p = kt$p.value, n = nrow(d))
  }))
  kw$p_holm <- ajustar_holm(kw$p)
  list(kw = kw, permanova = permanova_grupo(individuos, resumo, "talhao"), talhoes = n_t)
}

# ---- Execução de uma fazenda -------------------------------------------------

analisar_fazenda <- function(raiz, cod, saida = file.path(raiz, "saida", cod)) {
  dir.create(saida, recursive = TRUE, showWarnings = FALSE)
  bruto <- ler_fazenda(raiz, cod)
  fustes <- bruto$fustes
  individuos <- agregar_individuos(fustes)
  wd <- densidade_madeira(individuos, cod)
  individuos <- calcular_biomassa(individuos, wd)
  inc <- propagar_wd(individuos)
  resumo <- resumir_parcelas(bruto$parcelas, individuos) %>%
    left_join(inc, by = "parcela") %>%
    mutate(
      agb_kg_wd_inf = ifelse(is.na(agb_kg_wd_inf) & is.finite(agb_kg), agb_kg, agb_kg_wd_inf),
      agb_kg_wd_sup = ifelse(is.na(agb_kg_wd_sup) & is.finite(agb_kg), agb_kg, agb_kg_wd_sup)
    )
  tem_regen <- any(!is.na(fustes$regeneracao))
  inv <- resumo %>% filter(tem_inventario)
  stopifnot(all(abs(inv$agb_plantio_kg + inv$agb_regen_kg - inv$agb_kg) < 1e-6))
  stopifnot(abs(sum(inv$agb_kg) - sum(individuos$agb_kg, na.rm = TRUE)) < 1e-6)
  stopifnot(sum(inv$n_sem_biomassa) == sum(!is.finite(individuos$agb_kg)))
  # Fustes cujo D30 não entrou na área seccional: ausente ou implausível.
  d30_suspeito <- fustes %>%
    filter(!(is.finite(d30_cm) & d30_cm > 0 & d30_cm <= D30_MAX_PLAUSIVEL_CM)) %>%
    mutate(motivo = ifelse(is.finite(d30_cm) & d30_cm > D30_MAX_PLAUSIVEL_CM, "D30 implausível", "D30 ausente")) %>%
    select(codigo_fazenda, parcela, id_ind, especie, d30_cm, dap_cm, altura_m, motivo)

  pp <- metricas_plantio(individuos, resumo)
  kw <- kw_tratamento(pp)
  at <- list(
    parcelas = pp,
    descricao = descrever_tratamentos(pp),
    prob_melhor = prob_melhor(pp),
    um_contra_resto = um_contra_resto(pp),
    kw = kw,
    modelos = modelos_tratamento(pp, individuos, resumo),
    mde = mde_tratamento(pp),
    robustez = robustez_lider(pp)
  )
  at$posicao_media <- at$descricao %>%
    group_by(tratamento) %>%
    summarise(posicao_media = mean(posicao), n_primeiro = sum(posicao == 1), .groups = "drop") %>%
    arrange(posicao_media)

  perm_trat <- permanova_grupo(individuos, resumo, "tratamento")
  perm_trat_plantio <- permanova_grupo(individuos %>% filter(origem == "plantio"), resumo, "tratamento")
  talhoes <- analise_talhoes(resumo, individuos)

  pol <- bruto$parcelas %>%
    left_join(
      resumo %>% select(parcela, agb_kg, agb_kg_wd_inf, agb_kg_wd_sup, agb_mg_ha, n_individuos, sobrevivencia, riqueza, tem_inventario),
      by = "parcela"
    )

  readr::write_csv(resumo, file.path(saida, "resumo_parcelas.csv"), na = "")
  readr::write_csv(
    individuos %>% select(
      codigo_fazenda, parcela, id_ind, linha, especie, n_especies_id, n_fustes, sa_m2, dap_max_cm,
      altura_m, muda, fito, origem, wd, wd_sd, wd_nivel, agb_kg, carbono_kg
    ),
    file.path(saida, "individuos.csv"), na = ""
  )
  readr::write_csv(resumo_especies(individuos), file.path(saida, "especies.csv"), na = "")
  readr::write_csv(wd, file.path(saida, "densidade_madeira.csv"), na = "")
  readr::write_csv(pp, file.path(saida, "tratamentos_parcelas.csv"), na = "")
  readr::write_csv(at$descricao, file.path(saida, "tratamentos_descricao.csv"), na = "")
  readr::write_csv(at$prob_melhor, file.path(saida, "tratamentos_prob_melhor.csv"), na = "")
  readr::write_csv(at$um_contra_resto, file.path(saida, "tratamentos_um_contra_resto.csv"), na = "")
  readr::write_csv(at$kw, file.path(saida, "tratamentos_kruskal.csv"), na = "")
  readr::write_csv(at$modelos$tabela, file.path(saida, "tratamentos_modelos.csv"), na = "")
  readr::write_csv(at$mde, file.path(saida, "tratamentos_diferenca_detectavel.csv"), na = "")
  readr::write_csv(at$robustez, file.path(saida, "tratamentos_robustez_lider.csv"), na = "")
  readr::write_csv(d30_suspeito, file.path(saida, "d30_suspeito.csv"), na = "")
  gpkg <- file.path(saida, sprintf("parcelas_%s.gpkg", cod))
  if (file.exists(gpkg)) unlink(gpkg)
  sf::st_write(pol, gpkg, layer = "parcelas", quiet = TRUE, delete_dsn = TRUE)

  res <- list(
    cod = cod, nome = NOMES_FAZENDA[[cod]], arquivos = bruto$arquivos,
    parcelas = pol, fustes = fustes, individuos = individuos, wd = wd,
    resumo = resumo, especies = resumo_especies(individuos),
    tratamentos = at, kw = at$kw,
    permanova_tratamento = perm_trat, permanova_tratamento_plantio = perm_trat_plantio,
    talhoes = talhoes, tem_regeneracao = tem_regen, d30_suspeito = d30_suspeito
  )
  saveRDS(res, file.path(saida, "resultados.rds"))
  res
}

# ---- Gráficos (os mesmos para BH e ER) --------------------------------------

tema_relatorio <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "#444444", size = 9),
      legend.position = "bottom",
      axis.title = element_text(size = 10),
      strip.text = element_text(face = "bold")
    )
}

grafico_mapa <- function(res) {
  pol <- res$parcelas
  pol$rotulo <- paste0("P", pol$parcela, "\n", pol$tratamento)
  lim <- range(pol$agb_kg, na.rm = TRUE)
  pol$texto_claro <- !is.na(pol$agb_kg) & (pol$agb_kg - lim[1]) / diff(lim) < 0.55
  cen <- sf::st_as_sf(sf::st_drop_geometry(pol) %>% select(parcela, rotulo, texto_claro, lon, lat), coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  caixa <- sf::st_bbox(res$parcelas)
  mx <- 0.1 * (caixa["xmax"] - caixa["xmin"])
  my <- 0.1 * (caixa["ymax"] - caixa["ymin"])
  ggplot() +
    geom_sf(data = res$parcelas, aes(fill = agb_kg), color = "white", linewidth = 0.3) +
    geom_sf_text(data = cen, aes(label = rotulo, color = texto_claro), size = 2.2, fontface = "bold", lineheight = 0.9, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_viridis_c(
      option = "cividis", name = "Biomassa (kg/parcela)",
      labels = label_number(decimal.mark = ",", big.mark = "."), na.value = "grey80"
    ) +
    coord_sf(
      xlim = c(caixa["xmin"] - mx, caixa["xmax"] + mx),
      ylim = c(caixa["ymin"] - my, caixa["ymax"] + my), expand = FALSE
    ) +
    labs(
      title = paste("Parcelas de", res$nome),
      subtitle = "Posição e tratamento do GeoPackage. Cinza: parcela sem inventário.",
      x = NULL, y = NULL
    ) +
    tema_relatorio() +
    theme(legend.position = "right", axis.text = element_text(size = 7))
}

grafico_biomassa <- function(res) {
  d <- res$resumo %>% filter(tem_inventario)
  niveis <- niveis_tratamento(d$tratamento)
  d <- d %>%
    mutate(tratamento = factor(tratamento, levels = niveis)) %>%
    arrange(tratamento, parcela) %>%
    mutate(rotulo = factor(paste0("P", parcela), levels = unique(paste0("P", parcela))))
  med <- d %>% group_by(tratamento) %>% summarise(mediana = median(agb_kg), .groups = "drop")
  ggplot(d, aes(rotulo, agb_kg, fill = tratamento)) +
    geom_col(width = 0.8) +
    geom_hline(data = med, aes(yintercept = mediana), linetype = "dashed", linewidth = 0.4, color = "#222222") +
    geom_errorbar(aes(ymin = agb_kg_wd_inf, ymax = agb_kg_wd_sup), width = 0.25, linewidth = 0.3, color = "#222222") +
    facet_grid(. ~ tratamento, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(values = PAL_TRAT[niveis], guide = "none") +
    scale_y_continuous(labels = label_number(decimal.mark = ",", big.mark = "."), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "Biomassa aérea por tratamento",
      subtitle = "Parcelas do mesmo tratamento, lado a lado. Tracejado: mediana do tratamento. Traço: incerteza da densidade da madeira.",
      x = NULL, y = "Biomassa aérea (kg/parcela)"
    ) +
    tema_relatorio() +
    theme(
      axis.text.x = element_text(size = 6.5, angle = 90, vjust = 0.5, hjust = 1),
      strip.placement = "outside",
      panel.spacing.x = unit(0.25, "lines"),
      strip.text = element_text(face = "bold", size = 10)
    )
}

grafico_especies <- function(res) {
  top <- res$especies %>% slice_head(n = 12) %>% mutate(especie_nome = factor(especie_nome, levels = rev(especie_nome)))
  a <- ggplot(top, aes(n, especie_nome)) +
    geom_col(fill = "#0072B2", width = 0.72) +
    labs(title = "Abundância", x = "Indivíduos", y = NULL) +
    tema_relatorio()
  b <- ggplot(top, aes(agb_kg, especie_nome)) +
    geom_col(fill = "#E69F00", width = 0.72) +
    scale_x_continuous(labels = label_number(decimal.mark = ",", big.mark = ".")) +
    labs(title = "Biomassa", x = "kg", y = NULL) +
    tema_relatorio()
  a + b +
    plot_annotation(
      title = "Doze espécies mais abundantes",
      subtitle = "A ordem é a da abundância. A biomassa não segue a mesma ordem.",
      theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, color = "#444444"))
    )
}

grafico_tratamentos <- function(res) {
  pp <- res$tratamentos$parcelas
  d <- pp %>%
    mutate(sobrev_plantio = 100 * sobrev_plantio) %>%
    select(tratamento, agb_plantio_kg, sobrev_plantio, riqueza_plantio, altura_plantio_m, prop_sadio_plantio) %>%
    pivot_longer(-tratamento, names_to = "variavel", values_to = "valor") %>%
    mutate(
      variavel = recode(
        variavel,
        agb_plantio_kg = "Biomassa (kg)", sobrev_plantio = "Sobrevivência (%)",
        riqueza_plantio = "Riqueza", altura_plantio_m = "Altura média (m)", prop_sadio_plantio = "Proporção sadia"
      ),
      variavel = factor(variavel, levels = c("Biomassa (kg)", "Sobrevivência (%)", "Riqueza", "Altura média (m)", "Proporção sadia"))
    )
  ggplot(d, aes(tratamento, valor)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, fill = "#56B4E9", alpha = 0.85) +
    geom_jitter(width = 0.08, size = 1.5, alpha = 0.85) +
    facet_wrap(~variavel, scales = "free_y", nrow = 1) +
    labs(
      title = "As parcelas de cada tratamento, só o plantio",
      subtitle = paste0("Cada ponto é uma parcela de ", res$nome, ". T1 a T7 são o fator replicado."),
      x = "Tratamento", y = NULL
    ) +
    tema_relatorio()
}

grafico_ranking_tratamentos <- function(res) {
  d <- res$tratamentos$descricao %>%
    mutate(
      variavel_nome = factor(variavel_nome, levels = rev(unname(VARS_TRAT))),
      rotulo = dplyr::case_when(
        variavel %in% VARS_PROPORCAO ~ paste0(fmt_num(100 * mediana, 0), "%"),
        variavel %in% c("altura_plantio_m", "sa_plantio_m2") ~ fmt_num(mediana, 2),
        TRUE ~ fmt_num(mediana, 0)
      )
    )
  ggplot(d, aes(tratamento, variavel_nome, fill = posicao)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = rotulo), size = 3.2) +
    scale_fill_gradient(low = "#1b7837", high = "#f7f7f7", name = "Posição\n(1 = maior)", breaks = 1:7) +
    scale_x_discrete(position = "top") +
    labs(
      title = "Posição de cada tratamento, só o plantio",
      subtitle = "Número: mediana das parcelas do tratamento. Cor: posição entre os sete. Uma posição não é uma diferença testada.",
      x = NULL, y = NULL
    ) +
    tema_relatorio() +
    theme(panel.grid = element_blank(), legend.position = "right")
}

grafico_prob_melhor <- function(res) {
  d <- res$tratamentos$prob_melhor %>%
    mutate(variavel_nome = factor(variavel_nome, levels = rev(unname(VARS_TRAT))))
  ggplot(d, aes(tratamento, variavel_nome, fill = prob_melhor)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(fmt_num(100 * prob_melhor, 0), "%")), size = 3.2) +
    scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", limits = c(0, 1), labels = label_percent(), name = "Probabilidade") +
    scale_x_discrete(position = "top") +
    labs(
      title = "Com que frequência cada tratamento fica em primeiro",
      subtitle = "2.000 reamostragens das parcelas de cada tratamento. Com 14% (1 em 7) ninguém se destaca.",
      x = NULL, y = NULL
    ) +
    tema_relatorio() +
    theme(panel.grid = element_blank(), legend.position = "right")
}

grafico_talhoes <- function(res) {
  tal <- res$talhoes$talhoes$talhao
  d <- res$resumo %>%
    filter(tem_inventario, talhao %in% tal) %>%
    select(talhao, agb_kg, n_individuos, riqueza, altura_media_m) %>%
    pivot_longer(-talhao, names_to = "variavel", values_to = "valor") %>%
    mutate(
      talhao_f = factor(talhao),
      variavel = recode(
        variavel,
        agb_kg = "Biomassa (kg)", n_individuos = "Indivíduos vivos",
        riqueza = "Riqueza", altura_media_m = "Altura média (m)"
      )
    )
  pal <- setNames(hue_pal()(length(tal)), as.character(tal))
  pal[names(pal) == "558"] <- "#0072B2"
  pal[names(pal) == "559"] <- "#E69F00"
  pal[names(pal) == "560"] <- "#009E73"
  ggplot(d, aes(talhao_f, valor, fill = talhao_f)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.08, size = 1.4, alpha = 0.8) +
    facet_wrap(~variavel, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      title = "Parcelas dentro de cada talhão",
      subtitle = "Cada ponto é uma parcela. O teste usa só os talhões com pelo menos duas parcelas inventariadas.",
      x = "Talhão", y = NULL
    ) +
    tema_relatorio()
}

# ---- Textos dos slides (gerados dos resultados, nunca fixos) -----------------

NOMES_CURTOS <- c(
  agb_plantio_kg = "biomassa", sobrev_plantio = "sobrevivência", riqueza_plantio = "riqueza",
  altura_plantio_m = "altura", sa_plantio_m2 = "área seccional",
  prop_sadio_plantio = "proporção sadia", prop_dap5_plantio = "proporção com DAP ≥ 5 cm"
)

lista_pt <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return("")
  if (length(x) == 1) return(x)
  paste0(paste(x[-length(x)], collapse = ", "), " e ", x[length(x)])
}

sig_texto <- function(tab) {
  nomes <- tab$variavel_nome[is.finite(tab$p_holm) & tab$p_holm < 0.05]
  if (length(nomes) == 0) "nenhuma variável passa de Holm" else paste(nomes, collapse = "; ")
}

leitura_perm <- function(p) {
  if (!is.finite(p)) return("sem teste")
  if (p < 0.05) "a composição de espécies difere entre os grupos" else "a composição de espécies não se separa entre os grupos"
}

# "T6: biomassa, sobrevivência e altura. T4: área seccional."
texto_lideres <- function(at) {
  d <- at$descricao %>%
    filter(posicao == 1) %>%
    mutate(nome = NOMES_CURTOS[variavel]) %>%
    group_by(tratamento) %>%
    summarise(vars = lista_pt(nome), n = n(), .groups = "drop") %>%
    arrange(desc(n), tratamento)
  paste0(d$tratamento, ": ", d$vars, ".", collapse = " ")
}

texto_firmeza <- function(at, limite = 0.95) {
  rb <- at$robustez %>%
    mutate(nome = NOMES_CURTOS[variavel], rotulo = paste0(nome, " (", lider, ")"))
  est <- rb %>% filter(prop_jackknife_lider >= limite)
  ins <- rb %>% filter(prop_jackknife_lider < limite)
  a <- if (nrow(est) > 0) {
    paste0("O primeiro lugar resiste a tirar qualquer parcela em ", lista_pt(est$rotulo), ".")
  } else {
    "Em nenhuma variável o primeiro lugar resiste a tirar qualquer parcela."
  }
  b <- if (nrow(ins) > 0) {
    paste0(
      "Troca conforme a parcela retirada em ",
      lista_pt(paste0(ins$rotulo, ", em ", fmt_num(100 * ins$prop_jackknife_lider, 0), "% das retiradas")), "."
    )
  } else {
    ""
  }
  bio <- rb %>% filter(variavel == "agb_plantio_kg")
  c_ <- if (bio$posicao_sem_melhor_parcela == 1) {
    paste0("Na biomassa, ", bio$lider, " segue em primeiro mesmo sem a parcela ", bio$parcela_retirada, ".")
  } else {
    paste0("Na biomassa, sem a parcela ", bio$parcela_retirada, ", ", bio$lider, " cai para a posição ", bio$posicao_sem_melhor_parcela, ".")
  }
  str_squish(paste(a, b, c_))
}

# ---- Execução direta ---------------------------------------------------------

if (sys.nframe() == 0) {
  candidatos <- c(Sys.getenv("MAXIC_RAIZ", unset = getwd()), file.path(getwd(), "capitulo-03"), dirname(getwd()))
  raiz <- normalizePath(candidatos[dir.exists(file.path(candidatos, "0-raw"))][1])
  for (cod in c("BH", "ER")) {
    r <- analisar_fazenda(raiz, cod)
    cat("\n====", r$nome, "====\n")
    print(r$resumo %>% summarise(
      parcelas = sum(tem_inventario), individuos = sum(n_individuos, na.rm = TRUE),
      agb_Mg = sum(agb_kg, na.rm = TRUE) / 1000, mediana_kg = median(agb_kg, na.rm = TRUE),
      mediana_Mg_ha = median(agb_mg_ha, na.rm = TRUE), sem_biomassa = sum(n_sem_biomassa, na.rm = TRUE)
    ))
  }
}
