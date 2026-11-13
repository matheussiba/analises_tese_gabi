# Inventário MaxiC: pontos das parcelas, biomassa e testes.
# A planilha qualitativa fornece os quatro vértices. A planilha MaxiC
# fornece os fustes. Em Belo Horizonte cada tratamento é uma parcela e o
# fator com repetição é o talhão. Em Entre Rios os tratamentos T1–T7 se repetem.

user_lib <- Sys.getenv(
  "R_LIBS_USER",
  unset = path.expand("~/R/x86_64-pc-linux-gnu-library/4.3")
)
.libPaths(c(user_lib, .libPaths()))

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

# Carbono da biomassa aérea: fração default do IPCC (2006), não medida em laboratório.
FRACAO_CARBONO <- 0.47
# Equação de Chave et al. (2014), eq. 4, calibrada sobretudo para D >= 5 cm.
LIMIAR_DAP_CM <- 5
LIMIAR_INCLUSAO_CM <- 1

ascii_lower <- function(x) {
  y <- stringi::stri_trans_general(trimws(as.character(x)), "Latin-ASCII")
  tolower(y)
}

data_excel <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "N/A")] <- NA_character_
  num <- suppressWarnings(as.numeric(y))
  out <- as.Date(rep(NA_real_, length(y)), origin = "1970-01-01")
  serial <- is.finite(num) & num > 20000 & num < 80000
  out[serial] <- as.Date(num[serial], origin = "1899-12-30")
  resto <- !serial & !is.na(y)
  if (any(resto)) {
    iso <- suppressWarnings(as.Date(y[resto]))
    out[resto] <- iso
  }
  out
}

numero_br <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "N/A", "na", "Na", "null")] <- NA_character_
  y <- gsub("[[:space:]]", "", y)
  y <- gsub(",", ".", y)
  suppressWarnings(as.numeric(y))
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

codigo_fazenda <- function(x) {
  y <- ascii_lower(x)
  dplyr::case_when(
    str_detect(y, "belo\\s*horizonte|^bh$") ~ "BH",
    str_detect(y, "entre\\s*rios|^er$") ~ "ER",
    TRUE ~ toupper(substr(trimws(as.character(x)), 1, 12))
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

fmt_num <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "—",
    formatC(x, format = "f", digits = digits, decimal.mark = ",", big.mark = "")
  )
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "—",
    ifelse(p < 0.001, "< 0,001", formatC(p, format = "f", digits = 3, decimal.mark = ","))
  )
}

id_parcela <- function(x) {
  y <- ascii_lower(x)
  y <- str_replace(y, "^p", "")
  as.integer(round(numero_br(y)))
}

epsg_utm_sul <- function(lon) {
  zona <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  zona <- min(max(as.integer(zona), 18L), 25L)
  31960 + zona
}

parse_latlon <- function(x) {
  y <- as.character(x)
  partes <- str_split(str_replace_all(y, ";", ","), ",")
  lat <- vapply(partes, function(p) {
    if (length(p) < 2) return(NA_real_)
    numero_br(p[1])
  }, numeric(1))
  lon <- vapply(partes, function(p) {
    if (length(p) < 2) return(NA_real_)
    numero_br(p[2])
  }, numeric(1))
  tibble(lat = lat, lon = lon)
}

anel_utm <- function(lat, lon) {
  ok <- is.finite(lat) & is.finite(lon)
  if (sum(ok) < 3) {
    return(NULL)
  }
  epsg <- epsg_utm_sul(lon[ok])
  pts <- st_as_sf(tibble(lat = lat[ok], lon = lon[ok]), coords = c("lon", "lat"), crs = 4326)
  utm <- st_transform(pts, epsg)
  xy <- st_coordinates(utm)[, 1:2, drop = FALSE]
  if (nrow(xy) < 3 || any(duplicated(round(xy, 2)))) {
    # vértices repetidos no mesmo canto: ainda tenta o anel com os pontos distintos
    xy <- xy[!duplicated(round(xy, 2)), , drop = FALSE]
  }
  if (nrow(xy) < 3) {
    return(NULL)
  }
  centro <- colMeans(xy)
  ang <- atan2(xy[, 2] - centro[2], xy[, 1] - centro[1])
  xy <- xy[order(ang), , drop = FALSE]
  lado <- function(a) {
    b <- rbind(a[-1, , drop = FALSE], a[1, , drop = FALSE])
    sqrt(rowSums((a - b)^2))
  }
  geom <- st_sfc(st_polygon(list(rbind(xy, xy[1, , drop = FALSE]))), crs = epsg)
  tibble(
    area_m2 = as.numeric(st_area(geom)),
    lado_max_m = max(lado(xy)),
    lado_min_m = min(lado(xy)),
    geometry = geom
  )
}

ler_planilhas <- function(pasta_dados) {
  arquivos <- list.files(pasta_dados, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  arquivos <- arquivos[!str_detect(basename(arquivos), "^~\\$")]
  if (length(arquivos) == 0) {
    stop("Nenhum .xlsx em ", pasta_dados)
  }
  quali <- list()
  quant <- list()
  for (arq in arquivos) {
    abas <- excel_sheets(arq)
    aba_q <- abas[str_detect(tolower(abas), "qualit")][1]
    if (is.na(aba_q)) aba_q <- abas[1]
    aba_n <- abas[str_detect(tolower(abas), "maxic|quant")][1]
    q <- read_excel(arq, sheet = aba_q, col_types = "text") %>%
      mutate(arquivo = basename(arq))
    quali[[length(quali) + 1]] <- q
    if (!is.na(aba_n)) {
      n <- read_excel(arq, sheet = aba_n, col_types = "text") %>%
        mutate(arquivo = basename(arq))
      quant[[length(quant) + 1]] <- n
    }
  }
  list(quali = bind_rows(quali), quant = bind_rows(quant), arquivos = basename(arquivos))
}

preparar_qualitativa <- function(q) {
  nomes <- names(q)
  padrao <- function(alvo) {
    hit <- nomes[ascii_lower(nomes) == ascii_lower(alvo)]
    if (length(hit) == 0) stop("Coluna ausente na planilha qualitativa: ", alvo)
    hit[1]
  }
  coords <- bind_rows(lapply(c("P1", "P2", "P3", "P4"), function(p) {
    ll <- parse_latlon(q[[padrao(p)]])
    ll$ponto <- p
    foto_nome <- paste0("F", substr(p, 2, 2))
    foto_hit <- nomes[ascii_lower(nomes) == ascii_lower(foto_nome)]
    ll$foto <- if (length(foto_hit) == 0) NA_character_ else as.character(q[[foto_hit[1]]])
    ll$i <- seq_len(nrow(q))
    ll
  }))
  base <- tibble(
    i = seq_len(nrow(q)),
    fazenda = str_squish(as.character(q[[padrao("FAZENDA")]])),
    codigo_fazenda = codigo_fazenda(q[[padrao("FAZENDA")]]),
    protocolo = as.character(q[[padrao("PROTOCOLO")]]),
    nome_exp = as.character(q[[padrao("NOME_EXP")]]),
    talhao = as.integer(numero_br(q[[padrao("TALHAO")]])),
    data_coleta = data_excel(q[[padrao("DATA")]]),
    responsavel = str_squish(as.character(q[[padrao("RESPONSAVEL_TECNICO")]])),
    parcela = id_parcela(q[[padrao("ID_PARCELA")]]),
    tratamento = str_squish(as.character(q[[padrao("TRATAMENTO")]])),
    solo = str_squish(as.character(q[[padrao("SOLO")]])),
    declividade = str_squish(as.character(q[[padrao("DECLIVIDADE")]])),
    erosao = sim_nao(q[[padrao("PRES_EROSAO")]]),
    cat_erosao = as.character(q[[padrao("CAT_EROSAO")]]),
    formigas = sim_nao(q[[padrao("PRES_FORMIGAS")]]),
    nivel_formigas = as.character(q[[padrao("NIVEL_ATAQUE_FORM")]]),
    cupim = sim_nao(q[[padrao("PRES_CUPIM")]]),
    cupim_atividade = as.character(q[[padrao("CUPIM_ATIVIDADE")]]),
    n_cupinzeiros = numero_br(q[[padrao("NUMERO_CUPINZEIROS")]]),
    guandu = sim_nao(q[[padrao("PRES_GUANDU")]]),
    perc_folha_larga = percentual(q[[padrao("PERC_FOLHA_LARGA")]]),
    perc_folha_estreita = percentual(q[[padrao("PERC_FOLHA_ESTREITA")]]),
    animais = sim_nao(q[[padrao("IND_ANIMAIS_INV")]]),
    fogo = sim_nao(q[[padrao("IND_FOGO")]]),
    obs = as.character(q[[padrao("OBS")]])
  )
  base <- base %>% filter(!is.na(parcela), !is.na(fazenda))
  if (any(duplicated(base[c("codigo_fazenda", "parcela")]))) {
    stop("Há parcela repetida na mesma fazenda. Confira os arquivos em dados/.")
  }
  pontos <- coords %>%
    inner_join(base, by = "i") %>%
    select(-i) %>%
    mutate(foto = na_if(foto, "NA"))
  pontos
}

montar_geometrias <- function(pontos) {
  chaves <- pontos %>%
    distinct(codigo_fazenda, parcela) %>%
    mutate(id = row_number())
  poligonos <- vector("list", nrow(chaves))
  for (k in seq_len(nrow(chaves))) {
    sub <- pontos %>%
      filter(codigo_fazenda == chaves$codigo_fazenda[k], parcela == chaves$parcela[k]) %>%
      arrange(ponto)
    anel <- anel_utm(sub$lat, sub$lon)
    if (is.null(anel)) {
      poligonos[[k]] <- tibble(
        codigo_fazenda = chaves$codigo_fazenda[k],
        parcela = chaves$parcela[k],
        area_m2 = NA_real_,
        lado_max_m = NA_real_,
        lado_min_m = NA_real_,
        geometry = st_sfc(st_polygon(), crs = 4326)[1]
      )
      next
    }
    geom4326 <- st_transform(anel$geometry, 4326)
    poligonos[[k]] <- tibble(
      codigo_fazenda = chaves$codigo_fazenda[k],
      parcela = chaves$parcela[k],
      area_m2 = anel$area_m2,
      lado_max_m = anel$lado_max_m,
      lado_min_m = anel$lado_min_m,
      geometry = geom4326
    )
  }
  pol <- bind_rows(poligonos) %>% st_as_sf(sf_column_name = "geometry", crs = 4326)
  assinatura <- pontos %>%
    group_by(codigo_fazenda, parcela) %>%
    summarise(
      assinatura = paste(
        sprintf("%.5f %.5f", sort(lat), sort(lon)),
        collapse = "|"
      ),
      .groups = "drop"
    )
  # sort(lat) com sort(lon) separado não emparelha o vértice. Assinar o par ordenado.
  assinatura <- pontos %>%
    mutate(par = sprintf("%.5f,%.5f", round(lat, 5), round(lon, 5))) %>%
    group_by(codigo_fazenda, parcela) %>%
    summarise(assinatura = paste(sort(par), collapse = "|"), .groups = "drop")
  repetidas <- assinatura %>%
    group_by(codigo_fazenda, assinatura) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(chave = paste(codigo_fazenda, assinatura, sep = "|")) %>%
    pull(chave)
  attr <- pontos %>%
    distinct(
      codigo_fazenda, fazenda, protocolo, nome_exp, talhao, parcela, tratamento,
      data_coleta, responsavel, solo, declividade, erosao, formigas, cupim,
      n_cupinzeiros, perc_folha_larga, perc_folha_estreita, obs
    )
  pol <- pol %>%
    left_join(assinatura, by = c("codigo_fazenda", "parcela")) %>%
    mutate(
      geometria_duplicada = paste(codigo_fazenda, assinatura, sep = "|") %in% repetidas,
      geometria_valida = is.finite(area_m2) &
        area_m2 >= 400 & area_m2 <= 1600 &
        lado_max_m <= 80 & lado_min_m >= 8
    ) %>%
    left_join(attr, by = c("codigo_fazenda", "parcela"))
  list(poligonos = pol, assinatura = assinatura)
}

preparar_fustes <- function(quant, pontos) {
  if (nrow(quant) == 0) {
    return(tibble())
  }
  nomes <- names(quant)
  col <- function(alvo) {
    hit <- nomes[ascii_lower(nomes) == ascii_lower(alvo)]
    if (length(hit) == 0) stop("Coluna ausente na planilha quantitativa: ", alvo)
    hit[1]
  }
  col_opc <- function(alvo) {
    hit <- nomes[ascii_lower(nomes) == ascii_lower(alvo)]
    if (length(hit) == 0) NULL else hit[1]
  }
  regen_col <- col_opc("REGENERACAO")
  fustes <- tibble(
    fazenda = str_squish(quant[[col("FAZENDA")]]),
    codigo_fazenda = codigo_fazenda(quant[[col("FAZENDA")]]),
    protocolo = quant[[col("PROTOCOLO")]],
    parcela = id_parcela(quant[[col("ID_PARCELA")]]),
    linha = as.integer(numero_br(quant[[col("N_LINHA")]])),
    id_ind = {
      y <- str_squish(as.character(quant[[col("ID_IND")]]))
      y[y %in% c("", "NA", "N/A", "na")] <- NA_character_
      y
    },
    fuste = quant[[col("FUSTE")]],
    especie = str_squish(quant[[col("ESPECIE")]]),
    muda_campo = sim_nao(quant[[col("MUDA")]]),
    d30_cm = numero_br(quant[[col("D30_CM")]]),
    dap_cm = numero_br(quant[[col("DAP_CM")]]),
    altura_m = numero_br(quant[[col("ALTURA_TOTAL_M")]]),
    dc_m = numero_br(quant[[col("DC1_M")]]),
    fito = classe_fito(quant[[col("ESTADO_FITOSSANITARIO")]]),
    regeneracao = if (is.null(regen_col)) NA_character_ else sim_nao(quant[[regen_col]]),
    obs = quant[[col("OBS")]]
  ) %>%
    filter(!is.na(parcela), !is.na(id_ind), especie != "")
  fustes
}

agregar_individuos <- function(fustes) {
  gravidade <- c(sadio = 1L, inicial = 2L, avancado = 3L, outro = 2L)
  fustes %>%
    group_by(codigo_fazenda, fazenda, parcela, id_ind) %>%
    summarise(
      n_especies_id = n_distinct(especie),
      especie = especie[which.max(tabulate(match(especie, unique(especie))))],
      linha = if (all(is.na(linha))) NA_integer_ else min(linha, na.rm = TRUE),
      n_fustes = n(),
      n_fustes_dap = sum(is.finite(dap_cm) & dap_cm > 0),
      deq_cm = {
        d <- dap_cm[is.finite(dap_cm) & dap_cm > 0]
        if (length(d) == 0) NA_real_ else sqrt(sum(d^2))
      },
      altura_m = if (all(!is.finite(altura_m))) NA_real_ else stats::median(altura_m[is.finite(altura_m)]),
      dc_m = if (all(!is.finite(dc_m))) NA_real_ else stats::median(dc_m[is.finite(dc_m)]),
      # Protocolo: muda é o indivíduo com DAP inferior a 1 cm.
      muda = sum(is.finite(dap_cm) & dap_cm >= LIMIAR_INCLUSAO_CM) == 0,
      regeneracao = if (all(is.na(regeneracao))) NA else any(regeneracao == "sim", na.rm = TRUE),
      fito = {
        g <- gravidade[fito]
        g <- g[!is.na(g)]
        if (length(g) == 0) NA_character_ else names(gravidade)[max(g)]
      },
      .groups = "drop"
    ) %>%
    mutate(
      deq_cm = ifelse(is.finite(deq_cm) & deq_cm > 0, deq_cm, NA_real_),
      altura_m = ifelse(is.finite(altura_m) & altura_m > 0, altura_m, NA_real_)
    )
}

densidade_de_um_sitio <- function(tab) {
  partes <- str_split(str_squish(tab$especie), "\\s+")
  genus <- vapply(partes, function(p) p[1], character(1))
  species <- vapply(partes, function(p) {
    if (length(p) < 2) return(NA_character_)
    paste(p[-1], collapse = " ")
  }, character(1))
  species[ascii_lower(species) %in% c("sp.", "sp", "spp.", "spp")] <- NA_character_
  sem_nome <- ascii_lower(genus) %in% c("ni", "n/i", "outra_especie", "outra")
  genus[sem_nome] <- NA_character_
  species[sem_nome] <- NA_character_
  # A média de fallback é a do próprio sítio. Misturar as duas floras
  # mudaria a biomassa das árvores sem nome na outra fazenda.
  madeira <- getWoodDensity(
    genus = genus,
    species = species,
    stand = rep(tab$codigo_fazenda[1], nrow(tab)),
    region = "SouthAmericaTrop",
    verbose = FALSE
  )
  tab %>%
    mutate(
      genus = genus,
      species = species,
      wd = madeira$meanWD,
      wd_sd = madeira$sdWD,
      wd_nivel = as.character(madeira$levelWD),
      wd_n = madeira$nInd
    ) %>%
    select(codigo_fazenda, especie, genus, species, wd, wd_sd, wd_nivel, wd_n)
}

densidade_madeira <- function(individuos) {
  tab <- individuos %>% distinct(codigo_fazenda, especie)
  bind_rows(lapply(split(tab, tab$codigo_fazenda), densidade_de_um_sitio))
}

biomassa_individuos <- function(individuos, fustes, wd) {
  ind <- individuos %>%
    left_join(wd %>% select(codigo_fazenda, especie, wd_ind = wd, wd_nivel), by = c("codigo_fazenda", "especie"))
  f <- fustes %>%
    left_join(
      ind %>% select(codigo_fazenda, parcela, id_ind, altura_ind = altura_m),
      by = c("codigo_fazenda", "parcela", "id_ind")
    ) %>%
    left_join(wd %>% select(codigo_fazenda, especie, wd, wd_nivel_fuste = wd_nivel), by = c("codigo_fazenda", "especie")) %>%
    mutate(
      altura_uso = ifelse(is.finite(altura_m), altura_m, altura_ind),
      g_m2_fuste = ifelse(
        is.finite(dap_cm) & dap_cm >= LIMIAR_DAP_CM,
        pi * (dap_cm / 200)^2,
        NA_real_
      )
    )
  f$agb_kg <- NA_real_
  ok <- is.finite(f$dap_cm) & is.finite(f$altura_uso) & is.finite(f$wd) & f$dap_cm >= LIMIAR_INCLUSAO_CM
  if (any(ok)) {
    f$agb_kg[ok] <- as.numeric(computeAGB(D = f$dap_cm[ok], WD = f$wd[ok], H = f$altura_uso[ok])) * 1000
  }
  # Diâmetro equivalente, só como sensibilidade: junta fustes pequenos e pode
  # cruzar o limiar de 5 cm sem que nenhum fuste esteja no domínio da equação.
  deq <- ind %>%
    mutate(
      agb5_deq_kg = {
        valor <- rep(NA_real_, n())
        usar <- is.finite(deq_cm) & is.finite(altura_m) & is.finite(wd_ind) & deq_cm >= LIMIAR_DAP_CM
        if (any(usar)) {
          valor[usar] <- as.numeric(computeAGB(D = deq_cm[usar], WD = wd_ind[usar], H = altura_m[usar])) * 1000
        }
        valor
      }
    ) %>%
    select(codigo_fazenda, parcela, id_ind, agb5_deq_kg)

  por_ind <- f %>%
    group_by(codigo_fazenda, parcela, id_ind) %>%
    summarise(
      agb5_kg = if (any(is.finite(agb_kg) & dap_cm >= LIMIAR_DAP_CM)) {
        sum(agb_kg[is.finite(agb_kg) & dap_cm >= LIMIAR_DAP_CM])
      } else {
        NA_real_
      },
      agb1_kg = if (any(is.finite(agb_kg))) sum(agb_kg, na.rm = TRUE) else NA_real_,
      g_m2 = if (any(is.finite(g_m2_fuste))) sum(g_m2_fuste, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  ind %>%
    left_join(por_ind, by = c("codigo_fazenda", "parcela", "id_ind")) %>%
    left_join(deq, by = c("codigo_fazenda", "parcela", "id_ind")) %>%
    mutate(
      wd = wd_ind,
      carbono5_kg = agb5_kg * FRACAO_CARBONO,
      carbono1_kg = agb1_kg * FRACAO_CARBONO
    ) %>%
    select(-wd_ind)
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
  nomes <- names(tab)[tab == max(tab)]
  paste(sort(nomes), collapse = "; ")
}

resumir_parcelas <- function(pontos, poligonos, individuos) {
  base <- pontos %>%
    st_drop_geometry() %>%
    distinct(
      codigo_fazenda, fazenda, talhao, parcela, tratamento, data_coleta,
      solo, declividade, formigas, cupim, n_cupinzeiros,
      perc_folha_larga, perc_folha_estreita, erosao, obs
    )
  area <- poligonos %>%
    st_drop_geometry() %>%
    select(codigo_fazenda, parcela, area_m2, lado_max_m, lado_min_m, geometria_valida, geometria_duplicada)
  if (nrow(individuos) == 0) {
    return(base %>% left_join(area, by = c("codigo_fazenda", "parcela")))
  }
  sintese <- individuos %>%
    group_by(codigo_fazenda, parcela) %>%
    summarise(
      n_individuos = n(),
      n_mudas = sum(muda, na.rm = TRUE),
      n_dap5 = sum(is.finite(agb5_kg) & agb5_kg > 0),
      n_dap1 = sum(is.finite(agb1_kg) & agb1_kg > 0),
      riqueza = n_distinct(especie),
      shannon = shannon_nat(as.numeric(table(especie))),
      especie_pred_n = predominante(nome_exibicao(especie), rep(1, n())),
      especie_pred_agb = predominante(
        nome_exibicao(especie),
        ifelse(is.finite(agb5_kg), agb5_kg, 0)
      ),
      agb5_kg = sum(agb5_kg, na.rm = TRUE),
      agb1_kg = sum(agb1_kg, na.rm = TRUE),
      carbono5_kg = sum(carbono5_kg, na.rm = TRUE),
      area_basal_m2 = sum(g_m2, na.rm = TRUE),
      altura_media_m = mean(altura_m[is.finite(altura_m)], na.rm = TRUE),
      dap_medio_cm = mean(deq_cm[is.finite(deq_cm) & deq_cm >= LIMIAR_INCLUSAO_CM], na.rm = TRUE),
      prop_sadio = mean(fito == "sadio", na.rm = TRUE),
      prop_inicial = mean(fito == "inicial", na.rm = TRUE),
      prop_avancado = mean(fito == "avancado", na.rm = TRUE),
      prop_regen = if (all(is.na(regeneracao))) NA_real_ else mean(regeneracao, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      across(c(altura_media_m, dap_medio_cm), ~ ifelse(is.finite(.x), .x, NA_real_))
    )
  base %>%
    left_join(area, by = c("codigo_fazenda", "parcela")) %>%
    left_join(sintese, by = c("codigo_fazenda", "parcela")) %>%
    mutate(
      agb5_mg_ha = ifelse(geometria_valida, agb5_kg / 1000 / (area_m2 / 10000), NA_real_),
      agb1_mg_ha = ifelse(geometria_valida, agb1_kg / 1000 / (area_m2 / 10000), NA_real_),
      carbono5_mg_ha = agb5_mg_ha * FRACAO_CARBONO,
      densidade_ha = ifelse(geometria_valida, n_individuos / (area_m2 / 10000), NA_real_),
      g_m2_ha = ifelse(geometria_valida, area_basal_m2 / (area_m2 / 10000), NA_real_),
      # Taxa padronizada para a área de projeto (0,10 ha). Não substitui a área medida.
      agb5_mg_ha_010 = agb5_kg / 1000 / 0.10,
      tem_inventario = !is.na(n_individuos)
    )
}

resumo_especies <- function(individuos) {
  individuos %>%
    group_by(especie) %>%
    summarise(
      n = n(),
      agb5_kg = sum(agb5_kg, na.rm = TRUE),
      agb1_kg = sum(agb1_kg, na.rm = TRUE),
      g_m2 = sum(g_m2, na.rm = TRUE),
      wd = wd[1],
      wd_nivel = wd_nivel[1],
      .groups = "drop"
    ) %>%
    mutate(
      especie_nome = nome_exibicao(especie),
      pct_n = n / sum(n),
      pct_agb5 = agb5_kg / sum(agb5_kg)
    ) %>%
    arrange(desc(n))
}

composicao_parcela <- function(individuos) {
  individuos %>%
    group_by(codigo_fazenda, parcela, especie) %>%
    summarise(n = n(), agb5_kg = sum(agb5_kg, na.rm = TRUE), .groups = "drop") %>%
    group_by(codigo_fazenda, parcela) %>%
    mutate(pct_n = n / sum(n)) %>%
    ungroup()
}

teste_kruskal <- function(dados, variavel, grupo = "talhao", min_n = 2) {
  d <- dados %>%
    filter(tem_inventario, !is.na(.data[[variavel]]), !is.na(.data[[grupo]])) %>%
    mutate(grupo = as.character(.data[[grupo]]))
  n_por <- d %>% count(grupo)
  usar <- n_por %>% filter(n >= min_n) %>% pull(grupo)
  d <- d %>% filter(grupo %in% usar)
  if (n_distinct(d$grupo) < 2 || nrow(d) < 5) {
    return(tibble(
      variavel = variavel, teste = "Kruskal-Wallis", estatistica = NA_real_,
      gl = NA_real_, p = NA_real_, efeito = NA_real_, n = nrow(d),
      grupos = paste(sort(unique(d$grupo)), collapse = ", ")
    ))
  }
  kt <- kruskal.test(d[[variavel]] ~ factor(d$grupo))
  n <- nrow(d)
  eps <- unname(kt$statistic) * (n + 1) / (n^2 - 1)
  tibble(
    variavel = variavel,
    teste = "Kruskal-Wallis",
    estatistica = unname(kt$statistic),
    gl = unname(kt$parameter),
    p = kt$p.value,
    efeito = unname(eps),
    n = n,
    grupos = paste(sprintf("%s (n=%d)", n_por$grupo[match(usar, n_por$grupo)], n_por$n[match(usar, n_por$grupo)]), collapse = "; ")
  )
}

teste_wilcoxon_binario <- function(dados, variavel, fator) {
  d <- dados %>%
    filter(tem_inventario, !is.na(.data[[variavel]]), .data[[fator]] %in% c("sim", "nao"))
  n_por <- d %>% count(.data[[fator]])
  if (nrow(n_por) < 2 || any(n_por$n < 3)) {
    return(tibble(variavel = variavel, fator = fator, teste = "Wilcoxon", estatistica = NA_real_, p = NA_real_, n = nrow(d)))
  }
  wt <- wilcox.test(d[[variavel]] ~ factor(d[[fator]]), exact = FALSE)
  tibble(
    variavel = variavel,
    fator = fator,
    teste = "Wilcoxon",
    estatistica = unname(wt$statistic),
    p = wt$p.value,
    n = nrow(d)
  )
}

spearman_teste <- function(dados, x, y) {
  d <- dados %>% filter(tem_inventario, is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(d) < 6) {
    return(tibble(x = x, y = y, rho = NA_real_, p = NA_real_, n = nrow(d)))
  }
  ct <- cor.test(d[[x]], d[[y]], method = "spearman", exact = FALSE)
  tibble(x = x, y = y, rho = unname(ct$estimate), p = ct$p.value, n = nrow(d))
}

pairwise_grupo <- function(dados, variavel, grupo) {
  d <- dados %>%
    filter(tem_inventario, is.finite(.data[[variavel]])) %>%
    mutate(grupo = as.character(.data[[grupo]]))
  n_por <- d %>% count(grupo)
  d <- d %>% filter(grupo %in% n_por$grupo[n_por$n >= 2])
  if (dplyr::n_distinct(d$grupo) < 2) {
    return(tibble(variavel = character(), grupo1 = character(), grupo2 = character(), p_holm = numeric()))
  }
  pw <- pairwise.wilcox.test(
    d[[variavel]],
    factor(d$grupo),
    p.adjust.method = "holm",
    exact = FALSE
  )
  tab <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
  names(tab) <- c("grupo1", "grupo2", "p_holm")
  tab %>%
    filter(!is.na(p_holm)) %>%
    mutate(variavel = variavel, .before = 1)
}

teste_fazendas <- function(dados, variavel) {
  d <- dados %>%
    filter(tem_inventario, is.finite(.data[[variavel]]), !is.na(codigo_fazenda))
  if (dplyr::n_distinct(d$codigo_fazenda) < 2) {
    return(tibble(variavel = variavel, teste = "Wilcoxon", estatistica = NA_real_, p = NA_real_, n = nrow(d)))
  }
  wt <- wilcox.test(d[[variavel]] ~ factor(d$codigo_fazenda), exact = FALSE)
  med <- d %>%
    group_by(codigo_fazenda) %>%
    summarise(mediana = median(.data[[variavel]], na.rm = TRUE), .groups = "drop")
  tibble(
    variavel = variavel,
    teste = "Wilcoxon",
    estatistica = unname(wt$statistic),
    p = wt$p.value,
    n = nrow(d),
    mediana_bh = med$mediana[med$codigo_fazenda == "BH"][1],
    mediana_er = med$mediana[med$codigo_fazenda == "ER"][1]
  )
}

ajustar_holm <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "holm")
  out
}

permanova_grupo <- function(individuos, resumo, coluna, sitio = NULL) {
  # `sitio` não pode se chamar fazenda: a tabela tem essa coluna e o
  # filtro passaria a comparar o código com o nome por extenso.
  mapa <- resumo %>% filter(tem_inventario)
  if (!is.null(sitio)) mapa <- mapa %>% filter(.data$codigo_fazenda == sitio)
  mapa <- mapa %>%
    mutate(
      grupo_f = as.character(.data[[coluna]]),
      id = paste(codigo_fazenda, parcela, sep = "_")
    )
  n_g <- mapa %>% count(grupo_f)
  mapa <- mapa %>% filter(grupo_f %in% n_g$grupo_f[n_g$n >= 2])
  comp <- individuos %>%
    mutate(id = paste(codigo_fazenda, parcela, sep = "_")) %>%
    filter(id %in% mapa$id) %>%
    count(id, especie)
  vazio <- list(r2 = NA_real_, p = NA_real_, n = 0L, n_permutacoes = 999L)
  if (nrow(comp) == 0 || n_distinct(mapa$grupo_f) < 2) return(vazio)
  mat <- comp %>%
    pivot_wider(names_from = especie, values_from = n, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat$id
  mat$id <- NULL
  mat <- as.matrix(mat)
  grupos <- mapa$grupo_f[match(rownames(mat), mapa$id)]
  set.seed(260319)
  fit <- adonis2(mat ~ factor(grupos), method = "hellinger", permutations = 999)
  tibble_fit <- as.data.frame(fit)
  list(
    r2 = tibble_fit$R2[1],
    p = tibble_fit$`Pr(>F)`[1],
    n = nrow(mat),
    n_permutacoes = 999L
  )
}

rotulo_variavel <- function(x) {
  c(
    agb5_kg = "Biomassa aérea, D ≥ 5 cm (kg/parcela)",
    agb1_kg = "Biomassa aérea, D ≥ 1 cm (kg/parcela)",
    n_individuos = "Indivíduos vivos",
    riqueza = "Riqueza de espécies",
    altura_media_m = "Altura média (m)",
    area_basal_m2 = "Área basal, D ≥ 5 cm (m²/parcela)",
    prop_sadio = "Proporção de indivíduos sadios",
    densidade_ha = "Densidade (ind./ha)",
    agb5_mg_ha = "Biomassa, D ≥ 5 cm (Mg/ha)"
  )[x] %>% unname()
}

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

cores_talhao <- function(talhoes) {
  pal <- c("558" = "#0072B2", "559" = "#E69F00", "560" = "#009E73", "561" = "#CC79A7")
  fora <- setdiff(as.character(talhoes), names(pal))
  if (length(fora)) {
    extra <- setNames(hue_pal()(length(fora)), fora)
    pal <- c(pal, extra)
  }
  pal
}

filtrar_fazenda <- function(res, cod) {
  res$pontos <- res$pontos %>% filter(codigo_fazenda == cod)
  res$poligonos <- res$poligonos %>% filter(codigo_fazenda == cod)
  res$individuos <- res$individuos %>% filter(codigo_fazenda == cod)
  res$resumo <- res$resumo %>% filter(codigo_fazenda == cod)
  res$especies <- resumo_especies(res$individuos)
  res$composicao <- composicao_parcela(res$individuos)
  res
}

ordem_tratamento <- function(x) {
  n <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  n[is.na(n)] <- 10^6
  order(n, x)
}

grafico_mapa <- function(res, titulo = "Parcelas com geometria consistente") {
  pol <- res$poligonos %>% filter(geometria_valida)
  if (nrow(pol) == 0) {
    return(
      ggplot() +
        labs(title = titulo, subtitle = "Nenhuma parcela com geometria aceita.") +
        tema_relatorio()
    )
  }
  pts <- st_as_sf(res$pontos, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  pts <- pts %>%
    filter(paste(codigo_fazenda, parcela) %in% paste(pol$codigo_fazenda, pol$parcela))
  # Extensão só das geometrias válidas, para o vértice errado da parcela 5 não achatar o mapa.
  caixa <- st_bbox(pol)
  ggplot() +
    geom_sf(data = pol, aes(fill = agb5_kg), color = "white", linewidth = 0.3) +
    geom_sf(data = pts, size = 0.6, color = "#222222") +
    geom_sf_text(
      data = pol,
      aes(label = parcela),
      size = 2.5, color = "black", fontface = "bold", check_overlap = TRUE
    ) +
    scale_fill_viridis_c(
      option = "cividis",
      name = "Biomassa (kg/parcela)",
      labels = label_number(decimal.mark = ",", big.mark = "."),
      na.value = "grey80"
    ) +
    coord_sf(
      xlim = c(caixa["xmin"] - 0.0004, caixa["xmax"] + 0.0004),
      ylim = c(caixa["ymin"] - 0.0004, caixa["ymax"] + 0.0004),
      expand = FALSE
    ) +
    labs(
      title = titulo,
      subtitle = "Biomassa dos fustes com DAP ≥ 5 cm. Cinza: sem inventário quantitativo.",
      x = NULL, y = NULL
    ) +
    tema_relatorio() +
    theme(legend.position = "right")
}

paleta_grupo <- function(d) {
  um_talhao <- dplyr::n_distinct(d$talhao) <= 1
  if (um_talhao) {
    niveis <- unique(d$tratamento)[ordem_tratamento(unique(d$tratamento))]
    d$grupo <- factor(d$tratamento, levels = niveis)
    pal <- c(
      "T1" = "#0072B2", "T2" = "#E69F00", "T3" = "#009E73", "T4" = "#CC79A7",
      "T5" = "#D55E00", "T6" = "#56B4E9", "T7" = "#F0E442"
    )
    list(d = d, pal = pal[niveis], legenda = "Tratamento", um_talhao = TRUE)
  } else {
    d$grupo <- factor(d$talhao)
    list(d = d, pal = cores_talhao(d$talhao), legenda = "Talhão", um_talhao = FALSE)
  }
}

grafico_biomassa <- function(res) {
  d <- res$resumo %>%
    filter(tem_inventario) %>%
    mutate(rotulo = paste0(tratamento, "\nP", parcela))
  g <- paleta_grupo(d)
  d <- g$d
  subtitulo <- if (g$um_talhao) {
    "Soma dos fustes com DAP ≥ 5 cm. A cor é o tratamento."
  } else {
    "Soma dos fustes com DAP ≥ 5 cm. Nesta fazenda cada tratamento é uma parcela."
  }
  ggplot(d, aes(reorder(rotulo, parcela), agb5_kg, fill = grupo)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = g$pal, name = g$legenda) +
    scale_y_continuous(labels = label_number(decimal.mark = ",", big.mark = "."), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "Biomassa aérea por parcela",
      subtitle = subtitulo,
      x = "Tratamento / parcela",
      y = "Biomassa aérea (kg)"
    ) +
    tema_relatorio() +
    theme(axis.text.x = element_text(size = 7))
}

grafico_individuos <- function(res) {
  d <- res$resumo %>%
    filter(tem_inventario) %>%
    mutate(rotulo = paste0("P", parcela))
  g <- paleta_grupo(d)
  d <- g$d
  ggplot(d, aes(reorder(rotulo, parcela), n_individuos, fill = grupo)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = g$pal, name = g$legenda) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "Indivíduos vivos registrados",
      subtitle = "A planilha não lista mortos. Este número é o estoque vivo, não a taxa de sobrevivência.",
      x = "Parcela",
      y = "Indivíduos"
    ) +
    tema_relatorio() +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5))
}

grafico_talhoes <- function(res) {
  d <- res$resumo %>%
    filter(tem_inventario) %>%
    select(talhao, agb5_kg, n_individuos, riqueza, altura_media_m) %>%
    pivot_longer(-talhao, names_to = "variavel", values_to = "valor") %>%
    mutate(
      talhao_f = factor(talhao),
      variavel = recode(
        variavel,
        agb5_kg = "Biomassa (kg)",
        n_individuos = "Indivíduos vivos",
        riqueza = "Riqueza",
        altura_media_m = "Altura média (m)"
      )
    )
  pal <- cores_talhao(d$talhao)
  ggplot(d, aes(talhao_f, valor, fill = talhao_f)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.08, size = 1.4, alpha = 0.8) +
    facet_wrap(~variavel, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      title = "Parcelas dentro de cada talhão",
      subtitle = "Cada ponto é uma parcela. O teste usa só os talhões com pelo menos duas parcelas.",
      x = "Talhão", y = NULL
    ) +
    tema_relatorio()
}

grafico_especies <- function(res) {
  top <- res$especies %>% slice_head(n = 12) %>% mutate(especie_nome = factor(especie_nome, levels = rev(especie_nome)))
  a <- ggplot(top, aes(n, especie_nome)) +
    geom_col(fill = "#0072B2", width = 0.72) +
    labs(title = "Abundância", x = "Indivíduos", y = NULL) +
    tema_relatorio()
  b <- ggplot(top, aes(agb5_kg, especie_nome)) +
    geom_col(fill = "#E69F00", width = 0.72) +
    scale_x_continuous(labels = label_number(decimal.mark = ",", big.mark = ".")) +
    labs(title = "Biomassa (D ≥ 5 cm)", x = "kg", y = NULL) +
    tema_relatorio()
  a + b +
    plot_annotation(
      title = "Doze espécies mais abundantes",
      subtitle = "A ordem é a da abundância. A biomassa não segue a mesma ordem.",
      theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, color = "#444444"))
    )
}

grafico_composicao <- function(res) {
  top <- res$especies %>% slice_head(n = 6) %>% pull(especie)
  d <- res$composicao %>%
    mutate(grupo = ifelse(especie %in% top, nome_exibicao(especie), "Outras")) %>%
    group_by(codigo_fazenda, parcela, grupo) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(codigo_fazenda, parcela) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()
  ordem_sp <- c(nome_exibicao(top), "Outras")
  ordem_p <- res$resumo %>%
    filter(tem_inventario) %>%
    arrange(talhao, parcela) %>%
    mutate(rotulo = paste0("P", parcela)) %>%
    pull(rotulo)
  d <- d %>%
    mutate(
      grupo = factor(grupo, levels = ordem_sp),
      rotulo = factor(paste0("P", parcela), levels = ordem_p)
    )
  pal <- c(
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#999999"
  )
  names(pal) <- ordem_sp
  ggplot(d, aes(rotulo, pct, fill = grupo)) +
    geom_col(width = 0.86) +
    scale_y_continuous(labels = label_percent(decimal.mark = ",", suffix = "%", accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(
      title = "Composição de cada parcela",
      subtitle = "Seis espécies mais abundantes no conjunto e o restante agrupado.",
      x = "Parcela", y = "Indivíduos"
    ) +
    tema_relatorio() +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5))
}

grafico_sanidade <- function(res) {
  d <- res$individuos %>%
    filter(!is.na(fito)) %>%
    count(codigo_fazenda, parcela, fito) %>%
    mutate(
      fito = recode(fito, sadio = "Sadio", inicial = "Deterioração inicial", avancado = "Deterioração avançada", outro = "Outro"),
      rotulo = paste0("P", parcela)
    )
  ordem <- res$resumo %>% filter(tem_inventario) %>% arrange(parcela) %>% mutate(rotulo = paste0("P", parcela)) %>% pull(rotulo)
  d$rotulo <- factor(d$rotulo, levels = ordem)
  d$fito <- factor(d$fito, levels = c("Sadio", "Deterioração inicial", "Deterioração avançada", "Outro"))
  ggplot(d, aes(rotulo, n, fill = fito)) +
    geom_col(width = 0.86, position = "fill") +
    scale_y_continuous(labels = label_percent(decimal.mark = ",", suffix = "%", accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = c("Sadio" = "#009E73", "Deterioração inicial" = "#E69F00", "Deterioração avançada" = "#D55E00", "Outro" = "#999999"), name = NULL) +
    labs(
      title = "Estado fitossanitário",
      subtitle = "Proporção dos indivíduos registrados. Não é mortalidade: árvore morta não entrou na planilha.",
      x = "Parcela", y = "Indivíduos"
    ) +
    tema_relatorio() +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5))
}

grafico_tratamentos <- function(res) {
  d <- res$resumo %>%
    filter(tem_inventario) %>%
    select(tratamento, agb5_kg, n_individuos, riqueza, altura_media_m, prop_sadio) %>%
    pivot_longer(-tratamento, names_to = "variavel", values_to = "valor") %>%
    mutate(
      tratamento = factor(tratamento, levels = unique(tratamento)[ordem_tratamento(unique(tratamento))]),
      variavel = recode(
        variavel,
        agb5_kg = "Biomassa (kg)",
        n_individuos = "Indivíduos vivos",
        riqueza = "Riqueza",
        altura_media_m = "Altura média (m)",
        prop_sadio = "Proporção sadia"
      )
    )
  ggplot(d, aes(tratamento, valor)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, fill = "#56B4E9", alpha = 0.85) +
    geom_jitter(width = 0.08, size = 1.5, alpha = 0.85) +
    facet_wrap(~variavel, scales = "free_y", nrow = 1) +
    labs(
      title = "Quatro parcelas em cada tratamento",
      subtitle = "Cada ponto é uma parcela de Entre Rios. T1 a T7 são o fator replicado.",
      x = "Tratamento", y = NULL
    ) +
    tema_relatorio()
}

grafico_fazendas <- function(res) {
  d <- res$resumo %>%
    filter(tem_inventario) %>%
    transmute(
      fazenda = recode(codigo_fazenda, BH = "Belo Horizonte", ER = "Entre Rios"),
      agb5_kg, n_individuos, riqueza, altura_media_m
    ) %>%
    pivot_longer(-fazenda, names_to = "variavel", values_to = "valor") %>%
    mutate(
      variavel = recode(
        variavel,
        agb5_kg = "Biomassa (kg)",
        n_individuos = "Indivíduos vivos",
        riqueza = "Riqueza",
        altura_media_m = "Altura média (m)"
      )
    )
  ggplot(d, aes(fazenda, valor, fill = fazenda)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.55) +
    facet_wrap(~variavel, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Belo Horizonte" = "#0072B2", "Entre Rios" = "#E69F00"), guide = "none") +
    labs(
      title = "As duas fazendas, lado a lado",
      subtitle = "Cada ponto é uma parcela inventariada. As floras e as regiões não são as mesmas.",
      x = NULL, y = NULL
    ) +
    tema_relatorio()
}

grafico_estrutura <- function(res) {
  d <- res$individuos %>% filter(is.finite(deq_cm), deq_cm >= LIMIAR_INCLUSAO_CM, is.finite(altura_m))
  ggplot(d, aes(deq_cm, altura_m)) +
    geom_point(alpha = 0.25, size = 1, color = "#0072B2") +
    geom_vline(xintercept = LIMIAR_DAP_CM, linetype = "dashed", color = "#D55E00") +
    scale_x_continuous(labels = label_number(decimal.mark = ",")) +
    scale_y_continuous(labels = label_number(decimal.mark = ",")) +
    labs(
      title = "Altura e diâmetro equivalente",
      subtitle = "Cada ponto é um indivíduo. A linha marca 5 cm, limite usado na biomassa principal.",
      x = "Diâmetro equivalente (cm)", y = "Altura (m)"
    ) +
    tema_relatorio()
}

grafico_covariavel <- function(res) {
  d <- res$resumo %>% filter(tem_inventario, is.finite(perc_folha_larga), is.finite(agb5_kg))
  ggplot(d, aes(perc_folha_larga, agb5_kg)) +
    geom_point(aes(color = factor(talhao)), size = 2.4) +
    scale_x_continuous(labels = label_percent(decimal.mark = ",", suffix = "%", accuracy = 1)) +
    scale_y_continuous(labels = label_number(decimal.mark = ",", big.mark = ".")) +
    scale_color_manual(values = cores_talhao(d$talhao), name = "Talhão") +
    labs(
      title = "Folha larga e biomassa da parcela",
      subtitle = "Classes de cobertura anotadas na ficha qualitativa.",
      x = "Cobertura de folha larga", y = "Biomassa, D ≥ 5 cm (kg)"
    ) +
    tema_relatorio()
}

escrever_gpkg <- function(pontos, poligonos, pasta) {
  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
  campos_ponto <- c(
    "codigo_fazenda", "fazenda", "talhao", "parcela", "tratamento", "ponto",
    "lat", "lon", "foto", "data_coleta", "protocolo", "solo", "declividade",
    "formigas", "cupim", "erosao", "perc_folha_larga", "perc_folha_estreita"
  )
  campos_pol <- c(
    "codigo_fazenda", "fazenda", "talhao", "parcela", "tratamento", "data_coleta",
    "area_m2", "lado_max_m", "lado_min_m", "geometria_valida", "geometria_duplicada",
    "agb5_kg", "agb1_kg", "n_individuos", "riqueza", "solo", "declividade",
    "formigas", "cupim", "perc_folha_larga", "perc_folha_estreita", "obs"
  )
  pts <- st_as_sf(pontos, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  arquivos <- character()
  for (cod in sort(unique(pts$codigo_fazenda))) {
    arq <- file.path(pasta, sprintf("parcelas_%s.gpkg", cod))
    if (file.exists(arq)) unlink(arq)
    p1 <- pts %>% filter(codigo_fazenda == cod)
    p1 <- p1[, intersect(campos_ponto, names(p1))]
    g1 <- poligonos %>% filter(codigo_fazenda == cod)
    g1 <- g1[, intersect(c(campos_pol, attr(g1, "sf_column")), names(g1))]
    st_write(p1, arq, layer = "pontos", quiet = TRUE, delete_dsn = TRUE)
    st_write(g1, arq, layer = "parcelas", quiet = TRUE, append = TRUE)
    arquivos <- c(arquivos, arq)
  }
  arquivos
}

pasta_entrada <- function(raiz) {
  candidatos <- c(file.path(raiz, "0-raw"), file.path(raiz, "dados"))
  hit <- candidatos[dir.exists(candidatos)]
  if (length(hit) == 0) stop("Não há pasta 0-raw nem dados em ", raiz)
  hit[1]
}

run_analise <- function(raiz) {
  dados <- pasta_entrada(raiz)
  saida <- file.path(raiz, "saida")
  dir.create(saida, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(saida, "figuras"), recursive = TRUE, showWarnings = FALSE)

  bruto <- ler_planilhas(dados)
  pontos <- preparar_qualitativa(bruto$quali)
  geom <- montar_geometrias(pontos)
  fustes <- preparar_fustes(bruto$quant, pontos)
  individuos <- agregar_individuos(fustes)
  wd <- densidade_madeira(individuos)
  individuos <- biomassa_individuos(individuos, fustes, wd)
  agb_fuste <- sum(individuos$agb5_deq_kg, na.rm = TRUE)
  resumo <- resumir_parcelas(pontos, geom$poligonos, individuos)
  # Colar métricas no polígono para o GeoPackage.
  metricas <- resumo %>%
    select(codigo_fazenda, parcela, agb5_kg, agb1_kg, n_individuos, riqueza)
  poligonos <- geom$poligonos %>%
    left_join(metricas, by = c("codigo_fazenda", "parcela"))
  especies <- resumo_especies(individuos)
  composicao <- composicao_parcela(individuos)

  variaveis_teste <- c(
    "agb5_kg", "agb1_kg", "n_individuos", "riqueza",
    "altura_media_m", "area_basal_m2", "prop_sadio"
  )
  pares_de <- function(tab, dados, grupo) {
    alvos <- tab$variavel[is.finite(tab$p) & tab$p < 0.05]
    if (length(alvos) == 0) {
      return(tibble(variavel = character(), grupo1 = character(), grupo2 = character(), p_holm = numeric()))
    }
    bind_rows(lapply(alvos, function(v) pairwise_grupo(dados, v, grupo)))
  }
  bh <- resumo %>% filter(codigo_fazenda == "BH")
  er <- resumo %>% filter(codigo_fazenda == "ER")
  kw <- bind_rows(lapply(variaveis_teste, function(v) teste_kruskal(bh, v, "talhao")))
  kw$p_holm <- ajustar_holm(kw$p)
  kw$variavel_nome <- rotulo_variavel(kw$variavel)
  kw$fazenda <- "BH"
  pares <- pares_de(kw, bh, "talhao")

  kw_trat <- bind_rows(lapply(variaveis_teste, function(v) teste_kruskal(er, v, "tratamento")))
  kw_trat$p_holm <- ajustar_holm(kw_trat$p)
  kw_trat$variavel_nome <- rotulo_variavel(kw_trat$variavel)
  kw_trat$fazenda <- "ER"
  pares_trat <- pares_de(kw_trat, er, "tratamento")

  testes_fazenda <- bind_rows(lapply(variaveis_teste, function(v) teste_fazendas(resumo, v)))
  testes_fazenda$p_holm <- ajustar_holm(testes_fazenda$p)
  testes_fazenda$variavel_nome <- rotulo_variavel(testes_fazenda$variavel)

  exploratorios <- bind_rows(lapply(c("BH", "ER"), function(faz) {
    sub <- resumo %>% filter(codigo_fazenda == faz)
    out <- bind_rows(
      teste_wilcoxon_binario(sub, "agb5_kg", "formigas"),
      teste_wilcoxon_binario(sub, "n_individuos", "formigas"),
      teste_wilcoxon_binario(sub, "agb5_kg", "cupim"),
      teste_wilcoxon_binario(sub, "n_individuos", "cupim")
    )
    out$p_holm <- ajustar_holm(out$p)
    out$codigo_fazenda <- faz
    out
  }))
  spears <- bind_rows(lapply(c("BH", "ER"), function(faz) {
    sub <- resumo %>% filter(codigo_fazenda == faz)
    out <- bind_rows(
      spearman_teste(sub, "perc_folha_larga", "agb5_kg"),
      spearman_teste(sub, "perc_folha_larga", "n_individuos"),
      spearman_teste(sub, "perc_folha_estreita", "agb5_kg"),
      spearman_teste(sub, "perc_folha_estreita", "n_individuos"),
      spearman_teste(sub, "n_cupinzeiros", "agb5_kg")
    )
    out$p_holm <- ajustar_holm(out$p)
    out$codigo_fazenda <- faz
    out
  }))

  perm <- permanova_grupo(individuos, resumo, "talhao", "BH")
  perm_trat <- permanova_grupo(individuos, resumo, "tratamento", "ER")
  perm_fazenda <- permanova_grupo(individuos, resumo, "codigo_fazenda")

  med_talhao <- resumo %>%
    filter(tem_inventario) %>%
    group_by(codigo_fazenda, talhao) %>%
    summarise(
      n_parcelas = n(),
      mediana_agb5_kg = median(agb5_kg, na.rm = TRUE),
      mediana_agb1_kg = median(agb1_kg, na.rm = TRUE),
      mediana_individuos = median(n_individuos, na.rm = TRUE),
      mediana_riqueza = median(riqueza, na.rm = TRUE),
      mediana_altura_m = median(altura_media_m, na.rm = TRUE),
      mediana_prop_sadio = median(prop_sadio, na.rm = TRUE),
      .groups = "drop"
    )
  med_trat <- resumo %>%
    filter(tem_inventario) %>%
    group_by(codigo_fazenda, tratamento) %>%
    summarise(
      n_parcelas = n(),
      mediana_agb5_kg = median(agb5_kg, na.rm = TRUE),
      mediana_individuos = median(n_individuos, na.rm = TRUE),
      mediana_riqueza = median(riqueza, na.rm = TRUE),
      mediana_altura_m = median(altura_media_m, na.rm = TRUE),
      mediana_prop_sadio = median(prop_sadio, na.rm = TRUE),
      mediana_regen = median(prop_regen, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_parcelas >= 2)

  problemas <- resumo %>%
    mutate(problema = case_when(
      !geometria_valida & geometria_duplicada ~ "geometria inválida e vértices repetidos em outra parcela",
      !geometria_valida ~ "geometria inválida (vértice fora do retângulo da parcela ou área implausível)",
      geometria_duplicada ~ "os quatro vértices coincidem com os de outra parcela",
      !tem_inventario ~ "sem planilha quantitativa",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(problema)) %>%
    select(codigo_fazenda, talhao, parcela, tratamento, area_m2, problema)

  if (any(resumo$tem_inventario) && all(is.na(resumo$especie_pred_agb[resumo$tem_inventario]))) {
    stop("A espécie predominante por biomassa ficou vazia. O cálculo por parcela não usou o vetor de biomassa.")
  }
  gpkg <- escrever_gpkg(pontos, poligonos, saida)

  readr::write_csv(resumo, file.path(saida, "resumo_parcelas.csv"), na = "")
  readr::write_csv(
    individuos %>% select(
      codigo_fazenda, fazenda, parcela, id_ind, linha, especie, n_fustes,
      deq_cm, altura_m, muda, fito, regeneracao, wd, wd_nivel, agb5_kg, agb1_kg, carbono5_kg, g_m2
    ),
    file.path(saida, "individuos.csv"),
    na = ""
  )
  readr::write_csv(especies, file.path(saida, "especies.csv"), na = "")
  readr::write_csv(wd, file.path(saida, "densidade_madeira.csv"), na = "")
  readr::write_csv(kw, file.path(saida, "testes_talhao.csv"), na = "")
  readr::write_csv(pares, file.path(saida, "testes_pares_talhao.csv"), na = "")
  readr::write_csv(kw_trat, file.path(saida, "testes_tratamento.csv"), na = "")
  readr::write_csv(pares_trat, file.path(saida, "testes_pares_tratamento.csv"), na = "")
  readr::write_csv(testes_fazenda, file.path(saida, "testes_fazenda.csv"), na = "")
  readr::write_csv(exploratorios, file.path(saida, "testes_covariaveis.csv"), na = "")
  readr::write_csv(spears, file.path(saida, "correlacoes.csv"), na = "")
  readr::write_csv(med_trat, file.path(saida, "medianas_tratamento.csv"), na = "")

  comparacao_deq <- tibble(
    agb5_fuste_kg = sum(individuos$agb5_kg, na.rm = TRUE),
    agb5_deq_kg = agb_fuste
  )

  res <- list(
    arquivos = bruto$arquivos,
    pontos = pontos,
    poligonos = poligonos,
    fustes = fustes,
    individuos = individuos,
    wd = wd,
    resumo = resumo,
    especies = especies,
    composicao = composicao,
    kw = kw,
    pares = pares,
    kw_trat = kw_trat,
    pares_trat = pares_trat,
    testes_fazenda = testes_fazenda,
    exploratorios = exploratorios,
    spears = spears,
    permanova = perm,
    permanova_tratamento = perm_trat,
    permanova_fazenda = perm_fazenda,
    med_talhao = med_talhao,
    med_trat = med_trat,
    problemas = problemas,
    gpkg = gpkg,
    comparacao_deq = comparacao_deq,
    n_id_conflito = sum(individuos$n_especies_id > 1)
  )
  saveRDS(res, file.path(saida, "resultados.rds"))
  res
}

if (sys.nframe() == 0) {
  raiz <- Sys.getenv("MAXIC_RAIZ", unset = normalizePath(file.path(getwd(), "..")))
  if (!dir.exists(file.path(raiz, "dados"))) raiz <- getwd()
  res <- run_analise(raiz)
  cat("Arquivos:", paste(res$arquivos, collapse = ", "), "\n")
  cat("Parcelas qualitativas:", nrow(distinct(res$pontos, codigo_fazenda, parcela)), "\n")
  cat("Pontos:", nrow(res$pontos), "\n")
  cat("Indivíduos:", nrow(res$individuos), "\n")
  por_fazenda <- res$resumo %>%
    group_by(codigo_fazenda) %>%
    summarise(
      parcelas = sum(tem_inventario),
      agb5 = round(sum(agb5_kg, na.rm = TRUE), 1),
      validas = sum(geometria_valida, na.rm = TRUE),
      .groups = "drop"
    )
  print(por_fazenda)
  cat("AGB D>=5 kg:", round(sum(res$individuos$agb5_kg, na.rm = TRUE), 1), "\n")
  cat("AGB D>=1 kg:", round(sum(res$individuos$agb1_kg, na.rm = TRUE), 1), "\n")
  cat("AGB por fuste kg:", round(res$comparacao_deq$agb5_fuste_kg, 1), "\n")
  cat("GeoPackage:", paste(res$gpkg, collapse = ", "), "\n")
  print(res$kw)
  print(res$problemas)
}
