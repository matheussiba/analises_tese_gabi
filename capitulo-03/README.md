# Inventário MaxiC — capítulo 3

Inventário e biomassa das parcelas do experimento MaxiC nas fazendas Belo Horizonte (BH) e Entre Rios (ER). Cada fazenda tem uma apresentação em PowerPoint própria, com os mesmos slides e os mesmos gráficos.

## Como rodar

É preciso ter R e [Quarto](https://quarto.org/).

```bash
cd capitulo-03
quarto render relatorio/apresentacao-fazenda-BH.qmd
quarto render relatorio/apresentacao-fazenda-ER.qmd
```

Os PowerPoints saem em `capitulo-03/relatorio/`. As tabelas e os GeoPackages saem em `capitulo-03/saida/BH/` e `capitulo-03/saida/ER/`. Para gerar só as tabelas: `Rscript R/analise_fazendas.R`.

Pacotes de R: `dplyr`, `tidyr`, `stringr`, `readxl`, `sf`, `ggplot2`, `patchwork`, `scales`, `BIOMASS`, `vegan`, `readr`, `lme4`.

## Fontes de dados

- `0-raw/parcelas_BH.gpkg` e `0-raw/parcelas_ER.gpkg`: **única fonte de coordenada e de tratamento**. Um polígono por parcela, com o campo de parcela e o de tratamento.
- `0-raw/planilha_coleta_quali-quantitativa_fazenda_BH.xlsx` e `..._ER.xlsx`: indivíduos, fustes, altura, D30, espécie, estado e data. As colunas de coordenada e de tratamento das planilhas não são usadas.

## O que o cálculo faz

- Cada parcela tem 30 × 30 m, ou 900 m², e 150 mudas plantadas. O Mg/ha usa 900 m² fixos. A área desenhada no GeoPackage varia por causa do GPS e só é mostrada como referência.
- Biomassa por indivíduo pela equação de Ferez et al. (2015): ln(AGBW) = 6,039 + 0,945 ln(SA) + 0,961 ln(Ht) + 1,022 ln(ρ).
  - SA é a soma das áreas seccionais dos fustes a 30 cm (D30), em **m²**.
  - Ht é a mediana das alturas do indivíduo, em m.
  - ρ é a densidade da madeira do pacote BIOMASS, em **g cm⁻³**.
  - Com SA em cm² e ρ em kg m⁻³ a equação dá valores absurdos. As unidades acima são as da calibração original.
- D30 ausente ou maior que 50 cm fica sem biomassa e é listado em `saida/<fazenda>/d30_suspeito.csv`.
- A incerteza mostrada nas barras vem só da densidade da madeira (Monte Carlo, 1000 simulações, intervalo de 95%).
- Sobrevivência é indivíduos vivos de plantio ÷ 150. A planilha de ER marca regeneração, que fica fora das comparações entre tratamentos. A de BH não tem essa coluna.
- Comparação entre tratamentos T1 a T7: Kruskal–Wallis com correção de Holm, reamostragem da posição de cada tratamento, robustez do primeiro lugar, Cliff's delta, modelos, diferença detectável e PERMANOVA de Hellinger. Em BH há ainda a análise por talhão.

## Arquivos

| Caminho | Conteúdo |
| --- | --- |
| `0-raw/` | GeoPackages das parcelas, planilhas de campo e croquis |
| `R/analise_fazendas.R` | Leitura, biomassa, testes e gráficos (usado pelos dois PowerPoints) |
| `relatorio/apresentacao-fazenda-BH.qmd` | Apresentação de Belo Horizonte |
| `relatorio/apresentacao-fazenda-ER.qmd` | Apresentação de Entre Rios |
| `relatorio/_slides_*.qmd` | Slides compartilhados pelas duas apresentações |
| `saida/<fazenda>/resumo_parcelas.csv` | Uma linha por parcela |
| `saida/<fazenda>/parcelas_<fazenda>.gpkg` | Polígonos com a biomassa de cada parcela |
| `saida/<fazenda>/tratamentos_*.csv` | Posição, reamostragem, modelos e diferença detectável |

## Código antigo

`R/analise.R` e `relatorio/inventario-maxic-bh.qmd` (PDF) continuam no fluxo anterior, com Chave et al. (2014) e as planilhas antigas. Não foram atualizados para o Ferez nem para os GeoPackages.
