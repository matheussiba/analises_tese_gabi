# Inventário MaxiC — capítulo 2

Leitura das planilhas quali-quantitativas do experimento MaxiC nas fazendas Belo Horizonte e Entre Rios. O script gera um GeoPackage com os quatro vértices de cada parcela e um PDF, feito no Quarto, com biomassa, composição e testes por parcela, tratamento e talhão.

## Como rodar

É preciso ter R e [Quarto](https://quarto.org/).

```bash
cd capitulo-02
quarto render relatorio/inventario-maxic-bh.qmd
quarto render relatorio/apresentacao-maxic.qmd
```

O PDF e o PowerPoint saem em `capitulo-02/relatorio/`. As tabelas e os GeoPackages saem em `capitulo-02/saida/`.

As planilhas ficam em `capitulo-02/0-raw/` ou, se essa pasta não existir, em `capitulo-02/dados/`.

Pacotes de R: `dplyr`, `tidyr`, `stringr`, `readxl`, `sf`, `ggplot2`, `patchwork`, `scales`, `BIOMASS`, `vegan`, `readr`.

## O que o cálculo faz

- A aba `planilha_qualitativa` vira a camada `pontos` (quatro vértices por parcela) e a camada `parcelas`, em WGS 84. A área é calculada em SIRGAS 2000 / UTM, na zona do sítio (24S em Belo Horizonte, 23S em Entre Rios).
- A biomassa aérea soma os fustes com DAP ≥ 5 cm pela equação 4 de Chave et al. (2014). A densidade da madeira vem do pacote BIOMASS (base de Zanne et al., América do Sul tropical).
- Em Belo Horizonte cada tratamento é uma parcela. O teste com repetição é o Kruskal–Wallis entre talhões. Em Entre Rios os tratamentos T1 a T7 têm quatro parcelas cada, e o teste é entre tratamentos. A comparação entre fazendas é um Wilcoxon. A composição usa PERMANOVA com distância de Hellinger. Em Belo Horizonte, a PERMANOVA principal deixa de fora as parcelas que compartilham o mesmo retângulo.
- Em Entre Rios, a biomassa do plantio e a da regeneração saem separadas. O total que inclui a regeneração continua calculado. O intervalo da parcela propaga só o `wd_sd`. A fração de carbono principal continua 0,47.
- As planilhas não registram árvore morta nem o total plantado. O relatório não inventa taxa de sobrevivência: mostra indivíduos vivos e estado fitossanitário. A taxa só entra quando houver um arquivo `mudas_plantadas.csv` com o denominador por parcela.

## Arquivos

| Caminho | Conteúdo |
| --- | --- |
| `capitulo-02/0-raw/` ou `capitulo-02/dados/` | Planilhas de campo |
| `capitulo-02/R/analise.R` | Leitura, geometria, biomassa e testes |
| `capitulo-02/relatorio/inventario-maxic-bh.qmd` | Relatório em PDF |
| `capitulo-02/relatorio/apresentacao-maxic.qmd` | Apresentação em PowerPoint |
| `capitulo-02/saida/parcelas_BH.gpkg` | Pontos e polígonos de Belo Horizonte |
| `capitulo-02/saida/parcelas_ER.gpkg` | Pontos e polígonos de Entre Rios |
| `capitulo-02/saida/resumo_parcelas.csv` | Uma linha por parcela, com biomassa, intervalo da densidade da madeira e, em Entre Rios, plantio e regeneração |
| `capitulo-02/saida/pendencias_campo.csv` | GPS pendente: retângulo compartilhado, anel colapsado, vértice discrepante |
| `capitulo-02/saida/permanova.csv` | PERMANOVA principal e a sensibilidade que trata retângulos repetidos como independentes |
