#### Figures - KARS ID project
#### Shadi Shahatit - RA, JUST 2025
# Libraries ---------------------------------------------------------------

## Note that file paths MUST change!

library(readxl)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(viridis)
library(stringr)
library(zoo)
library(patchwork)
library(fmsb)
library(scales)
library(plotly)

# ClinVar -----------------------------------------------

## Load dataset
clinvar_data <- read_excel("C:/Users/Shadi Shahatit/OneDrive/Desktop/ShinyClinVar_KARS1_2025-11-22.xlsx")

## Rename key columns
clinvar_data <- clinvar_data %>%
  dplyr::rename(
    clinical_significance = ClinicalSignificance,
    protein_position     = pos_aa,
    cadd_phred            = CADD_phred,
    variant_name          = Name
  )

clinvar_data$clinical_significance <- as.factor(clinvar_data$clinical_significance)
clinvar_data$protein_position <- as.numeric(clinvar_data$protein_position)

## Merge clinical significance into 3 groups
clinvar_data$clin_merge <- case_when(
  str_detect(clinvar_data$clinical_significance, "Uncertain|conflicting") ~ "Uncertain/conflicting",
  str_detect(clinvar_data$clinical_significance, "Benign") ~ "Benign group",
  str_detect(clinvar_data$clinical_significance, "Pathogenic") ~ "Pathogenic group",
  TRUE ~ "Uncertain/conflicting"
)

clinvar_data$clin_merge <- factor(
  clinvar_data$clin_merge,
  levels = c("Uncertain/conflicting", "Benign group", "Pathogenic group")
)

## Clean and convert CADD column to numeric
clinvar_data$cadd_phred <- suppressWarnings(as.numeric(str_replace(clinvar_data$cadd_phred, "[A-Za-z]", "")))

## Lollipop plot replacing your scatter plot
ggplot(clinvar_data, aes(x = protein_position, y = cadd_phred, color = clin_merge)) +
  geom_segment(aes(x = protein_position, xend = protein_position, y = 0, yend = cadd_phred), linewidth = 0.8) +
  geom_point(size = 4) +
  theme_minimal() +
  scale_color_viridis(discrete = T) + 
  labs(
    x = "Protein amino-acid position",
    y = "CADD phred score",
    color = "Clinical group"
  )

# Gene-level analysis - Figure 2 info (Domain, MTR, phyloP, AlphaMissense) ---------------------------------------------------------------------

## domain

# for uniport
# 126 206
# 222 575

domains <- tibble::tibble(
  domain = c("IDR","Nucleic acid-binding",".","Lysyl tRNA synthetase domain"
             # "Aminoacyl-tRNA syn. II"
  ),
  start  = c(34,99,155,246
             # 272
  ),
  end    = c(99,242,234,602
             # 603
  )
) %>% mutate(mid = (start + end)/2)

## MTR

mtr_tbl  <- read.delim("C:/Users/Shadi Shahatit/OneDrive/Desktop/MTR_data_ENST00000319410.csv",
                       sep = "\t", header = TRUE, quote = "", check.names = FALSE)

mtr_full <- read.delim("C:/Users/Shadi Shahatit/OneDrive/Desktop/MTR_data_ENST00000319410_v2.csv",
                       sep = "\t", header = TRUE, quote = "", check.names = FALSE)

mtr_window <- mtr_full %>%
  transmute(
    protein_position = as.numeric(protein_position),
    mtr = suppressWarnings(as.numeric(mtr))
  ) %>%
  arrange(protein_position) %>%
  group_by(protein_position) %>%
  summarise(mtr = mean(mtr, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(mtr))

w <- 31L

mtr_fw <- mtr_window %>% arrange(protein_position) %>%
  mutate(mtr_fw31 = rollapply(mtr, width = w, FUN = mean, align = "left", fill = NA, na.rm = TRUE))

## phyloP conservation

lines <- readLines("C:/Users/Shadi Shahatit/OneDrive/Desktop/KARS_exons_100way_phyloP")
phyloP_KARSexons <- read.table(text = grep("^[0-9]", lines, value = TRUE), col.names = c("pos","phyloP"))
phyloP_KARSexons_percodon <- phyloP_KARSexons %>%
  arrange(pos) %>%
  mutate(idx = row_number(), aa = ((idx - 1L) %/% 3L) + 1L) %>%
  group_by(aa) %>%
  summarise(codon_start = min(pos), codon_end = max(pos), phyloP_mean = mean(phyloP, na.rm = TRUE), .groups = "drop")

phy_fw <- phyloP_KARSexons_percodon %>% arrange(aa) %>%
  mutate(phyloP_fw31 = rollapply(phyloP_mean, width = w, FUN = mean, align = "left", fill = NA, na.rm = TRUE))

## Alphafold

alphaMissense_file <- "C:/Users/Shadi Shahatit/OneDrive/Desktop/AlphaMissense-Hotspot-Q15046/AlphaMissense-Hotspot-Q15046.tsv"

alphaMissense_raw <- read.table(
  alphaMissense_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(alphaMissense_raw) <- make.names(names(alphaMissense_raw))

alphaMissense_fw <- alphaMissense_raw %>%
  transmute(
    aa_position = as.integer(position),
    mean_AM = as.numeric(`mean.AM.score`)
  ) %>%
  group_by(aa_position) %>%
  summarise(mean_AM = mean(mean_AM, na.rm = TRUE), .groups = "drop") %>%
  arrange(aa_position) %>%
  mutate(AM_fw31 = rollapply(mean_AM, w, mean, align = "left", fill = NA, na.rm = TRUE))

## plot

protein_length <- max(
  alphaMissense_fw$aa_position,
  mtr_fw$protein_position,
  phy_fw$aa,
  na.rm = TRUE)
xb <- seq(0, protein_length, 100)
xb[1] <- 1L

variant_pos <- 330L

# thresholds
t_blue  <- 1.00
t_gray  <- quantile(mtr_window$mtr, probs = c(0.50), na.rm = TRUE)
t_orng  <- quantile(mtr_window$mtr, probs = c(0.25), na.rm = TRUE)
t_green <- quantile(mtr_window$mtr, probs = c(0.05), na.rm = TRUE)

p_domain <- ggplot() +
  geom_rect(aes(xmin = 1, xmax = protein_length, ymin = 0.0625, ymax = 0.1875),
            fill = "grey85", color = "grey40", linewidth = 0.4) +
  geom_rect(data = domains,
            aes(xmin = start, xmax = end, ymin = 0, ymax = 0.25, fill = domain),
            color = "black", linewidth = 0.4) +
  geom_text(data = domains, aes(x = mid, y = 0.125, label = domain),
            color = "white", size = 3, fontface = "bold") +
  geom_vline(xintercept = 330, color = "#DAA520", linewidth = 1) +
  # geom_vline(xintercept = 330, color = "#DAA520", linewidth = 1.2) +
  scale_fill_manual(values = c(
    "IDR"                     = "#7DAA92",
    "Nucleic acid-binding"    = "#f39c12",
    "."                        = "#8e44ad",
    "Lysyl tRNA synthetase domain"           = "#3A7CA5"
    # "Aminoacyl-tRNA syn. II"  = "#8E6BAF"
  ))+
  # scale_x_continuous(name = "Protein position", limits = c(1, protein_length), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, protein_length, 100))+
  scale_y_continuous(NULL, breaks = NULL, limits = c(0, 0.3), expand = c(0, 0)) +
  labs(x="Amino acid position")+
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "none",
        plot.margin = margin(10, 25, 10, 25))

p_domain <- p_domain +
  scale_x_continuous(limits = c(1, protein_length), breaks = xb, expand = c(0,0))

p_mtr <- ggplot(mtr_fw, aes(x = protein_position, y = mtr_fw31)) +
  geom_hline(yintercept = t_blue,  linetype = "dashed", color = "#377bd1") +
  geom_hline(yintercept = t_gray,  linetype = "dashed", color = "grey45") +
  geom_hline(yintercept = t_orng,  linetype = "dashed", color = "#f39c12") +
  geom_hline(yintercept = t_green, linetype = "dashed", color = "#2ca02c") +
  geom_line(linewidth = 0.9, color = "#d62728", na.rm = TRUE) +
  geom_point(aes(y = mtr), size = 0.6, alpha = 0.25) +
  geom_vline(xintercept = 330, color = "#DAA520", linewidth = 1) +
  labs(x = "", y = "MTR") +
  scale_x_continuous(limits = c(1, protein_length), breaks = xb, expand = c(0,0)) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())

p_phyloP <- ggplot(phy_fw, aes(x = aa, y = phyloP_fw31)) +
  geom_line(linewidth = 0.9, color = "#2c7bb6", na.rm = TRUE) +
  geom_point(aes(y = phyloP_mean), size = 0.5, alpha = 0.25) +
  geom_vline(xintercept = 330, color = "#DAA520", linewidth = 1) +
  labs(x = "", y = "phyloP") +
  scale_x_continuous(limits = c(1, protein_length), breaks = xb, expand = c(0,0)) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())

p_alphafold <- ggplot(alphaMissense_fw, aes(aa_position, AM_fw31)) +
  geom_line(linewidth = 0.9, color = "#6b6ecf", na.rm = TRUE) +
  geom_point(aes(y = mean_AM), size = 0.5, alpha = 0.25) +
  geom_vline(xintercept = variant_pos, color = "#DAA520", linewidth = 1) +
  labs(x = "", y = "AlphaMissense") +
  scale_x_continuous(limits = c(1, protein_length), breaks = xb, expand = c(0, 0)) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

p_clinvar <- ggplot(subset(clinvar_data), aes(x = protein_position, y = cadd_phred, color = clin_merge)) +
  geom_segment(aes(x = protein_position, xend = protein_position, y = 0, yend = cadd_phred), linewidth = 0.8) +
  geom_point(size = 4) +
  geom_vline(xintercept = 330, color = "#DAA520", linewidth = 1) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "", y = "CADD phred score", color = "Clinical group") +
  scale_x_continuous(limits = c(1, protein_length), breaks = xb, expand = c(0,0)) +
  theme_minimal()

## final plot

p_mtr / p_phyloP / p_alphafold / p_domain # p_clinvar

# Functional analysis - Pathogenicity prediction  -------------------------------------------------------------

## OLD

## PredictSNP

PredictSNP <- data.frame(
  PredictSNP1_2 = 87,
  MAPP = 88,
  PhD_SNP = 88,
  PolyPhen1 = 74,
  PolyPhen2 = 81,
  SIFT = 79,
  SNAP = 85,
  PANTHER = 87,
  # PredictSNP2 = 87,
  CADD       = 84,
  DANN       = 77,
  FATHMM     = 73,
  FunSeq2    = 62
)

# scale from 50 to 100
PredictSNP_rad <- rbind(rep(100, ncol(PredictSNP)), rep(0, ncol(PredictSNP)), PredictSNP)

radarchart(PredictSNP_rad,
           axistype = 1,
           seg = 10, 
           pcol = "firebrick",
           pfcol = scales::alpha("red", 0.4),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           caxislabels = seq(0, 100, 10),  # includes 100
           vlcex = 0.8)

PredictSNP_Neu <- data.frame(
  PredictSNP1_2 = 50, MAPP = 50, PhD_SNP = 50, PolyPhen1 = 50, PolyPhen2 = 50,
  SIFT = 50, SNAP = 50, PANTHER = 50,
  CADD = 50, DANN = 50, FATHMM = 50, FunSeq2 = 50)

PredictSNP_Neu_rad <- rbind(rep(100, ncol(PredictSNP_Neu)),  # max
                            rep(0,   ncol(PredictSNP_Neu)),  # min
                            PredictSNP_Neu)                  # 50% ring

radarchart(
  PredictSNP_Neu_rad,
  axistype = 1,
  seg = 10, 
  caxislabels = seq(0, 100, 10),
  pcol = NA,
  pfcol = alpha("grey50", 0.35),
  plwd = 1,
  cglcol = "grey60",
  axislabcol = "black",
  cglty = 1,
  vlcex = 0.8)

## FINAL

scores <- c(
  PredictSNP1_2 = 87, MAPP = 88, PhD_SNP = 88, PolyPhen1 = 74, PolyPhen2 = 81,
  SIFT = 79, SNAP = 85, PANTHER = 87, CADD = 84, DANN = 77, FATHMM = 73, FunSeq2 = 62
)

axes <- names(scores)
red  <- as.numeric(scores)
grey <- rep(50, length(scores))

plot_ly() |>
  add_trace(
    type = "scatterpolar", mode = "lines",
    r = red, theta = axes,
    fill = "toself",
    line = list(color = "rgba(220,20,60,1)", width = 1.5),
    fillcolor = "rgba(220,20,60,0.25)",   # softer red background
    name = "Scores"
  ) |>
  add_trace(
    type = "scatterpolar", mode = "lines",
    r = grey, theta = axes,
    fill = "toself",
    line = list(color = "rgba(40,40,40,1)", width = 1.5),   # darker grey line
    fillcolor = "rgba(60,60,60,0.45)",                    # dark & visible grey fill
    name = "Neutral 50%"
  ) |>
  layout(
    showlegend = FALSE,
    polar = list(
      radialaxis = list(range = c(0, 100), dtick = 10, tickfont = list(size = 10)),
      angularaxis = list(direction = "clockwise")
    )
  )

# Functional analysis - MutPred2 data ----------------------------------------------------------------

# read the sheet

dat <- read_excel("C:/Users/Shadi Shahatit/OneDrive/Desktop/MutPred2_data_41467_2020_19669_MOESM3_ESM.xlsx", sheet = 2)

# gather all property score columns

long <- dat %>%
  select(`Protein ID`,
         starts_with("Property"),
         starts_with("Property score")) %>%
  pivot_longer(
    cols = matches("Property score"),
    names_to = "score_col",
    values_to = "Score") %>%
  mutate(
    Property = rep(c(dat$`Property...5`, dat$`Property...9`, dat$`Property...13`), each = 1)[1:n()],
    Protein  = rep(dat$`Protein ID`, times = 3))

# stats

summary(dat$`MutPred2 score`)
sd(dat$`MutPred2 score`)
quantile(dat$`MutPred2 score`, probs = c(0.989), na.rm = TRUE)
length(unique(dat$`Protein ID`))

# plot

ggplot(long, aes(x = Protein, y = Score)) +
  geom_point(size = 2) +
  # coord_flip() +
  theme_bw() +
  labs(x = "Protein", y = "Property score", title = "MutPred2 Property Scores")

ggplot(long, aes(x = Property, y = Score)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "Property", y = "Score")

ggplot(dat, aes(x = reorder(`Protein ID`, `MutPred2 score`), 
                y = `MutPred2 score`)) +
  geom_point(color = "steelblue", size = 2) +
  # coord_flip() +
  theme_bw() +
  labs(x = "Protein ID", y = "MutPred2 score",
       title = "MutPred2 Scores per Protein")

ggplot(dat, aes(x = (`Protein ID`), 
                y = `MutPred2 score`)) +
  geom_point(color = "steelblue", size = 2) +
  # coord_flip() +
  theme_classic() +
  labs(x = "Protein ID", y = "MutPred2 score",
       title = "MutPred2 Scores per Protein")

ggplot(dat, aes(x = `Protein ID`, y = `MutPred2 score`)) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0.971, color = "red", linetype = "dashed", linewidth = 0.7) +
  theme_classic() +
  labs(x = "Protein ID", y = "MutPred2 score",
       title = "MutPred2 Scores per Protein")

ggplot(dat, aes(x = `Protein ID`, y = `MutPred2 score`)) +
  geom_point(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0.971, color = "red", linetype = "dashed", linewidth = 0.7) +
  theme_classic() +
  theme(
    axis.text.x = element_blank()  # hides labels only
  ) +
  labs(x = "Gene", y = "MutPred2 score",
       title = "")


