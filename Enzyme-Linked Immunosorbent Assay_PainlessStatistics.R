# ======================================================================================
# Author: Mateusz Glenszczyk
# Email: mateusz.glenszczyk@gmail.com
# Date: 2026-03-19
# Description: ELISA-Pipeline
# ======================================================================================

#=======================================================================================
# My PhD script to quickly roll-over the pipeline of ELISA results analysis.
# For this script to work, you first need to have another R file with data matrix.
# My data matrices did not have the structure (e.g., 12 x 12), since I was adding
# more and more results on the go. Therefore less tweaks are required to adapt the code.
# Code may have some parts of polish language. Sorry for that! :D
#=======================================================================================

# =========================================
# LIBRARIES
# =========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(grid)
  library(scales)
})

source("E:/Naukowe/Awaria Laptopa - Kopia Zapasowa TU PRACUJEMY/Publikacyjne/Wersje/Artykuł 4 - Mikro i Fizjo/Supplementary Material/Codes/protein_matrices.R")

# =========================================
# THINGIES TO TWEAK
# =========================================
TEST <- "wilcox"     # "wilcox" albo "t"
ADJ  <- "holm"       # korekta wielokrotnych porównań
ALT  <- "greater"    # test jednostronny: czy grupa > BLANK

UNIT_FACTOR <- 10    # konwersja z ug/100 ul -> ug/ml

SPIDER_COLS <- c(
  PRD  = "#274DF5",
  PTSD = "#0E7005"
)

SPIDER_SHAPES <- c(
  PRD  = 21,
  PTSD = 21
)

SPIDER_LABELS <- c(
  PRD  = expression(italic("Pardosa lugubris")),
  PTSD = expression(italic("Parasteatoda tepidariorum"))
)

SAMPLE_LEVELS <- c("EGGS", "SILK")
SAMPLE_LABELS <- c(
  EGGS = "EGGS",
  SILK = "SILK"
)

PROTEIN_LABELS <- c(
  LIZO  = "Lysozyme-like peptides",
  DEFE  = "Defensin-like peptides",
  CECRO = "Cecropin-like peptides"
)

#Y_LABEL_RAW  <- expression("Blank corrected protein concentration (µg/ml)")
#Y_LABEL_CORR <- expression("Protein concentration above blank (µg/ml)")

Y_LABEL_RAW  <- expression("Blank-corrected protein concentration (µg/mL"^{-1}*")")
Y_LABEL_CORR <- expression("Blank-corrected protein concentration (µg/mL"^{-1}*")")

p_to_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# =========================================
# LONG FORMAT
# =========================================
blank_vals <- as.numeric(BLANK) * UNIT_FACTOR

make_df <- function(x, object_name){
  tibble(
    value = as.numeric(x) * UNIT_FACTOR,
    object = object_name
  ) |>
    mutate(
      spider = case_when(
        str_detect(object, "^PRDE_|^PRDS_")   ~ "PRD",
        str_detect(object, "^PTSDE_|^PTSDS_") ~ "PTSD",
        TRUE ~ NA_character_
      ),
      sample = case_when(
        str_detect(object, "^PRDE_|^PTSDE_") ~ "EGGS",
        str_detect(object, "^PRDS_|^PTSDS_") ~ "SILK",
        TRUE ~ NA_character_
      ),
      protein = case_when(
        str_detect(object, "_LIZO$")  ~ "LIZO",
        str_detect(object, "_DEFE$")  ~ "DEFE",
        str_detect(object, "_CECRO$") ~ "CECRO",
        TRUE ~ NA_character_
      )
    ) |>
    select(value, sample, spider, protein)
}

protein_objects <- list(
  PRDE_LIZO   = PRDE_LIZO,
  PRDE_DEFE   = PRDE_DEFE,
  PRDE_CECRO  = PRDE_CECRO,
  PRDS_LIZO   = PRDS_LIZO,
  PRDS_DEFE   = PRDS_DEFE,
  PRDS_CECRO  = PRDS_CECRO,
  PTSDE_LIZO  = PTSDE_LIZO,
  PTSDE_DEFE  = PTSDE_DEFE,
  PTSDE_CECRO = PTSDE_CECRO,
  PTSDS_LIZO  = PTSDS_LIZO,
  PTSDS_DEFE  = PTSDS_DEFE,
  PTSDS_CECRO = PTSDS_CECRO
)

df <- bind_rows(
  lapply(names(protein_objects), function(nm) make_df(protein_objects[[nm]], nm))
) |>
  mutate(
    sample  = factor(sample, levels = SAMPLE_LEVELS),
    spider  = factor(spider, levels = c("PRD", "PTSD")),
    protein = factor(protein, levels = c("LIZO", "DEFE", "CECRO"))
  )

# =========================================
# BLANK: poziom tła metody po przeliczeniu
# =========================================
blank_ref <- tibble(
  protein = factor(c("LIZO", "DEFE", "CECRO"), levels = c("LIZO", "DEFE", "CECRO")),
  blank_mean   = mean(blank_vals, na.rm = TRUE),
  blank_median = median(blank_vals, na.rm = TRUE),
  blank_sd     = sd(blank_vals, na.rm = TRUE),
  blank_q1     = quantile(blank_vals, 0.25, na.rm = TRUE),
  blank_q3     = quantile(blank_vals, 0.75, na.rm = TRUE)
)

# =========================================
# TESTY: EGGS/SILK vs BLANK
# jednostronnie: grupa > BLANK
# =========================================
compare_vs_blank <- function(d, ctrl_vals, test = "wilcox", alt = "greater"){
  d <- d |>
    filter(!is.na(value))
  
  groups <- d |>
    distinct(sample) |>
    pull(sample)
  
  res <- lapply(groups, function(smp){
    x <- d |>
      filter(sample == smp) |>
      pull(value)
    
    if (length(x) < 1 || length(ctrl_vals) < 1) return(NULL)
    
    if (test == "t") {
      tt <- t.test(x, ctrl_vals, alternative = alt)
      tibble(
        comparison = paste0(as.character(smp), " vs BLANK"),
        test = paste0("one-sided t-test (", alt, ")"),
        statistic = unname(tt$statistic),
        p_value = tt$p.value,
        n_treat = length(x),
        n_blank = length(ctrl_vals),
        median_treat = median(x),
        median_blank = median(ctrl_vals),
        mean_treat = mean(x),
        mean_blank = mean(ctrl_vals),
        diff_median = median(x) - median(ctrl_vals),
        diff_mean = mean(x) - mean(ctrl_vals)
      )
    } else {
      wt <- wilcox.test(x, ctrl_vals, alternative = alt, exact = FALSE)
      tibble(
        comparison = paste0(as.character(smp), " vs BLANK"),
        test = paste0("one-sided Wilcoxon rank-sum (", alt, ")"),
        statistic = unname(wt$statistic),
        p_value = wt$p.value,
        n_treat = length(x),
        n_blank = length(ctrl_vals),
        median_treat = median(x),
        median_blank = median(ctrl_vals),
        mean_treat = mean(x),
        mean_blank = mean(ctrl_vals),
        diff_median = median(x) - median(ctrl_vals),
        diff_mean = mean(x) - mean(ctrl_vals)
      )
    }
  })
  
  bind_rows(res)
}

tests_vs_blank <- df |>
  group_by(protein, spider) |>
  group_modify(~ compare_vs_blank(.x, ctrl_vals = blank_vals, test = TEST, alt = ALT)) |>
  ungroup() |>
  group_by(protein, spider) |>
  mutate(
    p_adj = p.adjust(p_value, method = ADJ),
    sig   = p_to_stars(p_adj)
  ) |>
  ungroup() |>
  arrange(protein, spider, comparison)

tests_vs_blank

# =========================================
# TABELA WYNIKÓW
# =========================================
tests_table <- tests_vs_blank |>
  mutate(
    sample = gsub(" vs BLANK", "", comparison),
    sample = factor(sample, levels = c("EGGS", "SILK")),
    sample_label = recode(as.character(sample), !!!SAMPLE_LABELS),
    spider_label = recode(as.character(spider),
                          PRD  = "Pardosa lugubris",
                          PTSD = "Parasteatoda tepidariorum"),
    protein_label = recode(as.character(protein),
                           LIZO  = "Lysozyme peptides",
                           DEFE  = "Defensins peptides",
                           CECRO = "Cecropin-like peptides"),
    direction = case_when(
      diff_median > 0 ~ "higher than blank",
      diff_median < 0 ~ "lower than blank",
      TRUE ~ "equal to blank"
    )
  ) |>
  arrange(protein, spider, sample) |>
  transmute(
    protein = protein_label,
    spider = spider_label,
    comparison = paste0(sample_label, " vs BLANK"),
    test,
    n_blank,
    n_treat,
    median_blank = round(median_blank, 3),
    median_treat = round(median_treat, 3),
    diff_median = round(diff_median, 3),
    mean_blank = round(mean_blank, 3),
    mean_treat = round(mean_treat, 3),
    diff_mean = round(diff_mean, 3),
    direction,
    statistic = round(statistic, 3),
    p_value = signif(p_value, 4),
    p_adj = signif(p_adj, 4),
    sig
  )

tests_table

write.csv(
  tests_table,
  "protein_statistics_vs_blank_one_sided.csv",
  row.names = FALSE
)

# =========================================
# WYKRES 1: SUROWE WARTOSCI PO PRZELICZENIU. 
# LINIA = MEDIANA BLANKU.
# PASMO = IQR BLANKU.
# =========================================
y_max_data  <- max(df$value, na.rm = TRUE)
y_max_blank <- max(blank_ref$blank_q3, na.rm = TRUE)
y_upper     <- max(y_max_data, y_max_blank) + 1.0
y_lower     <- min(0, min(df$value, na.rm = TRUE) - 0.5)

dodge_w  <- 0.78
jitter_w <- 0.18

p_raw <- ggplot(df, aes(x = sample, y = value, fill = spider)) +
  geom_rect(
    data = blank_ref,
    aes(xmin = -Inf, xmax = Inf, ymin = blank_q1, ymax = blank_q3),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.6
  ) +
  geom_hline(
    data = blank_ref,
    aes(yintercept = blank_median),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey30"
  ) +
  geom_boxplot(
    aes(group = interaction(sample, spider)),
    position = position_dodge(width = dodge_w),
    width = 0.62,
    notch = FALSE,
    outlier.shape = NA,
    alpha = 0.90,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_point(
    aes(shape = spider),
    position = position_jitterdodge(
      jitter.width = jitter_w,
      dodge.width  = dodge_w
    ),
    size = 2.5,
    alpha = 0.45,
    colour = "black",
    stroke = 0.35
  ) +
  facet_wrap(
    ~ protein,
    nrow = 1,
    scales = "fixed",
    labeller = as_labeller(PROTEIN_LABELS)
  ) +
  scale_x_discrete(labels = SAMPLE_LABELS, drop = FALSE) +
  scale_y_continuous(
    limits = c(y_lower, y_upper),
    breaks = pretty_breaks(n = 5),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  scale_fill_manual(values = SPIDER_COLS, labels = SPIDER_LABELS, drop = FALSE) +
  scale_shape_manual(values = SPIDER_SHAPES, labels = SPIDER_LABELS, drop = FALSE) +
  labs(
    x = NULL,
    y = Y_LABEL_RAW,
    fill = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.ticks.length = unit(3, "pt"),
    legend.position = "bottom",
    plot.margin = margin(8, 10, 8, 10),
    panel.spacing = unit(8, "pt")
  ) +
  guides(
    fill = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

stars_df <- df |>
  group_by(protein, spider, sample) |>
  summarise(y = max(value, na.rm = TRUE), .groups = "drop") |>
  left_join(
    tests_vs_blank |>
      mutate(sample = gsub(" vs BLANK", "", comparison)) |>
      select(protein, spider, sample, sig),
    by = c("protein", "spider", "sample")
  ) |>
  mutate(y = y + 0.5)

final_plot_raw <- p_raw +
  geom_text(
    data = stars_df,
    aes(x = sample, y = y, label = sig, group = spider),
    position = position_dodge(width = dodge_w),
    vjust = 0,
    size = 4
  )

print(final_plot_raw)

ggsave(
  filename = "protein_concentration_raw_with_blank_line.tiff",
  plot     = final_plot_raw,
  width    = 12,
  height   = 6.5,
  dpi      = 600,
  units    = "in",
  compression = "lzw"
)

# =========================================
# WYKRES 2: KORYGOWANIE DANYCH O BLANK.
# value_corr = value - median(BLANK)
# =========================================
blank_median_global <- median(blank_vals, na.rm = TRUE)

df_corr <- df |>
  mutate(value_corr = value - blank_median_global)

y2_max <- max(df_corr$value_corr, na.rm = TRUE) + 1.0
y2_min <- min(df_corr$value_corr, na.rm = TRUE) - 0.5

p_corr <- ggplot(df_corr, aes(x = sample, y = value_corr, fill = spider)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey30"
  ) +
  geom_boxplot(
    aes(group = interaction(sample, spider)),
    position = position_dodge(width = dodge_w),
    width = 0.62,
    notch = FALSE,
    outlier.shape = NA,
    alpha = 0.90,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_point(
    aes(shape = spider),
    position = position_jitterdodge(
      jitter.width = jitter_w,
      dodge.width  = dodge_w
    ),
    size = 2.5,
    alpha = 0.45,
    colour = "black",
    stroke = 0.35
  ) +
  facet_wrap(
    ~ protein,
    nrow = 1,
    scales = "fixed",
    labeller = as_labeller(PROTEIN_LABELS)
  ) +
  scale_x_discrete(labels = SAMPLE_LABELS, drop = FALSE) +
  scale_y_continuous(
    limits = c(y2_min, y2_max),
    breaks = pretty_breaks(n = 5),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  scale_fill_manual(values = SPIDER_COLS, labels = SPIDER_LABELS, drop = FALSE) +
  scale_shape_manual(values = SPIDER_SHAPES, labels = SPIDER_LABELS, drop = FALSE) +
  labs(
    x = NULL,
    y = Y_LABEL_CORR,
    fill = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.ticks.length = unit(3, "pt"),
    legend.position = "bottom",
    plot.margin = margin(8, 10, 8, 10),
    panel.spacing = unit(8, "pt")
  ) +
  guides(
    fill = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

stars_df_corr <- df_corr |>
  group_by(protein, spider, sample) |>
  summarise(y = max(value_corr, na.rm = TRUE), .groups = "drop") |>
  left_join(
    tests_vs_blank |>
      mutate(sample = gsub(" vs BLANK", "", comparison)) |>
      select(protein, spider, sample, sig),
    by = c("protein", "spider", "sample")
  ) |>
  mutate(y = y + 0.5)

final_plot_corr <- p_corr +
  geom_text(
    data = stars_df_corr,
    aes(x = sample, y = y, label = sig, group = spider),
    position = position_dodge(width = dodge_w),
    vjust = 0,
    size = 4
  )

print(final_plot_corr)

ggsave(
  filename = "protein_concentration_blank_corrected.tiff",
  plot     = final_plot_corr,
  width    = 12,
  height   = 6.5,
  dpi      = 600,
  units    = "in",
  compression = "lzw"
)

# =========================================
# PODSUMOWANIE RZECZY
# =========================================
summary_sig <- tests_vs_blank |>
  mutate(
    sample = gsub(" vs BLANK", "", comparison),
    direction = case_when(
      diff_median > 0 ~ "higher than blank",
      diff_median < 0 ~ "lower than blank",
      TRUE ~ "equal to blank"
    )
  ) |>
  select(
    protein, spider, sample,
    median_blank, median_treat, diff_median,
    p_value, p_adj, sig, direction
  ) |>
  arrange(protein, spider, sample)

summary_sig

