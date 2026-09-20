# Prototype; not yet carefully vetted.
# Figure 6A candidate panels from existing q04_network tables.
# Claude Science frame cbb4e1ad-18c1-4945-a745-291b60aef9bb, 2026-09-20.
# Inputs: ../data/*.csv, prepared from
#   revisions/2026-09-08_routeA_execution/q04_network/results/ (see ../data_provenance.md).
# Nothing here refits a network or touches primary data.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else dirname(normalizePath(sys.frame(1)$ofile %||% "."))
data_dir <- file.path(root, "data"); fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, colorblind safe
col_wt   <- "#0072B2"   # WT-associated: relative out-degree rank higher in WT
col_fer4 <- "#D55E00"   # fer-4-associated: relative out-degree rank higher in fer-4
col_scored <- "#D9D9D9"; col_one <- "#F2F2F2"; col_absent <- "white"
base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey30"), legend.position = "bottom",
        legend.text = element_text(size = 9), legend.title = element_text(size = 9))

# ------------------------------------------------------------------ data
M   <- read_csv(file.path(data_dir, "dof_context_status.csv"), show_col_types = FALSE)
per <- read_csv(file.path(data_dir, "dof_summary.csv"), show_col_types = FALSE)
en  <- read_csv(file.path(data_dir, "family_enrichment_paired_universe.csv"), show_col_types = FALSE)
ev  <- read_csv(file.path(data_dir, "tested_DOF_integrated_evidence_matrix.csv"), show_col_types = FALSE)
dirt <- read_csv(file.path(data_dir, "direction_by_tissue.csv"), show_col_types = FALSE)
yt  <- read_csv(file.path(data_dir, "yuc5_regulators_by_tissue.csv"), show_col_types = FALSE)
yu  <- read_csv(file.path(data_dir, "yuc5_universe_by_tissue.csv"), show_col_types = FALSE)
ys  <- fromJSON(file.path(data_dir, "yuc5_summary.json"))

zone_lab <- c("Proliferation Domain" = "P", "Transition Domain" = "T", "Elongation" = "E",
              "Maturation" = "M", "Distal" = "D", "Proximal" = "Px")
tissue_order <- c("Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Xylem Pole Pericycle",
                  "Phloem Pole Pericycle", "Procambium", "Protophloem", "Metaphloem & Companion Cell",
                  "Protoxylem", "Metaxylem", "Lateral Root Cap", "Columella")
tissue_short <- c("Atrichoblast" = "Atrichoblast", "Trichoblast" = "Trichoblast", "Cortex" = "Cortex",
                  "Endodermis" = "Endodermis", "Xylem Pole Pericycle" = "XPP", "Phloem Pole Pericycle" = "PPP",
                  "Procambium" = "Procambium", "Protophloem" = "Protophloem",
                  "Metaphloem & Companion Cell" = "Metaphloem/CC", "Protoxylem" = "Protoxylem",
                  "Metaxylem" = "Metaxylem", "Lateral Root Cap" = "LRC", "Columella" = "Columella")
status_levels <- c("WT-associated", "fer-4-associated", "scored, not nominated",
                   "scored in one genotype only", "not in network")

M <- M %>% mutate(zone_s = factor(zone_lab[zone], levels = zone_lab),
                  tissue_f = factor(tissue_short[tissue], levels = tissue_short[tissue_order]),
                  status = factor(status, levels = status_levels))
# AT3G45610 is DOF3.2 in the evidence matrix and DOF6 in the shared DOF list (q04 identifier note)
ev <- ev %>% mutate(symbol = ifelse(gene == "AT3G45610", "DOF6", symbol))
M   <- M   %>% mutate(symbol = ifelse(gene == "AT3G45610", "DOF6", symbol))
per <- per %>% mutate(symbol = ifelse(gene == "AT3G45610", "DOF6", symbol), tested = tested | gene == "AT3G45610")
tested_syms <- ev$symbol
nominated <- per %>% filter(n_WT + n_fer4 > 0) %>% arrange(desc(n_WT + n_fer4), desc(n_fer4)) %>% pull(symbol)
row_lab <- function(s) { s2 <- ifelse(s == "DOF6", "DOF6 (DOF3.2)", s); ifelse(s %in% tested_syms, paste0(s2, " \u2020"), s2) }

# ------------------------------------------------------------------ Panel 1: family enrichment forest
en2 <- en %>% filter(universe_n >= 10) %>%
  mutate(set = ifelse(genotype == "WT_specific", "WT-associated regulators", "fer-4-associated regulators"),
         is_dof = Family == "Dof", lab = sprintf("%d/%d", foreground_n, universe_n))
fam_order <- en2 %>% filter(genotype == "fer4_specific") %>% arrange(odds_ratio) %>% pull(Family)
en2 <- en2 %>% mutate(Family = factor(Family, levels = fam_order),
                      set = factor(set, levels = c("WT-associated regulators", "fer-4-associated regulators")),
                      sig = p_adj_BH < 0.05)
p_forest <- ggplot(en2, aes(x = odds_ratio, y = Family)) +
  geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +
  geom_errorbarh(aes(xmin = odds_lower, xmax = odds_upper, colour = is_dof), height = 0, linewidth = 0.5) +
  geom_point(aes(colour = is_dof, shape = sig), size = 2.4, fill = "white") +
  geom_text(data = filter(en2, is_dof), aes(label = lab), vjust = -1, size = 3, colour = "#000000") +
  facet_wrap(~set, nrow = 1) +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10)) +
  scale_colour_manual(values = c(`TRUE` = "#CC79A7", `FALSE` = "grey45"), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c("BH p >= 0.05", "BH p < 0.05"), name = NULL) +
  labs(x = "Odds ratio (one-sided Fisher, 95% CI)", y = NULL,
       title = "DOF is the most over-represented TF family among genotype-associated regulators",
       subtitle = sprintf("Universe: %d Kay TFs scored in both genotypes in at least one context; families with >= 10 members; label = DOFs nominated / DOFs in universe",
                          unique(en2$universe_total))) +
  base_theme + theme(panel.grid.major.y = element_blank())

# ------------------------------------------------------------------ Panel 2: DOF x context matrix
Mn <- M %>% filter(symbol %in% nominated) %>%
  mutate(symbol_f = factor(symbol, levels = rev(nominated)))
p_matrix <- ggplot(Mn, aes(x = zone_s, y = symbol_f)) +
  geom_tile(aes(fill = status), colour = "grey75", linewidth = 0.25) +
  geom_point(data = filter(Mn, de_called), shape = 8, size = 1.4, colour = "black") +
  facet_grid(~tissue_f, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = c(`WT-associated` = col_wt, `fer-4-associated` = col_fer4,
                               `scored, not nominated` = col_scored, `scored in one genotype only` = col_one,
                               `not in network` = col_absent), drop = FALSE, name = NULL) +
  scale_y_discrete(labels = row_lab) +
  labs(x = "Developmental zone within tissue (P proliferation, T transition, E elongation, M maturation, D distal, Px proximal)", y = NULL,
       title = "Where the basal networks nominate each DOF (26 of 37 family members, 44 cell type x zone contexts)",
       subtitle = paste("Blue: relative out-degree rank higher in WT than in fer-4 (WT-associated). Orange: higher in fer-4 (fer-4-associated).",
                        "Ranks are relative within a context and are not TF activity. Shown: the 26 DOFs nominated at least once.",
                        "* basal fer-4 vs WT DE call in that context (repaired current table). \u2020 tested in Fig. 6c to 6f / ED Fig. 10", sep = "\n")) +
  base_theme +
  theme(strip.placement = "outside", strip.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 8, face = "italic"),
        panel.spacing.x = unit(1.5, "pt"), panel.grid = element_blank(), legend.key.size = unit(10, "pt"))

cand1 <- p_forest / p_matrix + plot_layout(heights = c(1, 1.55)) +
  plot_annotation(tag_levels = "i", title = "Candidate 1 (preferred design): network nomination of DOFs, by family and by cellular context",
                  subtitle = "Exploratory sheet at full resolution; a drop-in 6A would use panel i restricted to significant families (see fig6a_candidate1_compact) plus a tissue-collapsed panel ii",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
ggsave(file.path(fig_dir, "fig6a_candidate1_enrichment_context.pdf"), cand1, width = 11, height = 12, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig6a_candidate1_enrichment_context.png"), cand1, width = 11, height = 12, dpi = 200, bg = "white")

# ------------------------------------------------------------------ Panel 3: DOF ledger (Venn replacement)
n_ctx <- n_distinct(M$context)
led <- M %>% count(symbol, status, name = "n") %>%
  complete(symbol, status, fill = list(n = 0)) %>%
  left_join(per %>% select(symbol, tested, n_DE), by = "symbol")
ord <- per %>% arrange(n_WT + n_fer4, n_scored) %>% pull(symbol)
led <- led %>% mutate(symbol_f = factor(symbol, levels = ord))
# Three statuses per assay: positive = reported positive result; tested = assayed, no reported induction/effect;
# not tested = not in the submitted package. Coverage and outcome are kept distinct.
st <- function(x, pos_pat, tested_pat) case_when(grepl(pos_pat, x) ~ "positive", grepl(tested_pat, x) ~ "tested, no reported induction", TRUE ~ "not tested")
marks <- ev %>% transmute(symbol,
  `pYUC5-LUC activation (protoplast)` = "positive",
  `RALF23 induction, root qPCR` = st(root_RALF_qPCR, "^WT induction shown", "No significant WT induction"),
  `Binds YUC5 promoter (EMSA)` = st(EMSA, "shown", "^Tested"),
  `Mutant RALF23 growth assay` = st(mutant_RALF_growth, "^Tested", "^Assayed"),
  `Endogenous YUC5 in mutant` = st(endogenous_YUC5_in_DOF_mutant, "shown", "^Tested")) %>%
  pivot_longer(-symbol, names_to = "assay", values_to = "status3") %>%
  mutate(status3 = factor(status3, levels = c("positive", "tested, no reported induction", "not tested")),
         done = status3 != "not tested") %>%
  mutate(symbol_f = factor(symbol, levels = ord),
         assay = factor(assay, levels = c("pYUC5-LUC activation (protoplast)", "RALF23 induction, root qPCR",
                                          "Binds YUC5 promoter (EMSA)", "Mutant RALF23 growth assay", "Endogenous YUC5 in mutant")))
p_led <- ggplot(led, aes(y = symbol_f, x = n, fill = status)) +
  geom_col(width = 0.8, colour = "grey60", linewidth = 0.2) +
  scale_fill_manual(values = c(`WT-associated` = col_wt, `fer-4-associated` = col_fer4,
                               `scored, not nominated` = col_scored, `scored in one genotype only` = col_one,
                               `not in network` = col_absent), drop = FALSE, name = NULL) +
  scale_y_discrete(labels = row_lab) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, n_ctx, 11)) +
  labs(x = sprintf("Number of network contexts (of %d)", n_ctx), y = NULL,
       title = "All 37 coauthor-listed DOFs: network status per context",
       subtitle = "26 nominated in at least one context; 11 never nominated\n(DOF4.2, DOF4.3, DOF4.4 and AT1G65935 absent from every network)") +
  base_theme + theme(axis.text.y = element_text(size = 8, face = "italic"), panel.grid.major.y = element_blank())
p_marks <- ggplot(marks, aes(x = assay, y = symbol_f)) +
  geom_tile(fill = "white", colour = "grey80") +
  geom_point(aes(shape = status3), size = 2.4, colour = "black", fill = "white") +
  scale_shape_manual(values = c(positive = 16, `tested, no reported induction` = 21, `not tested` = 4), name = NULL, drop = FALSE) +
  scale_y_discrete(labels = row_lab, drop = FALSE, limits = ord) +
  labs(x = NULL, y = NULL, title = "Assays in this paper",
       subtitle = "blank row = not selected\nmutant growth = coverage only") +
  base_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_blank(),
                     panel.grid = element_blank(), legend.direction = "vertical")
cand2 <- p_led + p_marks + plot_layout(widths = c(3, 1.2)) +
  plot_annotation(title = "Candidate 2 (Venn replacement): one row per DOF, network status and assay coverage",
                  subtitle = "Exploratory sheet; for a main-figure panel keep the 26 nominated rows or the 10 tested rows",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
ggsave(file.path(fig_dir, "fig6a_candidate2_dof_ledger.pdf"), cand2, width = 11, height = 9, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig6a_candidate2_dof_ledger.png"), cand2, width = 11, height = 9, dpi = 200, bg = "white")

# ------------------------------------------------------------------ Panel 4: additional network insight (Extended Data candidate)
# (i) direction of DOF nominations per tissue against the other-TF reference
d2 <- dirt %>% mutate(group = ifelse(is_dof, "DOF", "Other TF"),
                      tissue_f = factor(tissue_short[tissue], levels = tissue_short[tissue_order])) %>%
  select(tissue_f, group, WT_specific, fer4_specific) %>%
  pivot_longer(c(WT_specific, fer4_specific), names_to = "dir", values_to = "n") %>%
  mutate(dir = ifelse(dir == "WT_specific", "WT-associated", "fer-4-associated"),
         dir = factor(dir, levels = c("WT-associated", "fer-4-associated")))
tot <- d2 %>% group_by(tissue_f, group) %>% summarise(N = sum(n), .groups = "drop")
d2 <- d2 %>% left_join(tot, by = c("tissue_f", "group")) %>% mutate(frac = ifelse(N > 0, n / N, NA))
p_dir <- ggplot(d2, aes(x = group, y = frac, fill = dir)) +
  geom_col(width = 0.75, colour = "grey50", linewidth = 0.2) +
  geom_text(data = tot, aes(x = group, y = 1.04, label = paste0("n=", N)), inherit.aes = FALSE, size = 2.6) +
  facet_grid(~tissue_f, switch = "x") +
  scale_fill_manual(values = c(`WT-associated` = col_wt, `fer-4-associated` = col_fer4), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Share of nominations",
       title = "Direction of nominations by tissue: endodermal DOF records mostly rank higher in WT",
       subtitle = paste("Endodermis: 18 WT-associated vs 6 fer-4-associated DOF records (14 and 5 distinct DOFs); other TFs 144 vs 159.",
                        "Records are TF-by-context rows from one metric and are not independent replicates; no significance test is implied.",
                        "Same 18:6 count in all nine top-N/tail settings and after dropping the two affected vascular fits (coordinator check, 2026-09-20).", sep = "\n")) +
  base_theme + theme(strip.placement = "outside", strip.text.x = element_text(size = 8, angle = 90, hjust = 1),
                     axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5), panel.spacing.x = unit(2, "pt"))

# (ii) YUC5 regulator view with explicit non-estimable tissues
yu2 <- yu %>% filter(tissue != "all") %>%
  mutate(tissue_l = recode(tissue, atrichoblast = "Atrichoblast", trichoblast = "Trichoblast", cortex = "Cortex",
                           endodermis = "Endodermis", pericycle = "Pericycle", phloem = "Phloem", procambium = "Procambium",
                           xylem = "Xylem", lrc = "LRC", columella = "Columella")) %>%
  left_join(yt, by = "tissue") %>%
  mutate(uniq_src = replace_na(uniq_src, 0), uniq_dof_src = replace_na(uniq_dof_src, 0),
         other = uniq_src - uniq_dof_src,
         note = case_when(!YUC5_in_target_universe ~ "YUC5 not in target set (not estimable)",
                          uniq_src == 0 ~ "YUC5 in target set, no retained edge",
                          TRUE ~ sprintf("%d regulators, %d DOF", uniq_src, uniq_dof_src)))
yl <- yu2 %>% select(tissue_l, DOF = uniq_dof_src, `Other TF` = other, note, YUC5_in_target_universe) %>%
  pivot_longer(c(DOF, `Other TF`), names_to = "cls", values_to = "n") %>%
  mutate(tissue_l = factor(tissue_l, levels = c("Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle",
                                                "Phloem", "Procambium", "Xylem", "LRC", "Columella")))
p_yuc <- ggplot(yl, aes(x = tissue_l, y = n, fill = cls)) +
  geom_col(width = 0.7, colour = "grey50", linewidth = 0.2) +
  geom_text(data = distinct(yl, tissue_l, note, YUC5_in_target_universe) %>% filter(!YUC5_in_target_universe),
            aes(x = tissue_l, y = 4, label = note), inherit.aes = FALSE, angle = 90, hjust = 0, size = 2.6, colour = "grey20") +
  geom_text(data = yu2 %>% filter(YUC5_in_target_universe), aes(x = tissue_l, y = uniq_src + 6, label = note),
            inherit.aes = FALSE, size = 2.7, colour = "grey10") +
  scale_fill_manual(values = c(DOF = "#CC79A7", `Other TF` = "grey60"), name = "Retained YUC5 regulators") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, 250)) +
  labs(x = NULL, y = "Unique TFs with a retained YUC5 edge",
       title = "Fitted YUC5 regulators exist only in stele networks; DOF share is no higher than in the prior",
       subtitle = sprintf("Retained: %d DOF of %d unique YUC5 regulators (5.8%%), a subset of the supplied prior (%d DOF of %d, 5.0%%). Descriptive percentages; nested sets, no test.\nCortex, endodermis, epidermis and cap networks never contained YUC5 as a target, so a DOF to YUC5 edge there is not estimable, not absent.",
                          ys$retained_DOF_sources, ys$retained_unique_sources, ys$prior_DOF, ys$prior_YUC5_regulators)) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
cand3 <- p_dir / p_yuc + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "i", title = "Candidate 3 (Extended Data sheet): two additional readouts of the existing network framework",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
ggsave(file.path(fig_dir, "fig6_network_insight_candidate3.pdf"), cand3, width = 11, height = 10, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig6_network_insight_candidate3.png"), cand3, width = 11, height = 10, dpi = 200, bg = "white")


# ------------------------------------------------------------------ Compact panel: families significant in either set
sig_fams <- en2 %>% filter(sig) %>% pull(Family) %>% unique() %>% as.character()
enc <- en2 %>% filter(Family %in% sig_fams) %>% mutate(Family = factor(Family, levels = fam_order[fam_order %in% sig_fams]),
                                                      set = factor(set, levels = c("WT-associated regulators", "fer-4-associated regulators")))
p_compact <- ggplot(enc, aes(x = odds_ratio, y = Family, colour = set, group = set)) +
  geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +
  geom_errorbarh(aes(xmin = odds_lower, xmax = odds_upper), height = 0, linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = sig), size = 2.8, fill = "white", position = position_dodge(width = 0.6)) +
  geom_text(aes(label = lab, x = odds_upper * 1.15), position = position_dodge(width = 0.6), size = 2.7, hjust = 0, show.legend = FALSE) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8, 16), limits = c(0.5, 40)) +
  scale_colour_manual(values = c(`WT-associated regulators` = col_wt, `fer-4-associated regulators` = col_fer4), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c("BH p >= 0.05", "BH p < 0.05"), name = NULL) +
  labs(x = "Odds ratio (95% CI)", y = NULL,
       title = "TF families over-represented among\ngenotype-associated network regulators",
       subtitle = sprintf("Families with BH p < 0.05 in either set\nlabel = nominated / in universe (%d Kay TFs)", unique(en2$universe_total))) +
  base_theme + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"),
                     plot.title = element_text(size = 10), plot.subtitle = element_text(size = 7.5), panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "fig6a_candidate1_compact.pdf"), p_compact, width = 4.2, height = 3.6, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig6a_candidate1_compact.png"), p_compact, width = 4.2, height = 3.6, dpi = 300, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(root, "code", "sessionInfo.txt"))
cat("done\n")
