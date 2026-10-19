#!/usr/bin/env Rscript
## 02 — Fig. 6a: TF-family forest plot. Needs 01 to have run (reads TF_family_enrichment_all_tests.csv).
# 20260923 — Panel A: TF-family enrichment among genotype-specific TFs, ANALYSIS 1 (gene level).
# Odds ratio on a LOG scale (a ratio scale must be log or 4x and 1/4x look unequal), 95% CI whiskers,
# reference line at OR = 1, nominated/eligible on the right, DOF bold.
# 20260924 (user): genotype is encoded by COLOUR — WT blue, fer-4 orange. Significance is encoded by
# FILL: filled = significant (BH q < 0.05), open = not.
# 20261001 (user): ONE COMBINED PANEL to save space, in the format of Trevor's
# from_Trevor/combined_forest_prototype.pdf — WT and fer-4 share each family row (WT above, fer-4
# below), alternate rows shaded, two coloured nominated/eligible columns under a header on the right,
# colour + fill keys underneath. Replaces the two side-by-side genotype panels.
# Background = the 715 TFs that reach the top 100 by degree_centrality_out in >=1 genotype (the only
# genes that could be nominated). Families with >=10 eligible TFs; BH applied across ALL families x
# both genotypes upstream (84 tests), the >=10 filter AFTER (see 20260923_DOF_enrichment_report.R).
# ⚠ The WT and fer-4 sets overlap (16 of 29 DOFs are in both) — the user handles that in the legend.
suppressMessages({library(ggplot2);library(dplyr);library(patchwork);library(scales)})
source("00_paths.R")
COL_WT  <- "#2b6cb0"   # blue  = WT
COL_FER <- "#dd6b20"   # orange = fer-4
res <- read.csv(file.path(OUT_DIR,"TF_family_enrichment_all_tests.csv"), stringsAsFactors=FALSE) %>%
       filter(eligible >= 10)                 # display filter, applied after the BH correction in 01

# top 10 families by the larger of the two odds ratios (PI's selection rule)
top10 <- res %>% group_by(family) %>% summarise(mx = max(OR), .groups="drop") %>%
         arrange(desc(mx)) %>% head(10) %>% pull(family)
d <- res %>% filter(family %in% top10) %>%
     mutate(genotype = factor(genotype, levels=c("WT","fer-4")),
            lab = paste0(nominated,"/",eligible),
            sig = q < 0.05)
# 2026-10-02 (PI): background = the 1,292 TFs that could have entered the ranking. On this background an
# OR of 0 occurs (a family with no nominated member). A log axis cannot hold 0, so -- as in Trevor's
# prototype -- such a row is drawn as a LEFT-POINTING ARROW labelled "0" and gets no point or whisker.
# An infinite OR is not handled (no case yet): fail loudly rather than silently drop it.
stopifnot(all(is.finite(d$OR)), all(is.finite(d$hi)))
d$zero <- d$OR == 0
ord <- d %>% filter(genotype=="fer-4") %>% arrange(OR) %>% pull(family)   # bottom -> top
N   <- length(ord)
d$row <- match(d$family, ord)
DODGE <- 0.21
PTS   <- 1.8     # point size; 2.3 made the WT and fer-4 circles in one row overlap (user 2026-10-01)
d$y <- d$row + ifelse(d$genotype=="WT", DODGE, -DODGE)                      # WT above, fer-4 below
famlab <- gsub("_", " ", ord)                                               # MIKC_MADS -> MIKC MADS
faces  <- ifelse(ord=="Dof", "bold", "plain")
bands  <- data.frame(row = seq(N, 1, by=-2))                                # shade every other row from the top

YLIM <- c(0.5, N + 0.5)
# x range: the largest whisker is DOF fer-4 hi = 13.05, so 20 leaves headroom without the empty
# 16-62 stretch the two-panel version carried. ⚠ bands use XL, not -Inf/Inf: on a log axis +-Inf
# rects silently did not draw.
XL <- c(.15, 20)
# widen if a CI ever needs it (only non-zero ORs have a plotted whisker); keep the fixed breaks
XL <- c(min(XL[1], 0.8*min(d$lo[!d$zero])), max(XL[2], 1.25*max(d$hi[!d$zero])))
ARX <- c(0.36, XL[1]*1.12)   # zero-OR arrow: from ARX[1] leftwards to just inside the axis start
TOPM <- 20      # top margin (pt) in BOTH panels: room for the column header, and keeps rows aligned

# ---- forest panel ----
pF <- ggplot(d, aes(OR, y)) +
  geom_rect(data=bands, aes(xmin=XL[1], xmax=XL[2], ymin=row-.5, ymax=row+.5), inherit.aes=FALSE,
            fill="grey94", colour=NA) +
  geom_vline(xintercept=1, linetype=2, colour="grey55", linewidth=.35) +
  geom_errorbarh(data=filter(d, !zero), aes(xmin=lo, xmax=hi, colour=genotype), height=0, linewidth=.5) +
  # OR = 0: left arrow in the genotype colour (2026-10-02 user: NO "0" label above it)
  geom_segment(data=filter(d, zero), aes(x=ARX[1], xend=ARX[2], y=y, yend=y, colour=genotype),
               linewidth=.5, arrow=arrow(length=unit(1.3,"mm"), type="closed")) +
  # fill is a CONSTANT per layer, not an aesthetic: colour carries genotype, fill only significance
  geom_point(data=filter(d, !sig & !zero), aes(colour=genotype), shape=21, size=PTS, fill="white", stroke=.6) +
  geom_point(data=filter(d,  sig & genotype=="WT"),    shape=21, size=PTS, fill=COL_WT,  colour=COL_WT,  stroke=.6) +
  geom_point(data=filter(d,  sig & genotype=="fer-4"), shape=21, size=PTS, fill=COL_FER, colour=COL_FER, stroke=.6) +
  scale_colour_manual(values=c("WT"=COL_WT, "fer-4"=COL_FER), guide="none") +
  # ⚠ lower limit must sit BELOW the smallest CI bound (Trihelix fer-4 lo = 0.19) or ggplot
  # silently DROPS the whisker -- it warned "Removed 1 row" on the first build.
  scale_x_log10(limits=XL, breaks=c(.25,1,4,16), labels=c("0.25","1","4","16"),
                expand=expansion(0,0)) +
  scale_y_continuous(limits=YLIM, breaks=seq_len(N), labels=famlab, expand=expansion(0,0)) +
  labs(x="Odds ratio (95% CI)", y=NULL) +
  theme_classic(base_size=8) +
  theme(axis.text.y=element_text(face=faces, size=7.6, colour="grey10"),
        axis.text.x=element_text(size=7.2),
        axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),      # user 20260924: no horizontal grid lines (the bands are row shading)
        plot.margin=margin(TOPM,2,3,3))

# ---- nominated / eligible columns, on the same y rows ----
num <- d %>% mutate(x = ifelse(genotype=="WT", 1, 2))
pN <- ggplot(num, aes(x, row)) +
  geom_text(aes(label=lab, colour=genotype), size=2.6) +
  scale_colour_manual(values=c("WT"=COL_WT, "fer-4"=COL_FER), guide="none") +
  annotate("text", x=1.5, y=N+1.15, label="Nominated / eligible", size=2.5, colour="grey30") +
  annotate("text", x=1, y=N+0.62, label="WT", size=2.6, colour=COL_WT) +
  # the mutant allele is always ITALIC (user 2026-09-25)
  annotate("text", x=2, y=N+0.62, label='italic("fer-4")', parse=TRUE, size=2.6, colour=COL_FER) +
  scale_x_continuous(limits=c(.45,2.55), expand=expansion(0,0)) +
  # ⚠ ylim must go in coord_cartesian, NOT scale limits: the header sits ABOVE the rows, and scale
  # limits would DROP it ("Removed 1 row") instead of just drawing it outside the panel.
  scale_y_continuous(expand=expansion(0,0)) +
  coord_cartesian(ylim=YLIM, clip="off") +
  theme_void() + theme(plot.margin=margin(TOPM,10,3,0))

# ---- bottom strip: colour key and fill key ----
# ⚠⚠ dot size and text size are ABSOLUTE (pt/mm) while x positions are FRACTIONS of the width, so an
# offset tuned at one width breaks at another (the 2026-09-25 legend overlapped at 110 mm). The strip
# is rebuilt PER OUTPUT WIDTH and every offset scales as 1/width via `k`. Two rows, two items each,
# sharing two text left edges, so the keys line up as columns.
TXTX <- 0.07; COL2 <- 0.35; SZ <- 2.45   # COL2 widened 0.27 -> 0.35 for the longer "Filled: BH q < 0.05" (2026-10-01)
Y1 <- 0.78; Y2 <- 0.22
build <- function(figw) {
  k   <- 110/figw
  gap <- 0.030 * k                 # text left edge -> dot centre, constant in mm
  x1  <- TXTX; x2 <- TXTX + COL2 * k
  strip <- ggplot() +
    geom_point(aes(x=x1-gap, y=Y1), shape=21, size=PTS, fill=COL_WT,  colour=COL_WT,  stroke=.6) +
    annotate("text", x=x1, y=Y1, label="WT-specific", hjust=0, size=SZ, colour="grey20") +
    geom_point(aes(x=x2-gap, y=Y1), shape=21, size=PTS, fill=COL_FER, colour=COL_FER, stroke=.6) +
    annotate("text", x=x2, y=Y1, label='italic("fer-4")*"-specific"', parse=TRUE, hjust=0, size=SZ, colour="grey20") +
    geom_point(aes(x=x1-gap, y=Y2), shape=21, size=PTS, fill="grey25", colour="grey20", stroke=.6) +
    annotate("text", x=x1, y=Y2, label="Filled: BH q < 0.05", hjust=0, size=SZ, colour="grey25") +
    geom_point(aes(x=x2-gap, y=Y2), shape=21, size=PTS, fill="white", colour="grey20", stroke=.6) +
    annotate("text", x=x2, y=Y2, label="Open: BH q \u2265 0.05", hjust=0, size=SZ, colour="grey25") +
    scale_x_continuous(limits=c(0,1), expand=expansion(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=expansion(0,0)) +
    coord_cartesian(clip="off") +
    theme_void() + theme(plot.margin=margin(5,3,4,3))
  ((pF + pN + plot_layout(widths=c(1, 0.36))) / strip) + plot_layout(heights=c(1, 0.17))
}

ggsave(file.path(OUT_DIR,"Fig6a_TF_family_forest.pdf"), build(110), width=110, height=78, units="mm", device=cairo_pdf)
ggsave(file.path(OUT_DIR,"Fig6a_TF_family_forest.png"), build(110), width=110, height=78, units="mm", dpi=400, bg="white")
ggsave(file.path(OUT_DIR,"Fig6a_TF_family_forest_large.png"), build(176), width=176, height=125, units="mm", dpi=300, bg="white")
cat("wrote panelA forest (pdf, png, large png)\n\n")
print(d %>% arrange(desc(row), genotype) %>%
      mutate(across(c(OR,lo,hi),~round(.,2)), q=signif(q,3)) %>%
      select(family,genotype,lab,OR,lo,hi,q,sig), row.names=FALSE)
nz <- d[!d$zero, ]
stopifnot(all(nz$lo >= XL[1]), all(nz$hi <= XL[2]))   # every plotted CI is fully inside the panel
stopifnot("an OR of 0 can never be significant" = !any(d$zero & d$sig),
          "DOF must be in the plot" = "Dof" %in% d$family)
cat(sprintf("\nall CIs inside x range %.2f-%.1f; %d zero-OR arrow(s); filled points: %s\n", XL[1], XL[2], sum(d$zero),
    paste(paste(d$family[d$sig], d$genotype[d$sig]), collapse=", ")))
