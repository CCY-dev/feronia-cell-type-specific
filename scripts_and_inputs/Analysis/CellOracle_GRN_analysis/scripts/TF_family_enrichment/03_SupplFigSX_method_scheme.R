#!/usr/bin/env Rscript
## 03 — Supplementary Fig SX: graphical method scheme (base version; the published figure was simplified by hand).
## Network glyphs are schematic; the Δrank density, thresholds and all counts are computed. Needs 01.
# 20260924 — GRAPHICAL method scheme: real network pictures with arrows showing what is extracted
# at each step, ending at the TF-family enrichment figure.
#
# What is DRAWN vs what is REAL:
#   * the mini-network glyphs (stages 1-3) are SCHEMATIC — deterministic (seeded) toy graphs whose
#     only job is to show the operation being performed. They are not the real 20k-link networks,
#     which cannot be drawn legibly at thumbnail size.
#   * the Delta-rank distribution (stage 4) and every COUNT in the figure are REAL, computed from the
#     source tables at build time.
suppressMessages({library(ggplot2);library(dplyr)})
source("00_paths.R")
CW <- "#2b6cb0"; CF <- "#dd6b20"; CN <- "#16324f"; CH <- "#c0392b"; CG <- "grey55"
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)

## ---------------- real numbers ----------------
fam <- rd(PLANTTFDB); famv <- setNames(fam$Family, toupper(fam$Gene_ID))
w   <- rd(WIDE_TABLE)
w$gene <- toupper(w$gene)
o   <- w[w$metric=="degree_centrality_out" & w$Is_TF=="Yes", ]
RDC <- "fer4_untreated-WT_untreated: rank_diff"
ELIG<- unique(o$gene); NCT <- length(unique(o$cell_type))
gWT <- unique(o$gene[o$Label=="WT_specific"]); gFER <- unique(o$gene[o$Label=="fer4_specific"])
DOFg<- ELIG[!is.na(famv[ELIG]) & famv[ELIG]=="Dof"]
res <- read.csv(file.path(OUT_DIR,"TF_family_enrichment_all_tests.csv"), stringsAsFactors=FALSE) %>% filter(eligible >= 10)
dr  <- res %>% filter(genotype=="fer-4", family=="Dof")
dw  <- res %>% filter(genotype=="WT",    family=="Dof")
## 2026-10-02 (PI): enrichment background = the 1,292 TFs that could enter the ranking
BGT <- read.csv(file.path(OUT_DIR,"enrichment_background_TFs.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(BGT) == 1292)
## 2026-10-02: percentiles over ALL out-degree rows (TFs + the 10 non-TF gene rows), exactly as the pipeline
## computed them. TF rows only gave 90th = 132.7, printed as "133"; the applied cut-off is 132.0 (same calls,
## since Delta-rank is an integer).
rdv <- w[w$metric=="degree_centrality_out", RDC]; rdv <- rdv[is.finite(rdv)]
stopifnot(length(rdv) == 5534)
q10 <- quantile(rdv,.10); q90 <- quantile(rdv,.90)
cat(sprintf("cell types %d | eligible %d | fer4 %d | WT %d | DOF %d | q10 %.0f q90 %.0f\n",
    NCT, length(ELIG), length(gFER), length(gWT), length(DOFg), q10, q90))

## ---------------- helper: a schematic mini-network ----------------
mininet <- function(cx, cy, r, nn=9, ne=13, seed=1, hub=NA) {
  set.seed(seed)
  a <- sort(runif(nn, 0, 2*pi)); rr <- r*sqrt(runif(nn, .18, 1))
  nd <- data.frame(id=1:nn, x=cx+rr*cos(a), y=cy+rr*sin(a))
  ed <- data.frame(from=sample(nn, ne, TRUE), to=sample(nn, ne, TRUE))
  ed <- ed[ed$from != ed$to, ]
  if (!is.na(hub)) ed <- rbind(data.frame(from=hub, to=setdiff(1:nn, hub)[1:5]), ed)
  ed$x <- nd$x[ed$from]; ed$y <- nd$y[ed$from]
  ed$xe<- nd$x[ed$to];   ed$ye<- nd$y[ed$to]
  ed$ishub <- !is.na(hub) & ed$from == hub
  nd$ishub <- !is.na(hub) & nd$id == hub
  list(n=nd, e=ed)
}
NODES <- list(); EDGES <- list(); k <- 0
addnet <- function(g, tag) { k <<- k+1; g$n$tag <- tag; g$e$tag <- tag
  NODES[[k]] <<- g$n; EDGES[[k]] <<- g$e }

## ---- stage 1: a grid of per-context networks, WT row and fer-4 row ----
for (j in 0:2) { addnet(mininet(6+j*9.5, 47.5, 3.1, seed=10+j), "WT")
                 addnet(mininet(6+j*9.5, 38.5, 3.1, seed=20+j), "FER") }
## ---- stage 2: one enlarged network with a hub highlighted ----
addnet(mininet(44, 43, 6.2, nn=12, ne=16, seed=7, hub=1), "HUB")

nodes <- bind_rows(NODES); edges <- bind_rows(EDGES)

## ---- stage 3: two rank ladders, one gene tracked between them ----
set.seed(3); NL <- 9
lad <- bind_rows(
  data.frame(g="WT",   pos=1:NL, x=66, y=seq(49.5, 37.5, length.out=NL)),
  data.frame(g="fer4", pos=1:NL, x=74, y=seq(49.5, 37.5, length.out=NL)))
trackW <- 7; trackF <- 2                       # the tracked gene climbs from rank 7 to rank 2
lad$hit <- (lad$g=="WT" & lad$pos==trackW) | (lad$g=="fer4" & lad$pos==trackF)
# the top-100 screen: only the top of each ladder is kept. The tracked gene is BELOW the cut in WT
# but ABOVE it in fer-4 -- exactly the case the UNION rule keeps, which is the point of the panel.
KEEP <- 4
lad$kept <- lad$pos <= KEEP
YCUT <- mean(lad$y[lad$g=="WT" & lad$pos %in% c(KEEP, KEEP+1)])
trk <- data.frame(x=66, y=lad$y[lad$g=="WT" & lad$pos==trackW],
                  xe=74, ye=lad$y[lad$g=="fer4" & lad$pos==trackF])

## ---- stage 4: the REAL Delta-rank distribution ----
# Delta-rank is extremely peaked at 0, so the full support renders as a spike with the tails
# invisible. Clip to the central 99% for DISPLAY; the 10th/90th percentiles are still the real ones.
# Delta-rank is tight: IQR is only -30..33 while the 0.5/99.5 pct reach -609/+687, so a quantile
# window leaves the 10th/90th thresholds crammed together in the middle. A symmetric +/-400 window
# covers ~95% of contexts and puts the thresholds where they can be read.
WIN <- c(-400, 400)
dn <- density(rdv, n=1024, from=WIN[1], to=WIN[2])
dd <- data.frame(x=dn$x, y=dn$y); SC <- 9.0 / max(dd$y)
dd$yy <- 6 + dd$y*SC
mapx <- function(v) 12 + 22*(v - WIN[1])/(WIN[2] - WIN[1])
dd$x2 <- mapx(dd$x); xq10 <- mapx(q10); xq90 <- mapx(q90)

## ---- stage 5: the forest-plot glyph ----
## 2026-10-02: the glyph is now the REAL fer-4 column (top 7 families by OR, OR > 0, with real CIs and
## real BH fill), so it can no longer imply a single significant family.
NF <- 7
fr <- res %>% filter(genotype=="fer-4", OR > 0) %>% arrange(desc(OR)) %>% head(NF)
stopifnot(fr$family[1] == "Dof")
fo <- data.frame(y=seq(20.5, 13.5, length.out=NF), or=fr$OR, sig=fr$q < 0.05)
fo$x <- 72 + 3.0*log2(fo$or); fo$lo <- 72 + 3.0*log2(fr$lo); fo$hi <- 72 + 3.0*log2(fr$hi)
stopifnot(all(fo$lo > 60), all(fo$hi < 92))   # glyph stays inside the stage-5 area

# 2026-09-25 (user): the mutant allele is ALWAYS italic. `pm=TRUE` renders the label as plotmath so
# only "fer-4" is italicised. ⚠ plotmath ignores `fontface`, so bold must be written INTO the
# expression (bold(...) / bolditalic(...)) or those labels render lighter than their neighbours.
## 2026-09-25: every text call goes through lbl(), so recording here captures ALL labels for the
## editable PPTX export with no risk of the two versions drifting apart.
.LAB <- list()
lbl <- function(x, y, t, s=2.5, col=CN, f="bold", h=0.5, pm=FALSE) {
  .LAB[[length(.LAB)+1]] <<- data.frame(x=x, y=y, text=t, size=s, colour=col,
                                        face=f, hjust=h, plotmath=pm, stringsAsFactors=FALSE)
  # ⚠ braces required: a bare `else` on its own line does not parse in R
  if (pm) annotate("text", x=x, y=y, label=t, size=s, colour=col, hjust=h, parse=TRUE)
  else    annotate("text", x=x, y=y, label=t, size=s, colour=col, fontface=f, hjust=h)
}
ar <- function(x, xe, y, ye=y) annotate("segment", x=x, xend=xe, y=y, yend=ye,
       arrow=arrow(length=unit(1.7,"mm"), type="closed"), colour=CG, linewidth=.42)

p <- ggplot() +
  ## stage 1 + 2 networks
  geom_segment(data=edges %>% filter(!ishub), aes(x,y,xend=xe,yend=ye), colour="grey72", linewidth=.22) +
  geom_segment(data=edges %>% filter(ishub), aes(x,y,xend=xe,yend=ye), colour=CH, linewidth=.5,
               arrow=arrow(length=unit(1.0,"mm"), type="closed")) +
  geom_point(data=nodes %>% filter(!ishub),
             aes(x, y, fill=ifelse(tag=="WT", CW, ifelse(tag=="FER", CF, CN))),
             shape=21, size=1.5, colour="white", stroke=.25) +
  geom_point(data=nodes %>% filter(ishub), aes(x,y), shape=21, size=3.6, fill=CH, colour="white", stroke=.4) +
  scale_fill_identity() +
  lbl(1.2, 47.5, "WT", 2.5, CW, "bold", 0.5) + lbl(1.2, 38.5, 'bolditalic("fer-4")', 2.5, CF, "bold", 0.5, pm=TRUE) +
  lbl(6, 52.0, "cell type 1", 2.1, CG, "plain") +
  lbl(15.5, 52.0, "cell type 2", 2.1, CG, "plain") +
  lbl(25, 52.0, "cell type 3", 2.1, CG, "plain") +
  lbl(31.5, 43, "…", 4, CG, "plain") +
  lbl(15.5, 31.4, sprintf("%d cell types × 2 genotypes = %d networks", NCT, 2*NCT), 2.2, CG, "plain") +
  lbl(15.5, 55.5, "1 · infer a GRN in every context", 2.7, CN) +
  lbl(44, 55.5, "2 · out-degree centrality", 2.7, CN) +
  lbl(44, 33.0, "count each gene's OUTGOING links", 2.2, CH, "plain") +
  lbl(44, 30.6, "(how many targets it regulates)", 2.2, CG, "plain") +
  ar(34, 37, 43) +
  ## stage 3 ladders
  ar(51.5, 60.5, 43) +
  lbl(70, 55.5, "3 · rank, screen, then Δrank", 2.7, CN) +
  lbl(66, 51.6, "WT", 2.4, CW) + lbl(74, 51.6, 'bolditalic("fer-4")', 2.4, CF, pm=TRUE) +
  geom_point(data=lad, aes(x, y,
                           fill=ifelse(!kept, "grey80", ifelse(g=="WT", CW, CF)),
                           size=ifelse(hit, 2.4, 1.3)),
             shape=21, colour="white", stroke=.3) +
  scale_size_identity() +
  annotate("segment", x=64.2, xend=75.8, y=YCUT, yend=YCUT, colour="grey45", linetype=2, linewidth=.34) +
  lbl(76.4, YCUT, "top 100", 2.0, CG, "plain", 0) +
  geom_segment(data=trk, aes(x,y,xend=xe,yend=ye), colour=CH, linewidth=.5, linetype=2) +
  lbl(70, 36.2, "keep only the top 100 of each network;", 2.1, CG, "plain") +
  lbl(70, 34.0, sprintf("the UNION enters the ranking = %d TFs", length(ELIG)), 2.15, CN) +
  lbl(70, 31.8, 'paste(Delta, "rank = rank(", italic("fer-4"), ") \u2212 rank(WT)")', 2.05, CH, pm=TRUE) +
  ## connector into row 2 -- plain line, no arrowheads (user 20260924)
  annotate("segment", x=70, xend=70, y=29.9, yend=27.9, colour=CG, linewidth=.42) +
  annotate("segment", x=70, xend=23, y=27.9, yend=27.9, colour=CG, linewidth=.42) +
  annotate("segment", x=23, xend=23, y=27.9, yend=26.2, colour=CG, linewidth=.42) +
  ## stage 4: real Delta-rank density with the two tails
  lbl(23, 24.2, "4 · threshold the Δrank distribution", 2.7, CN) +
  geom_ribbon(data=dd %>% filter(x <= q10), aes(x=x2, ymin=6, ymax=yy), fill=CF, alpha=.85) +
  geom_ribbon(data=dd %>% filter(x >= q90), aes(x=x2, ymin=6, ymax=yy), fill=CW, alpha=.85) +
  geom_line(data=dd, aes(x2, yy), colour=CN, linewidth=.42) +
  annotate("segment", x=12, xend=34, y=6, yend=6, colour="grey40", linewidth=.3) +
  annotate("segment", x=xq10, xend=xq10, y=5.2, yend=6.6, colour=CF, linewidth=.4) +
  annotate("segment", x=xq90, xend=xq90, y=5.2, yend=6.6, colour=CW, linewidth=.4) +
  lbl(xq10 - 1.6, 4.3, sprintf("10th pct\n%.0f", q10), 1.95, CF, "plain") +
  lbl(xq90 + 1.6, 4.3, sprintf("90th pct\n%.0f", q90), 1.95, CW, "plain") +
  lbl(23, 2.0, sprintf("\u0394rank   (axis clipped to \u00b1400; %.0f%% of contexts)",
                       100*mean(abs(rdv)<=400)), 1.95, CG, "plain") +
  lbl(14.6, 14.4, 'paste("more central in ", italic("fer-4"))', 1.95, CG, pm=TRUE) +
  lbl(14.6, 12.2, sprintf('paste(bold("%d "), bolditalic("fer-4"), bold("-specific"))', length(gFER)), 2.2, CF, pm=TRUE) +
  lbl(31.4, 14.4, "more central in WT", 1.95, CG, "plain") +
  lbl(31.4, 12.2, sprintf("%d WT-specific", length(gWT)), 2.2, CW) +
  ar(37, 51, 17.5) +
  ## stage 5 (was 6): forest glyph
  lbl(72, 24.2, "5 · family enrichment", 2.7, CH) +
  annotate("segment", x=72, xend=72, y=12.6, yend=21.6, colour="grey60", linetype=2, linewidth=.3) +
  geom_segment(data=fo, aes(x=lo, xend=hi, y=y, yend=y), colour="grey45", linewidth=.34) +
  geom_point(data=fo, aes(x, y, fill=ifelse(sig, CH, "white")), shape=21, size=2.1,
             colour="grey30", stroke=.35) +
  lbl(70.4, 20.5, "DOF", 2.1, CH, "bold", 1) +
  lbl(72, 10.6, sprintf('paste(bold("DOF:  OR %.2f (WT),  %.2f ("), bolditalic("fer-4"), bold(")"))', dw$OR, dr$OR), 2.2, CH, pm=TRUE) +
  lbl(72, 8.2, sprintf("Fisher's exact per family (BH) vs %s TFs scored in both genotypes",
                       format(nrow(BGT), big.mark=",")), 2.0, CG, "plain") +
  coord_cartesian(xlim=c(-2, 98), ylim=c(3, 57), clip="off") +
  # NOTE: the "glyphs are schematic" disclaimer now lives ONLY in the report caption
  # (20260923_DOF_enrichment_report.R) -- keep it there if this figure is reused elsewhere.
  labs(title="From per cell type GRN to TF-family enrichment") +
  theme_void(base_size=8) +
  theme(plot.title=element_text(size=9, face="bold", hjust=0),
        legend.position="none", plot.margin=margin(6,5,4,5))

ggsave(file.path(OUT_DIR,"SupplFigSX_method_scheme.pdf"), p, width=180, height=112, units="mm", device=cairo_pdf)
ggsave(file.path(OUT_DIR,"SupplFigSX_method_scheme.png"), p, width=180, height=112, units="mm", dpi=400, bg="white")
write.csv(data.frame(n_ct=NCT, n_ctx=2*NCT, n_elig=length(ELIG), n_fer=length(gFER), n_wt=length(gWT),
                     n_dof=length(DOFg), q10=unname(q10), q90=unname(q90), dof_OR=dr$OR, dof_q=dr$q),
          file.path(OUT_DIR,"SupplFigSX_method_scheme_stats.csv"), row.names=FALSE)
cat("wrote 20260924_method_scheme_graphical.{pdf,png} + stats\n")
