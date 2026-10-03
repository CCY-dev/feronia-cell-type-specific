#!/usr/bin/env Rscript
## 05 — Fig X2: bipartite DOF → auxin-gene network for the seven DOFs named in the text. Needs 04 (edge tables).
# 20260925 — Panel D, redesigned: YUC5 as the anchor, the auxin pathway as context.
#
# ⚠ THIS PANEL SITS BEFORE 6D/6E, so it is a PREDICTION figure: nothing may reference the LUC or EMSA
#   results shown later. No experimental-status encoding.
# ⚠ READABILITY AT SMALL SIZE (PI): 31 nodes vs 458 edges — edge colour is the first thing to die when
#   the figure shrinks, so the genotype signal lives in the NODE FILL and edges are plain structure.
#   One variable per channel:
#       DOF SIZE   = out-degree (auxin targets reached) -- the quantity the analysis rests on
#       DOF FILL   = its genotype-specific label        -- ties straight to panel A
#       target FILL= auxin pathway class                -- intrinsic, not a result
#       YUC5       = darkened EDGES + bold label        -- the anchor (no size/outline)
# ⚠ ggplot permits ONE fill scale, and two different node types need different palettes, so every
#   node carries a precomputed hex under scale_fill_identity() and the legend is drawn by hand.
# ⚠ YUC5 is the best-connected BIOSYNTHESIS gene (21 vs YUC9 10, TAA1 3) but NOT the best-connected
#   target overall (IAA13/IAA2/IAA4 have 28; it is 10th of 32). The title must not imply otherwise.
suppressMessages({library(ggplot2); library(dplyr)})
source("00_paths.R")
LAB_COL <- c(`fer-4-specific`="#dd6b20", `WT-specific`="#2b6cb0", both="#6b7f95", `not nominated`="#c8ced6")
CLS     <- c(`YUC/TAA`="#e8b30b", GH3="#d6529b", `Aux/IAA`="#5aa7dd", SAUR="#1f9d74")
has_repel <- requireNamespace("ggrepel", quietly=TRUE)

e  <- read.csv(file.path(OUT_DIR,"DOF_auxin_edges_all_contexts.csv"), stringsAsFactors=FALSE)
pr <- read.csv(file.path(OUT_DIR,"DOF_auxin_pair_states.csv"), stringsAsFactors=FALSE)
g  <- function(v,k) { x<-setNames(e[[v]],e[[k]]); x[!duplicated(names(x))] }
pr$ssym <- g("source_sym","source")[pr$source]
pr$tsym <- g("target_sym","target")[pr$target]
pr$class<- g("target_class","target")[pr$target]

w <- read.delim(WIDE_TABLE,
                check.names=FALSE, stringsAsFactors=FALSE)
w$gene <- toupper(w$gene); o <- w[w$metric=="degree_centrality_out" & w$Is_TF=="Yes", ]
gWT <- unique(o$gene[o$Label=="WT_specific"]); gFER <- unique(o$gene[o$Label=="fer4_specific"])

## ⚠⚠ 2026-09-25 (user): the panel now shows ONLY THE DOFs NAMED IN THE MANUSCRIPT TEXT, so the
## figure and the prose reference the same genes. Keyed by AGI, NOT by symbol: the paper uses COMMON
## names and features.tsv uses SYSTEMATIC ones, so a symbol join silently drops COG1 and OBP1.
NAMED <- c(AT4G24060="DOF4.6", AT3G61850="DAG1",  AT3G50410="OBP1", AT1G07640="OBP2",
           AT3G55370="OBP3",   AT2G34140="CDF4",  AT1G29160="COG1")
ALIAS <- c(AT3G50410="OBP1", AT1G29160="COG1")   # the two whose systematic name differs in our tables
lab_of <- function(v) ifelse(v %in% gWT & v %in% gFER, "both",
                      ifelse(v %in% gFER, "fer-4-specific",
                      ifelse(v %in% gWT,  "WT-specific", "not nominated")))
pr$lab <- lab_of(pr$source)
pr_all <- pr
stopifnot("not every named DOF is in the auxin edge set" = all(names(NAMED) %in% pr$source))
pr <- pr %>% filter(source %in% names(NAMED))
NDROP  <- n_distinct(pr_all$source) - n_distinct(pr$source)
NEDROP <- nrow(pr_all) - nrow(pr)
stopifnot("a named DOF was not nominated in either genotype" = !any(pr$lab == "not nominated"),
          n_distinct(pr$source) == length(NAMED))
cat(sprintf("restricted to the %d DOFs named in the text (dropped %d other DOFs, %d of %d edges)\n",
            n_distinct(pr$source), NDROP, NEDROP, nrow(pr_all)))

ind  <- pr %>% count(tsym, class, name="indeg")
keep <- c(ind$tsym[ind$class %in% c("YUC/TAA","GH3")],
          (ind %>% filter(class=="Aux/IAA") %>% arrange(desc(indeg)) %>% head(5))$tsym,
          (ind %>% filter(class=="SAUR")    %>% arrange(desc(indeg)) %>% head(3))$tsym)
k <- pr %>% filter(tsym %in% keep)
cat(sprintf("kept %d of %d targets | %d of %d edges (%.0f%%)\n",
    length(keep), nrow(ind), nrow(k), nrow(pr), 100*nrow(k)/nrow(pr)))

ord <- c("YUC/TAA","GH3","Aux/IAA","SAUR")                 # pathway order, top to bottom
tars <- k %>% count(tsym, class, name="indeg") %>%
        mutate(class=factor(class, levels=ord)) %>% arrange(class, desc(indeg)) %>%
        mutate(y=rev(seq_along(tsym)), x=1, anchor=tsym=="YUC5", fillc=CLS[as.character(class)])
dofs <- k %>% count(source, ssym, name="outdeg") %>% arrange(outdeg) %>%
        mutate(y=seq_along(source), x=0, lab=lab_of(source), fillc=LAB_COL[lab],
               # ⚠ print BOTH names for the two aliases, or the reader cannot match this panel to the
               # LUC/EMSA panels, which use the common name only.
               disp=ifelse(source %in% names(ALIAS),
                           sprintf("%s (%s)", ssym, unname(ALIAS[source])), ssym))
stopifnot("a not-nominated DOF survived into the figure" = !any(dofs$lab == "not nominated"))
sc <- function(v,n) if (n<=1) rep(.5,length(v)) else (v-1)/(n-1)
dofs$yy <- sc(dofs$y, nrow(dofs)); tars$yy <- sc(tars$y, nrow(tars))
ed <- k %>% left_join(dofs[,c("source","yy")], by="source") %>% rename(y0=yy) %>%
            left_join(tars[,c("tsym","yy")],  by="tsym")   %>% rename(y1=yy) %>%
            mutate(anch = tsym=="YUC5")
cat(sprintf("DOF labels: %s\n", paste(names(table(dofs$lab)), table(dofs$lab), sep="=", collapse="  ")))
cat(sprintf("edges into YUC5: %d\n", sum(ed$anch)))

## hand-drawn legend (identity fill leaves no automatic one)
## ⚠ 2026-09-25 (user): the legend moved from UNDER the tiers to the RIGHT, one item per row, and the
## two tiers were pulled closer together to pay for it — the figure keeps its 168 mm width, the space
## just moves from the (empty) middle to the legend column.
TARX <- 0.52                  # target tier x, was 1.0
LEGX <- 0.86                  # legend column: dot at LEGX, its label at LEGX + LEGPAD
LEGPAD <- 0.030
tars$x <- TARX
ed$xend <- TARX

LABP <- LAB_COL[names(LAB_COL) %in% dofs$lab]      # only the labels still on the figure
sk   <- data.frame(o=unique(as.integer(round(quantile(dofs$outdeg, c(0,.5,1))))))
## ⚠ row positions are COMPUTED from the section sizes, so adding or dropping a legend entry can
## never silently overlap the next section.
ROWH <- 0.070; SECG <- 0.052
n1 <- length(LABP); n2 <- length(ord); n3 <- nrow(sk)
yh1 <- 1.03;                     ye1 <- yh1 - ROWH - (0:(n1-1))*ROWH
yh2 <- min(ye1) - ROWH - SECG;   ye2 <- yh2 - ROWH - (0:(n2-1))*ROWH
yh3 <- min(ye2) - ROWH - SECG;   ye3 <- yh3 - ROWH - (0:(n3-1))*ROWH
stopifnot("the legend is taller than the panel" = min(ye3) > -0.02)

# ⚠ capitalise the DISPLAY only -- LAB_COL's names are matched by lab_of(), so renaming the level
# itself would silently break the fill lookup. fer-4 keeps its lower case: it is an allele name.
lg <- rbind(
  data.frame(x=LEGX, y=ye1, f=unname(LABP),
             t=ifelse(names(LABP)=="both", "Both", names(LABP))),
  data.frame(x=LEGX, y=ye2, f=unname(CLS[ord]), t=ord))
# ⚠ the allele is ALWAYS italic (user, 2026-09-25). plotmath ignores plain strings' fontface, so the
# one entry naming it is drawn parsed and every other entry stays plain text.
lg$e <- ifelse(lg$t == "fer-4-specific", 'paste(italic("fer-4"), "-specific")', NA_character_)
sk$x <- LEGX; sk$y <- ye3

p <- ggplot() +
  geom_segment(data=ed %>% filter(!anch), aes(x=0, y=y0, xend=xend, yend=y1),
               # ⚠ tuning history (user 2026-09-25): grey72/.20 was set for the 221-edge version and
               # read as blank space after the trim; grey45/.55 was then too heavy and competed with
               # the gold YUC5 edges. This is the midpoint.
               colour="grey58", alpha=.38, linewidth=.19) +
  geom_segment(data=ed %>% filter(anch), aes(x=0, y=y0, xend=xend, yend=y1),
               colour="#8a6d1f", alpha=.70, linewidth=.34) +
  geom_point(data=dofs, aes(x, yy, size=outdeg, fill=fillc), shape=21, colour="white", stroke=.45) +
  # ⚠ YUC5 is NOT enlarged or ringed: the darkened edges plus the bold label are the only emphasis.
  geom_point(data=tars, aes(x, yy, fill=fillc), shape=21, size=2.7, colour="white", stroke=.4) +
  geom_text(data=dofs, aes(x=-0.035, y=yy, label=disp), hjust=1, size=2.6, colour="grey20") +
  geom_text(data=tars %>% filter(!anchor), aes(x=TARX+0.032, y=yy, label=tsym), hjust=0, size=2.0, colour="grey25") +
  geom_text(data=tars %>% filter(anchor),  aes(x=TARX+0.040, y=yy, label=tsym), hjust=0, size=2.6,
            fontface="bold", colour="#5a4a08") +
  annotate("text", x=0,    y=1.07, label="DOF regulators", size=2.5, fontface="bold", colour="grey20") +
  annotate("text", x=TARX, y=1.07, label="Auxin genes",    size=2.5, fontface="bold", colour="grey20") +
  ## ---- legend column, one item per row ----
  annotate("text", x=LEGX-0.022, y=yh1, label="DOF:",   hjust=0, size=2.2, fontface="bold", colour="grey20") +
  annotate("text", x=LEGX-0.022, y=yh2, label="Target:", hjust=0, size=2.2, fontface="bold", colour="grey20") +
  annotate("text", x=LEGX-0.022, y=yh3, label="Auxin targets reached", hjust=0, size=2.2,
           fontface="bold", colour="grey20") +
  geom_point(data=lg, aes(x, y, fill=f), shape=21, size=2.6, colour="white", stroke=.4) +
  geom_text(data=lg[is.na(lg$e),], aes(x+LEGPAD, y, label=t), hjust=0, size=2.0, colour="grey20") +
  geom_text(data=lg[!is.na(lg$e),], aes(x+LEGPAD, y, label=e), hjust=0, size=2.0, colour="grey20", parse=TRUE) +
  geom_point(data=sk, aes(x, y, size=o), shape=21, fill="#b9c2cc", colour="white", stroke=.4) +
  geom_text(data=sk, aes(x+LEGPAD, y, label=o), hjust=0, size=2.0, colour="grey35") +
  scale_fill_identity() +
  scale_size_continuous(range=c(1.6,5.0), breaks=c(5,10,15), name="auxin targets") +
  scale_x_continuous(limits=c(-0.54, 1.10)) + scale_y_continuous(limits=c(-0.03, 1.12)) +
  labs(title="DOF regulators are predicted to target the auxin pathway",
       subtitle=sprintf("the %d DOFs named in the text \u2192 %d auxin genes \u00b7 node size = auxin targets reached\nYUC5 highlighted: the best-connected auxin BIOSYNTHESIS gene, reached by %d of the %d",
                        nrow(dofs), nrow(tars), tars$indeg[tars$anchor], nrow(dofs))) +
  theme_void(base_size=8) +
  theme(plot.title=element_text(size=9, face="bold"),
        plot.subtitle=element_text(size=6.4, colour="grey35", lineheight=1.2),
        legend.position="none", plot.margin=margin(5,4,4,4))

ggsave(file.path(OUT_DIR,"FigX2_DOF_auxin_network.pdf"), p, width=168, height=96, units="mm", device=cairo_pdf)
ggsave(file.path(OUT_DIR,"FigX2_DOF_auxin_network.png"), p, width=168, height=96, units="mm", dpi=400, bg="white")
# stats for the report, so no caption number is ever typed by hand
lc <- table(factor(dofs$lab, levels=names(LAB_COL)))
write.csv(data.frame(n_dof=nrow(dofs), n_tar_kept=nrow(tars), n_tar_all=nrow(ind),
                     n_edge_kept=nrow(k), n_edge_all=nrow(pr), n_edge_everything=nrow(pr_all),
                     n_dof_dropped=NDROP, n_edge_dropped=NEDROP,
                     named=paste(dofs$disp[order(-dofs$outdeg)], collapse='; '),
                     yuc5_indeg=tars$indeg[tars$anchor], yuc9_indeg=tars$indeg[tars$tsym=="YUC9"],
                     taa1_indeg=tars$indeg[tars$tsym=="TAA1"],
                     max_target_indeg=max(tars$indeg),
                     # rank of YUC5 among ALL auxin targets (before the target thinning) -- guards the
                     # "best BIOSYNTHESIS gene, not best overall" caveat against going stale
                     yuc5_rank_overall=sum(ind$indeg > ind$indeg[ind$tsym=="YUC5"]) + 1,
                     n_fer=lc[["fer-4-specific"]], n_wt=lc[["WT-specific"]],
                     n_both=lc[["both"]], n_none=lc[["not nominated"]],
                     min_outdeg=min(dofs$outdeg), max_outdeg=max(dofs$outdeg), med_outdeg=median(dofs$outdeg)),
          file.path(OUT_DIR,"FigX2_DOF_auxin_network_stats.csv"), row.names=FALSE)
cat("wrote 20260925_panelD_YUC5_anchor.{pdf,png}\n")
print(as.data.frame(tars %>% select(tsym, class, indeg)), row.names=FALSE)
