#!/usr/bin/env Rscript
## 05 — Fig. 6b: bipartite network of the 10 DOFs tested in the pYUC5:LUC assay and their predicted auxin-gene
## targets. Needs 04 (edge tables). Reproduces the PI's Fig. 6b layout (2026-10-19) from the published edge table.
#
# Selection rules (all derived from the data, asserted below):
#   DOFs    = the 10 DOFs tested in the pYUC5:LUC assay ("selected 10 out of the 26 DOFs based on their expression
#             levels and cell types"), keyed by AGI. The paper uses COMMON names and the gene annotation uses
#             SYSTEMATIC ones (DOF6 = DOF3.2, OBP1 = DOF3.4, COG1 = DOF1.5), so never join on a symbol.
#   groups  = genotype-specific nomination (Δrank of out-degree), as in Fig. 6a: WT only / fer-4 only / Both.
#   targets = every YUC/TAA1 and GH3 gene, the 6 Aux/IAA and 3 SAUR genes reached by the most of the 10 DOFs.
#             Ties at either cut-off would make the choice arbitrary, so they stop the script.
#   order   = targets by pathway class, then by number of DOFs reaching them (desc), then name; DOFs within a
#             group with a YUC5 link first, then by name.
# Encoding: DOF node size = number of the displayed targets it reaches; YUC5 edges dark, all others light grey.
# An edge = the link was recovered in at least one cell type x genotype network (union over contexts).
suppressMessages({library(ggplot2); library(dplyr)})
source("00_paths.R")
CLS <- c(`YUC/TAA1`="#e1a62b", GH3="#d07fa8", `Aux/IAA`="#3a78b5", SAUR="#3aa17a")

e  <- read.csv(file.path(OUT_DIR,"DOF_auxin_edges_all_contexts.csv"), stringsAsFactors=FALSE)
pr <- read.csv(file.path(OUT_DIR,"DOF_auxin_pair_states.csv"), stringsAsFactors=FALSE)
sy <- e %>% distinct(target, tsym=target_sym, class=target_class) %>%
      mutate(class=ifelse(class=="YUC/TAA", "YUC/TAA1", class))
stopifnot("a target has more than one symbol/class" = !anyDuplicated(sy$target))

w <- read.delim(WIDE_TABLE, check.names=FALSE, stringsAsFactors=FALSE)
w$gene <- toupper(w$gene); o <- w[w$metric=="degree_centrality_out" & w$Is_TF=="Yes", ]
gWT <- unique(o$gene[o$Label=="WT_specific"]); gFER <- unique(o$gene[o$Label=="fer4_specific"])
grp_of <- function(v) ifelse(v %in% gWT & v %in% gFER, "Both",
                      ifelse(v %in% gFER, "fer-4 only", ifelse(v %in% gWT, "WT only", NA)))

## the 10 LUC-tested DOFs, with the names used in the manuscript
DOF10 <- c(AT2G34140="CDF4", AT3G45610="DOF6",   AT3G50410="OBP1",   AT1G29160="COG1", AT3G61850="DAG1",
           AT4G00940="DOF4.1", AT4G24060="DOF4.6", AT1G07640="OBP2", AT3G55370="OBP3", AT2G28810="DOF2.2")
k <- pr %>% filter(source %in% names(DOF10)) %>% left_join(sy, by="target")
stopifnot("not every tested DOF has an auxin edge" = setequal(unique(k$source), names(DOF10)),
          "a tested DOF was not nominated in either genotype" = !anyNA(grp_of(names(DOF10))))

ind <- k %>% count(tsym, class, name="indeg") %>% arrange(desc(indeg), tsym)
topn <- function(cl, n) {
  x <- ind %>% filter(class == cl)
  stopifnot("tie at a target cut-off" = nrow(x) <= n || x$indeg[n] > x$indeg[n + 1])
  head(x$tsym, n)
}
keep <- c(ind$tsym[ind$class %in% c("YUC/TAA1","GH3")], topn("Aux/IAA", 6), topn("SAUR", 3))
ks <- k %>% filter(tsym %in% keep)

## ---- layout: rows from the top, with a gap between groups; both tiers span the same height ----
GAP <- 0.9
ord <- c("YUC/TAA1","GH3","Aux/IAA","SAUR")
tars <- ind %>% filter(tsym %in% keep) %>% mutate(class=factor(class, levels=ord)) %>%
        arrange(class, desc(indeg), tsym) %>%
        mutate(r = row_number() - 1 + GAP * (as.integer(class) - 1))
GRP <- c("WT only","fer-4 only","Both")
yuc5_src <- unique(ks$source[ks$tsym == "YUC5"])
dofs <- ks %>% count(source, name="size") %>%
        mutate(name=unname(DOF10[source]), grp=factor(grp_of(source), levels=GRP), yuc5=source %in% yuc5_src) %>%
        arrange(grp, desc(yuc5), name) %>%
        mutate(r = row_number() - 1 + GAP * (as.integer(grp) - 1))
sc <- function(r) 1 - r / max(r)
tars$y <- sc(tars$r); dofs$y <- sc(dofs$r)

TX <- 1                                         # target tier x (DOF tier at 0)
ed <- ks %>% left_join(dofs %>% select(source, y0=y), by="source") %>%
             left_join(tars %>% select(tsym, y1=y), by="tsym") %>% mutate(anch = tsym == "YUC5")

## group brackets: a vertical bar beside each group, label to its left (fer-4 italic via plotmath)
PADB <- 0.025
gb <- dofs %>% group_by(grp) %>% summarise(lo=min(y) - PADB, hi=max(y) + PADB, mid=mean(y), .groups="drop") %>%
      mutate(lab=c(`WT only`='"WT only"', `fer-4 only`='italic("fer-4")~"only"', Both='"Both"')[as.character(grp)])
cb <- tars %>% group_by(class) %>% summarise(lo=min(y) - PADB, hi=max(y) + PADB, mid=mean(y), .groups="drop") %>%
      mutate(col=CLS[as.character(class)])

BARX <- -0.70; CBARX <- 1.36
p <- ggplot() +
  geom_segment(data=ed %>% filter(!anch), aes(x=0, y=y0, xend=TX, yend=y1), colour="grey78", linewidth=.22) +
  geom_segment(data=ed %>% filter(anch),  aes(x=0, y=y0, xend=TX, yend=y1), colour="grey25", linewidth=.45) +
  geom_point(data=dofs, aes(0, y, size=size), shape=21, fill="grey40", colour="white", stroke=.4) +
  geom_point(data=tars, aes(TX, y, fill=I(CLS[as.character(class)])), shape=21, size=1.9, colour="white", stroke=.3) +
  geom_text(data=dofs, aes(-0.07, y, label=name), hjust=1, size=2.7, fontface="italic", colour="grey15") +
  geom_text(data=tars, aes(TX + 0.05, y, label=tsym, fontface=ifelse(tsym == "YUC5", "bold.italic", "italic")),
            hjust=0, size=2.6, colour="grey15") +
  geom_segment(data=gb, aes(x=BARX, xend=BARX, y=lo, yend=hi), colour="grey65", linewidth=.35) +
  geom_text(data=gb, aes(BARX - 0.04, mid, label=lab), parse=TRUE, hjust=1, size=2.5, colour="grey15") +
  geom_segment(data=cb, aes(x=CBARX, xend=CBARX, y=lo, yend=hi, colour=I(col)), linewidth=.6) +
  geom_text(data=cb, aes(CBARX + 0.05, mid, label=class), hjust=0, size=2.5, colour="grey15") +
  annotate("text", x=-0.25, y=1.075, label="DOFs", size=2.9, fontface="bold", colour="grey10") +
  annotate("text", x=TX + 0.18, y=1.075, label="Predicted DOF targets", size=2.9, fontface="bold", colour="grey10") +
  scale_size_continuous(range=c(2.2, 4.2)) +
  scale_x_continuous(limits=c(-1.30, 1.85), expand=c(0, 0)) +
  scale_y_continuous(limits=c(-0.04, 1.11), expand=c(0, 0)) +
  coord_cartesian(clip="off") +
  theme_void(base_size=7) +
  # no size key, as in the PI's panel (it collides with the DOF column at this size); the caption states
  # that node size = number of displayed auxin genes reached (min/max in the stats CSV)
  theme(legend.position="none", plot.margin=margin(4, 4, 4, 4))

ggsave(file.path(OUT_DIR,"Fig6b_DOF_auxin_network.pdf"), p, width=125, height=78, units="mm", device=cairo_pdf)
ggsave(file.path(OUT_DIR,"Fig6b_DOF_auxin_network.png"), p, width=125, height=78, units="mm", dpi=400, bg="white")

## stats for captions/text, so no number is typed by hand
gc <- table(dofs$grp)
st <- data.frame(n_dof=nrow(dofs), n_targets=nrow(tars), n_edges_drawn=nrow(ks),
                 n_edges_all_targets=nrow(k), n_WT_only=gc[["WT only"]], n_fer4_only=gc[["fer-4 only"]],
                 n_both=gc[["Both"]], n_dof_to_YUC5=length(yuc5_src),
                 dof_without_YUC5=paste(sort(dofs$name[!dofs$yuc5]), collapse="; "),
                 dofs=paste(dofs$name, collapse="; "), targets=paste(tars$tsym, collapse="; "),
                 min_targets_reached=min(dofs$size), max_targets_reached=max(dofs$size))
write.csv(st, file.path(OUT_DIR,"Fig6b_DOF_auxin_network_stats.csv"), row.names=FALSE)
print(t(st))
