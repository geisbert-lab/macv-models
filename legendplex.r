#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(limma))
theme_set(ggpubr::theme_pubr() +
          theme(legend.background=element_rect(fill="transparent"),
                legend.position="right"))

# helper variables
cols.dates <- c("Baseline"="white",
                "4"="#ffffb2",
                "7"="#fd8d3c", 
                "10"="#bd0026",
                "≥14/Term"="black")
cols.species <- c(AGM="#cab2d6", Cyno="#6a3d9a")
cols.strain <- c("9530537"="#b2df8a", "Malale"="#33a02c")
shapes.species <- c(AGM=21, Cyno=22)
shapes.strain <- c("9530537"=21, "Malale"=22)
shapes.outcome <- c(Succumbed=21, Survived=22)
dpi.dates <- c("Baseline"=0, 
               "4"=4, 
               "7"=7, 
               "10"=10, 
               "≥14/Term"=14)

# for volcano plots
size.reg.time <- c(Up=1, Down=1)
cols.reg.time <- c(Up="#e41a1c", Down="#377eb8")
size.reg.strain <- c("9530537"=1, "Malale"=1)

# figures
fig6 <- list()
supp12 <- list()
supp13 <- list()
supp14 <- list()
supp15 <- list()
supp16 <- list()
supp17 <- list()

# helper functions
plot.pca.study <- function(study, meta.mat, logconc.mat, cols=cols.dates) {
  # subset metadata
  meta.mat <- filter(meta.mat, Study==study)
  
  # log-transform, remove batch (plate) effect, and run PCA
  pca <- logconc.mat[, meta.mat$ID] %>%
         t() %>% 
         prcomp()
  
  # extract PC %s
  pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
  pcs <- round(100*pcs)
  pcs <- paste0(names(pcs), " (", pcs, "%)")
  
  # plot it
  pca$x %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(meta.mat, by="ID") %>%
    ggplot(aes(PC1, PC2)) +
    geom_vline(xintercept=0, linetype=3, col="lightgrey") +
    geom_hline(yintercept=0, linetype=3, col="lightgrey") +
    geom_point(aes(fill=Daterange), size=2, pch=21) +
    scale_fill_manual(NULL, values=cols) +
    labs(x=pcs[1], y=pcs[2],
         title=study) +
    guides(fill=guide_legend(override.aes=list(pch=21, size=5)),
           shape=guide_legend(override.aes=list(fill="black")))
}
plot.pca.species <- function(species, meta.mat, conc.mat, 
                             cols=cols.dates, shapes=shapes.strain) {
  # subset metadata
  meta.mat <- filter(meta.mat, Species==species)
  
  # log-transform, remove batch (plate) effect, and run PCA
  pca <- model.matrix(~0 + Daterange, data=meta.mat)
  pca <- log2(conc.mat[, meta.mat$ID] + 0.001) %>%
         removeBatchEffect(batch=meta.mat$Batch, design=pca) %>%
         t() %>% 
         prcomp()
  
  # extract PC %s
  pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
  pcs <- round(100*pcs)
  pcs <- paste0(names(pcs), " (", pcs, "%)")
  
  # plot it
  pca$x %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(meta.mat, by="ID") %>%
    ggplot(aes(PC1, PC2)) +
    geom_vline(xintercept=0, linetype=3, col="lightgrey") +
    geom_hline(yintercept=0, linetype=3, col="lightgrey") +
    geom_point(aes(shape=Strain, fill=Daterange), size=2) +
    scale_shape_manual(NULL, values=shapes.strain) + 
    scale_fill_manual(NULL, values=cols) +
    labs(x=pcs[1], y=pcs[2],
         title=species) +
    guides(fill=guide_legend(override.aes=list(pch=21, size=5)),
           shape=guide_legend(override.aes=list(fill="black", size=5)))
}
get.results <- function(coef, diffexpr, psig=0.05, fsig=1) {
  # get results and loosely format
  resmat <- topTable(diffexpr, coef=coef, number=Inf, sort.by="logFC") %>%
    rownames_to_column("Protein") %>%
    mutate(coef=coef,
           psig=(adj.P.Val < psig),
           fsig=(abs(logFC) > fsig),
           Significant=(psig & fsig)) %>%
    rename(lfc=logFC, padj=adj.P.Val) %>%
    select(Protein, lfc, padj, Significant, coef)
  # make a "regulation" column
  resmat$Regulation <- "None"
  resmat$Regulation[resmat$Significant & resmat$lfc > 0] <- "Up"
  resmat$Regulation[resmat$Significant & resmat$lfc < 0] <- "Down"
  # return the formatted data.frame
  return(resmat)
}
plot.spaghetti <- function(lfc.mat) {
  lfc.mat %>%
    ggplot(aes(DPI, lfc)) +
    geom_hline(yintercept=0, linetype=3, col="grey60") +
    geom_line(aes(group=NHP, col=NHP), lineend="round") +
    scale_color_brewer(NULL, palette="Dark2") +
    scale_x_continuous("Days postinfection", breaks=c(0, 7, 14, 21)) +
    facet_wrap(~Protein, ncol=4, scales="free_y") +
    labs(y="Fold change (log2)",
         title=paste0(unique(lfc.mat$Species), ", ", unique(lfc.mat$Strain))) +
    guides(color=guide_legend(direction="horizontal", nrow=2,
                              override.aes=list(linewidth=2))) +
    theme(legend.position=c(0.6, 0.05))
}

## load data -------------------------------------------------------------------
# metadata
meta <- readxl::read_excel("samplesheet.xlsx", "animals")
meta <- readxl::read_excel("samplesheet.xlsx", "samples.legendplex") %>%
        left_join(meta, by="NHP") %>%
        filter((Outcome=="Succumbed" | Daterange=="Baseline"),
               DPI!=3) %>% # remove 3 DPI samples
        mutate(Daterange=factor(Daterange, levels=names(cols.dates)),
               Strain=factor(Strain, levels=names(shapes.strain)),
               Species=factor(Species, levels=names(shapes.species))) %>%
        as.data.frame()
rownames(meta) <- meta$ID

# predicted concentration 
cmat <- readxl::read_excel("samplesheet.xlsx", "data.legendplex") %>%
        # average technical replicates
        pivot_longer(cols=-ID, names_to="Analyte", values_to="Concentration") %>%
        group_by(Analyte, ID) %>%
        summarise(Concentration=mean(Concentration),
                  .groups="drop") %>%
        # re-format as a matrix
        pivot_wider(names_from=ID, values_from=Concentration) %>%
        column_to_rownames("Analyte") %>%
        as.matrix()

# align rows and columns
x <- intersect(colnames(cmat), rownames(meta))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

## clustering ------------------------------------------------------------------
# plot by species (both strains)
supp12$a <- plot.pca.species("Cyno", meta, cmat)
supp12$b <- plot.pca.species("AGM", meta, cmat)

# plot by study (single strain, single species)
# use pooled log2 concentration matrix with batch removed; since there are
# insufficient "n" otherwise
log.cmat <- model.matrix(~0 + Daterange, data=meta)
log.cmat <- log2(cmat + 0.001) %>% # add pseudocount so the minimum log.conc=0
            removeBatchEffect(batch=meta$Batch, design=log.cmat)
supp12$c <- plot.pca.study("Cyno, Malale", meta, log.cmat)
supp12$d <- plot.pca.study("AGM, Malale", meta, log.cmat)
supp12$e <- plot.pca.study("Cyno, 9530537", meta, log.cmat)
supp12$f <- plot.pca.study("AGM, 9530537", meta, log.cmat)

## differential expression over time, by study ---------------------------------
# using log2 concentration matrix with batch effect regressed out

# for each study, we're comparing time points postinfection to baseline, 
# controlling for within-NHP variation (paired analysis)
# This means our groupings are the interaction between species, strain, and DPI
complevels <- c("Baseline"="Baseline", "4"="D4", "7"="D7", 
                "10"="D10", "≥14/Term"="D14Term")
meta <- meta %>%
        mutate(Daterange=factor(Daterange, 
                                levels=names(complevels), 
                                labels=complevels),
               Group=interaction(Species, Strain, Daterange))
modmat <- model.matrix(~0 + Group, data=meta)
colnames(modmat) <- str_remove(colnames(modmat), "^Group")

# define contrasts
ctrs <- meta %>%
        filter(DPI > 0) %>%
        mutate(Baseline=paste0(Species, ".", Strain, ".Baseline"),
               Contrast=paste(Group, "-", Baseline))
ctrs <- makeContrasts(contrasts=unique(ctrs$Contrast), 
                      levels=levels(meta$Group))

# get duplicate correlation since we're running paired analyses
dc <- duplicateCorrelation(log.cmat, modmat, block=meta$NHP)

# run limma
rmat <- log.cmat %>%
        lmFit(modmat, block=meta$NHP, correlation=dc$consensus.correlation) %>%
        contrasts.fit(ctrs) %>%
        eBayes()

# extract resuls, return, and annotate
rmat <- colnames(rmat) %>%
        lapply(get.results, diffexpr=rmat) %>%
        do.call(rbind, .) %>%
        # extract species, strain, and daterange
        mutate(Species=str_extract(coef, "^Cyno|^AGM"),
               Strain=str_extract(coef, "(?<=\\.)9530537|Malale"),
               Daterange=str_extract(coef, "(?<=\\.)[A-z0-9]+(?= )"),
               Daterange=factor(Daterange, 
                                levels=complevels, 
                                labels=names(complevels))) %>%
        select(-coef)
write.csv(rmat, "legendplex/results-study.csv", row.names=FALSE)
rm(dc, modmat, ctrs, complevels)

# plot AGM and cyno, by strain
proteins <- rmat %>%
            filter(Significant) %>%
            group_by(Species, Strain, Daterange, Regulation) %>%
            slice_min(order_by=padj, n=5, with_ties=FALSE) %>%
            ungroup()
supp13$a <- rmat %>%
            filter(Species=="AGM") %>%
            ggplot(aes(lfc, -log10(padj))) +
            geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
            ggrepel::geom_text_repel(data=filter(proteins, Species=="AGM"),
                                     aes(label=Protein),
                                     size=3, box.padding=0.4, force=2, 
                                     min.segment.length=0.2) +
            scale_size_manual(values=size.reg.time, na.value=0.25) +
            scale_color_manual(values=cols.reg.time, na.value="grey80") +
            scale_x_continuous("Fold change (log2)") +
            scale_y_continuous("FDR-adjusted p-value (-log10)") +
            facet_grid(Strain ~ Daterange) +
            labs(title = "AGM protein changes over time, by strain") +
            theme(legend.position="none")
supp13$b <- rmat %>%
            filter(Species=="Cyno") %>%
            ggplot(aes(lfc, -log10(padj))) +
            geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
            ggrepel::geom_text_repel(data=filter(proteins, Species=="Cyno"),
                                     aes(label=Protein),
                                     size=3, box.padding=0.4, force=2, 
                                     min.segment.length=0.2) +
            scale_size_manual(values=size.reg.time, na.value=0.25) +
            scale_color_manual(values=cols.reg.time, na.value="grey80") +
            scale_x_continuous("Fold change (log2)") +
            scale_y_continuous("FDR-adjusted p-value (-log10)") +
            facet_grid(Strain ~ Daterange) +
            labs(title = "Cyno protein changes over time, by strain") +
            theme(legend.position="none")
rm(proteins)

## fold change over time by NHP ------------------------------------------------
# using all upregulation proteins, regardless of species or strain
proteins <- rmat %>%
            filter(Regulation=="Up") %>%
            select(Protein) %>%
            unlist() %>%
            unique()

# first, subset and melt the data
lfcmat <- log.cmat %>% 
          as.data.frame() %>%
          rownames_to_column("Protein") %>%
          filter(Protein %in% proteins) %>%
          pivot_longer(-Protein, names_to="ID", values_to="log.conc") %>%
          left_join(meta, by="ID")
# extract baseline
bl <- lfcmat %>%
      filter(Daterange=="Baseline") %>%
      rename(Baseline=log.conc) %>%
      select(NHP, Protein, Baseline)
# add baseline and calculate lfc (take difference since we're operating in log)
lfcmat <- lfcmat %>%
          left_join(bl, by=c("NHP", "Protein")) %>%
          mutate(lfc=log.conc-Baseline)

# plot each study separately 
supp14 <- plot.spaghetti(filter(lfcmat, Species=="AGM", Strain=="9530537"))
supp15 <- plot.spaghetti(filter(lfcmat, Species=="AGM", Strain=="Malale"))
supp16 <- plot.spaghetti(filter(lfcmat, Species=="Cyno", Strain=="9530537"))
supp17 <- plot.spaghetti(filter(lfcmat, Species=="Cyno", Strain=="Malale"))

# clean up
rm(proteins, bl)

## plot average FC over time for selected analytes -----------------------------
selected <- c("IP-10", "MCP-1", "MRP8/14", "tPA", "Eotaxin", "sST2", "IL-6")

# build baseline zeros
bl <- expand.grid(Species=names(cols.species),
                  Strain=names(cols.strain),
                  Protein=selected) %>%
      mutate(Daterange="Baseline",
             lfc=0)

# add baseline zeros to results matrix and plot
fig6$d <- rmat %>%
          filter(Protein %in% selected) %>%
          select(colnames(bl)) %>%
          rbind(bl) %>%
          mutate(Study=paste0(Species, ", ", Strain)) %>%
          ggplot(aes(Daterange, lfc)) +
          geom_hline(yintercept=0, linetype=3, color="lightgrey") + 
          geom_line(aes(group=Protein, col=Protein)) +
          scale_color_brewer(NULL, palette="Dark2") +
          scale_y_continuous("Fold change (log2)", limits=c(-2, 10)) +
          facet_wrap(~Study, nrow=1) +
          labs(x=NULL, y="Fold change (log2)",
               title="Circulating inflammatory mediators") +
          theme(axis.text.x=element_text(angle=45, hjust=1))

# clean up
rm(selected, bl)

# DE proteins between species --------------------------------------------------
# since this is a species-specific question, use peak lfc per NHP and pool
# strains
lfcmat <- lfcmat %>%
          group_by(Protein, NHP, Species, Strain) %>%
          slice_max(n=1, order_by=lfc, with_ties=FALSE) %>%
          ungroup()

# compare peak LFC via mann-whitney U test
pval <- lfcmat %>%
        group_by(Protein) %>%
        rstatix::wilcox_test(lfc ~ Species) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(scales="free_y") %>%
        filter(p < 0.05)

# plot it
fig6$e <- lfcmat %>%
          filter(Protein %in% pval$Protein) %>%
          ggplot(aes(Species, lfc)) +
          geom_hline(yintercept=0, linetype=3, col="lightgrey") +
          geom_boxplot(aes(fill=Species), col="black", alpha=0.5, 
                       outliers=FALSE) +
          geom_jitter(aes(shape=Strain, fill=Species), 
                      width=0.2, height=0, size=0.75) +
          ggpubr::stat_pvalue_manual(data=pval) +
          scale_shape_manual(NULL, values=shapes.strain) +
          scale_fill_manual(NULL, values=cols.species) +
          scale_y_continuous("Fold change (log2)", limits=c(0, NA),
                             expand=expansion(c(0.05, 0.15))) +
          facet_wrap(~Protein, nrow=1, scales="free_y") +
          labs(x=NULL, title="Proteins with significant differences in peak levels between species") +
          guides(shape=guide_legend(override.aes=list(fill="black", size=3)),
                 fill=guide_legend(override.aes=list(shape=21))) +
          theme(axis.text.x=element_text(angle=30, hjust=1))

# clean up
rm(pval, lfcmat)

## assemble figures ------------------------------------------------------------
# main figure d-e
cowplot::plot_grid(plotlist=fig6, labels=c("d", "e"), 
                   ncol=1, rel_heights=c(5, 4))
ggsave("figure-main06-de.png", units="in", width=7.5, height=4.5)

# supplemental figure 12
cowplot::plot_grid(plotlist=supp12, labels="auto", ncol=2)
ggsave("figure-supplementary12.png", units="in", width=7.5, height=8)

# supplemental figure 13
cowplot::plot_grid(plotlist=supp13, labels="auto", ncol=1)
ggsave("figure-supplementary13.png", units="in", width=7.5, height=8)

# supplemental figure 14-17
ggsave("figure-supplementary14.png", plot=supp14, 
       units="in", width=7.5, height=10)
ggsave("figure-supplementary15.png", plot=supp15, 
       units="in", width=7.5, height=10)
ggsave("figure-supplementary16.png", plot=supp16, 
       units="in", width=7.5, height=10)
ggsave("figure-supplementary17.png", plot=supp17, 
       units="in", width=7.5, height=10)

## fin! ------------------------------------------------------------------------
sessionInfo()
