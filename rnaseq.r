#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
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
cols.reg <- c(Up="#e41a1c", Down="#377eb8")
shapes.species <- c(AGM=21, Cyno=22)
shapes.strain <- c("9530537"=21, "Malale"=22)
size.reg <- c(Up=0.5, Down=0.5)
size.strain <- c("9530537"=0.5, Malale=0.5) 
size.species <- c(AGM=0.5, Cyno=0.5)
dpi.dates <- c("Baseline"=0, 
               "4"=4, 
               "7"=7, 
               "10"=10, 
               "≥14/Term"=14)

# helper functions
get.results <- function(group, de.obj, th.lfc=2, th.padj=0.05) {
  # get basic results
  r <- results(de.obj, name=group, 
               cooksCutoff=FALSE, alpha=0.05, independentFiltering=FALSE)
  # shrink results
  r <- lfcShrink(de.obj, coef=group, res=r, type="apeglm", quiet=TRUE) %>%
       # apply some formatting
       as.data.frame() %>%
       rownames_to_column("Gene") %>%
       dplyr::rename(lfc=log2FoldChange) %>%
       # remove genes with NA lfc
       filter(!is.na(lfc)) %>%
       # add group and a significance column
       mutate(Group=group, 
              psig=(!is.na(padj) & (padj < th.padj)),
              fsig=(abs(lfc) > th.lfc),
              Significant=(psig & fsig)) %>%
       select(Gene, lfc, padj, Group, Significant)
  # add regulation column
  r$Regulation <- "None"
  r$Regulation[r$Significant & (r$lfc < 0)] <- "Down"
  r$Regulation[r$Significant & (r$lfc > 0)] <- "Up"
  # return matrix
  return(r)
}

## load metadata ---------------------------------------------------------------
# sample info: succumbed only and must pass QC
meta <- readxl::read_excel("samplesheet.xlsx", "animals")
meta <- readxl::read_excel("samplesheet.xlsx", "samples") %>%
        left_join(meta, by="NHP") %>%
        filter(QC, 
               (Outcome=="Succumbed" | Daterange=="Baseline"),
               DPI!=3) %>%
        mutate(Daterange=factor(Daterange, levels=names(cols.dates)),
               Strain=factor(Strain, levels=names(shapes.strain)),
               Species=factor(Species, levels=names(shapes.species))) %>%
        as.data.frame()
rownames(meta) <- meta$ID

## determine best quantification -----------------------------------------------
# cyno
qunt <- read.delim("rnaseq/mfascicularis.tsv.summary") %>%
        filter(Status=="Assigned") %>%
        reshape2::melt(id.vars="Status",
                       variable.name="ID",
                       value.name="Reads") %>%
        right_join(meta, by="ID") %>%
        mutate(Genome="Cyno")

# AGM
qunt <- read.delim("rnaseq/csabaeus.tsv.summary") %>%
        filter(Status=="Assigned") %>%
        reshape2::melt(id.vars="Status",
                       variable.name="ID",
                       value.name="Reads") %>%
        right_join(meta, by="ID") %>%
        mutate(Genome="AGM") %>%
        rbind(qunt)

# for AGMs, which genome is better? -- cyno
agms <- qunt %>%
        filter(Species=="AGM") %>%
        rstatix::wilcox_test(Reads ~ Genome, paired=TRUE, id="ID") %>%
        rstatix::add_significance() %>%
        rstatix::add_x_position()
agms <- qunt %>%
        filter(Species=="AGM") %>%
        ggplot(aes(Genome, Reads)) +
        geom_violin(col="black", fill="lightgrey") +
        geom_jitter(col="black", height=0, width=0.2, size=0.5) +
        ggpubr::stat_pvalue_manual(data=agms, y.position=1e8) +
        scale_y_continuous("Quantified reads", limits=c(1e6, 2e8), 
                           transform="log10") +
        labs(x="Alignment genome", title="AGM quantification")
# for cynos, which genome is better? -- cyno
cyno <- qunt %>%
        filter(Species=="Cyno") %>%
        rstatix::wilcox_test(Reads ~ Genome, paired=TRUE, id="ID") %>%
        rstatix::add_significance() %>%
        rstatix::add_x_position()
cyno <- qunt %>%
        filter(Species=="Cyno") %>%
        ggplot(aes(Genome, Reads)) +
        geom_violin(col="black", fill="lightgrey") +
        geom_jitter(col="black", height=0, width=0.2, size=0.5) +
        ggpubr::stat_pvalue_manual(data=cyno, y.position=1e8) +
        scale_y_continuous("Quantified reads", limits=c(1e6, 2e8), 
                           transform="log10") +
        labs(x="Alignment genome", title="Cyno quantification")
cowplot::plot_grid(cyno, agms, nrow=1)
ggsave("rnaseq/genomes.png", units="in", width=7.5, height=2.5)

# clean up
rm(qunt, agms, cyno)

## load counts -----------------------------------------------------------------
cmat <- read.delim("rnaseq/mfascicularis.tsv") %>%
        reshape2::melt(id.vars="Gene",
                       variable.name="ID", 
                       value.name="Counts") %>%
        # remove rRNA, y-RNA, and non-named genes (e.g., ENSMFAG*)
        filter(!str_detect(Gene, "rRNA"),
               Gene != "Y_RNA",
               !str_detect(Gene, "^ENSMFAG"),
               !str_detect(Gene, "^SNO|^sno")) %>%
        # collapse viral genes to "MACV"
        mutate(Gene=if_else(Gene %in% c("Z", "L", "GPC", "N"), 
                            "MACV", Gene)) %>%
        # collapse duplicate gene names
        group_by(ID, Gene) %>%
        summarise(Counts=sum(Counts),
                  .groups="drop") %>%
        reshape2::dcast(Gene ~ ID, value.var="Counts") %>%
        # format as a matrix
        column_to_rownames("Gene") %>%
        as.matrix()

# align rows and columns
x <- intersect(rownames(meta), colnames(cmat))
meta <- meta[x, ]
cmat <- cmat[, x]
mode(cmat) <- "integer"
rm(x)

# remove genes with sparse counts: at least 5 samples with count > 10
keepgenes <- rowSums(cmat > 10)
keepgenes <- (keepgenes >= 5)
cmat <- cmat[keepgenes, ]
rm(keepgenes)

# save CPM (no log) for CIBERSORTx
cmat %>%
  edgeR::DGEList() %>%
  edgeR::normLibSizes() %>%
  edgeR::cpm() %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  write.table("rnaseq/cibersort-input.tsv", sep="\t")
# once run, load CIBERSORT output
cibx <- read.csv("rnaseq/cibersort-output.csv", check.names=FALSE) %>%
        filter(`P-value` < 0.05) %>%
        select(-`P-value`, -Correlation, -RMSE, 
               -`Absolute score (sig.score)`,
               # remove eosinophils, activated NK
               -Eosinophils, -`NK cells activated`,
               # remove macrophages (sould be tissue-resident),
               -starts_with("Macrophages")) %>%
        dplyr::rename(ID="Mixture") %>%
        reshape2::melt(id.vars="ID",
                       variable.name="Celltype",
                       value.name="Score") %>%
        right_join(meta, by="ID")

# samples per bin -- will need to remove cyno '95 D10 from comparisons, n<3
meta %>%
  filter(Outcome=="Succumbed") %>%
  mutate(Group=paste0(Species, ", MACV ", Strain)) %>%
  group_by(Group, Daterange) %>%
  summarise(Samples=n(),
            .groups="drop") %>%
  arrange(Samples)

## clustering and QC -----------------------------------------------------------
# PCA by study (single strain, single species) -- looks good
unique(meta$Study) %>%
  lapply(function(i) {
    m <- filter(meta, Study==i)
    pca <- DESeqDataSetFromMatrix(cmat[, rownames(m)], m, ~1) %>%
           varianceStabilizingTransformation() %>%
           plotPCA(intgroup="Daterange", returnData=TRUE)
    pcs <- round(100*attr(pca, "percentVar"))
    pcs <- paste0(c("PC1 (", "PC2 ("), pcs, "%)")
    pca %>%
      ggplot(aes(PC1, PC2)) +
      geom_vline(xintercept=0, linetype=3, col="lightgrey") +
      geom_hline(yintercept=0, linetype=3, col="lightgrey") +
      geom_point(aes(fill=Daterange), size=2, pch=21) +
      scale_fill_manual(NULL, values=cols.dates) +
      labs(x=pcs[1], y=pcs[2],
           title=i) +
      guides(fill=guide_legend(override.aes=list(pch=21, size=5)),
             shape=guide_legend(override.aes=list(fill="black")))
  }) %>%
  cowplot::plot_grid(plotlist=.)
ggsave("rnaseq/pca-split-study.png", units="in", width=7.5, height=4)

# what is the strain-specific effect for each species, if any? (combine strains)
unique(meta$Species) %>%
  lapply(function(i) {
    m <- filter(meta, Species==i)
    pca <- DESeqDataSetFromMatrix(cmat[, rownames(m)], m, ~1) %>%
           varianceStabilizingTransformation() %>%
           plotPCA(intgroup=c("Strain", "Daterange"), returnData=TRUE)
    pcs <- round(100*attr(pca, "percentVar"))
    pcs <- paste0(c("PC1 (", "PC2 ("), pcs, "%)")
    pca %>%
      ggplot(aes(PC1, PC2)) +
      geom_vline(xintercept=0, linetype=3, col="lightgrey") +
      geom_hline(yintercept=0, linetype=3, col="lightgrey") +
      geom_point(aes(shape=Strain, fill=Daterange), size=2) +
      scale_shape_manual(NULL, values=shapes.strain) +
      scale_fill_manual(NULL, values=cols.dates) +
      labs(x=pcs[1], y=pcs[2],
           title=i) +
      guides(fill=guide_legend(override.aes=list(pch=21, size=5)),
             shape=guide_legend(override.aes=list(fill="black", size=5)))
  }) %>%
  cowplot::plot_grid(plotlist=.)
ggsave("rnaseq/pca-split-species.png", units="in", width=7.5, height=2.5)

## DE #1: models ---------------------------------------------------------------
# compare each time point to baseline, subsetting by species and strain
# use all baseline samples for a given NHP species for higher N

# AGM, Malale 
de <- meta %>%
      filter(Species=="AGM", 
             Strain=="Malale" | Daterange=="Baseline")
de <- DESeqDataSetFromMatrix(cmat[, rownames(de)], de, ~Daterange) %>%
      DESeq(quiet=TRUE)
rmat <- resultsNames(de)[-1] %>%
        lapply(get.results, de.obj=de) %>%
        do.call(rbind, .) %>%
        mutate(Species="AGM",
               Strain="Malale") 

# AGM, 9530537
de <- meta %>%
      filter(Species=="AGM", 
             Strain=="9530537" | Daterange=="Baseline")
de <- DESeqDataSetFromMatrix(cmat[, rownames(de)], de, ~Daterange) %>%
      DESeq(quiet=TRUE)
rmat <- resultsNames(de)[-1] %>%
        lapply(get.results, de.obj=de) %>%
        do.call(rbind, .) %>%
        mutate(Species="AGM",
               Strain="9530537") %>%
        rbind(rmat)

# cyno, Malale
de <- meta %>%
      filter(Species=="Cyno", 
             Strain=="Malale" | Daterange=="Baseline")
de <- DESeqDataSetFromMatrix(cmat[, rownames(de)], de, ~Daterange) %>%
      DESeq(quiet=TRUE)
rmat <- resultsNames(de)[-1] %>%
        lapply(get.results, de.obj=de) %>%
        do.call(rbind, .) %>%
        mutate(Species="Cyno",
               Strain="Malale") %>%
        rbind(rmat)

# cyno, 9530537 -- remove D10 since n=2
de <- meta %>%
      filter(Species=="Cyno", 
             (Strain=="9530537" | Daterange=="Baseline"),
             Daterange!="10")
de <- DESeqDataSetFromMatrix(cmat[, rownames(de)], de, ~Daterange) %>%
      DESeq(quiet=TRUE)
rmat <- resultsNames(de)[-1] %>%
        lapply(get.results, de.obj=de) %>%
        do.call(rbind, .) %>%
        mutate(Species="Cyno",
               Strain="9530537") %>%
        rbind(rmat)

# extract daterange from group column
rmat <- rmat %>%
        mutate(Daterange=factor(Group,
                                levels=make.names(paste0("Daterange_", 
                                                         names(cols.dates), 
                                                         "_vs_Baseline")),
                                labels=names(cols.dates))) %>%
        select(-Group)

# save results and clean up
rm(de)
write.csv(rmat, "rnaseq/results-study.csv", row.names=FALSE)

# save for IPA
rmat %>%
  select(Gene, Daterange, Species, Strain, lfc, padj) %>%
  reshape2::melt(id.vars=c("Gene", "Daterange", "Species", "Strain")) %>%
  mutate(Group=interaction(Species, Strain, Daterange, variable)) %>%
  reshape2::dcast(Gene ~ Group, value.var="value") %>%
  write.csv("rnaseq/ipa-input-study.csv", row.names=FALSE)

# ALT OPTION: load pre-saved results matrix
#rmat <- read.csv("rnaseq/results-study.csv") %>%
#        mutate(Daterange=factor(Daterange, levels=names(cols.dates)))

# total DE genes per group
x <- expand.grid(Species=names(shapes.species),
            Strain=names(shapes.strain)) %>%
     mutate(Daterange="Baseline",
            Genes=0)
rmat %>%
  filter(Significant) %>%
  group_by(Species, Strain, Daterange) %>%
  summarize(Genes=n(),
            .groups="drop") %>%
  # add zeros
  rbind(x) %>%
  # set daterange to numeric to get an accurate curve
  mutate(DPI=factor(Daterange, levels=names(dpi.dates), labels=dpi.dates),
         DPI=as.integer(as.character(DPI))) %>%
  # plot it
  ggplot(aes(DPI, Genes)) +
  geom_line(aes(group=interaction(Species, Strain))) +
  geom_point(aes(shape=Strain, fill=Species), size=3) +
  scale_shape_manual(values=shapes.strain) +
  scale_fill_manual(values=cols.species) +
  scale_x_continuous("Days postinfection", breaks=dpi.dates) +
  labs(y="Total DE genes",
       title="mRNA response magnitude") +
  guides(fill=guide_legend(override.aes=list(pch=21, size=5)),
         shape=guide_legend(override.aes=list(fill="black", size=5)))
ggsave("rnaseq/degenes-study.png", units="in", width=3.75, height=2.5)
rm(x)

# capped volcano plots
rmat$padj[rmat$padj < 1e-50] <- 1e-50
genelist <- rmat %>%
            filter(Significant) %>%
            group_by(Strain, Species, Daterange, Regulation) %>%
            top_n(n=5, wt=abs(lfc)) %>%
            ungroup()
# cyno
rmat %>%
  filter(Species=="Cyno") %>% 
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation), alpha=0.5) +
  ggrepel::geom_text_repel(data=filter(genelist, Species=="Cyno"), 
                           aes(label=Gene), size=3) +
  scale_size_manual(values=size.reg, na.value=0.25) +
  scale_color_manual(values=cols.reg, na.value="lightgrey") +
  scale_x_continuous("Fold change (log2)", breaks=c(-10, -5, 0, 5, 10)) +
  scale_y_continuous("FDR-adjusted p-value (-log10)",
                     breaks=c(0, 10, 20, 30, 40, 50)) +
  facet_grid(Strain ~ Daterange) +
  labs(x="Fold change (log2)",
       title="Cyno mRNA changes over time, by strain") +
  guides(size="none", col="none")
ggsave("rnaseq/volcano-study-cyno.png", units="in", width=7.5, height=4)
# AGM
rmat %>%
  filter(Species=="AGM") %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation), alpha=0.5) +
  ggrepel::geom_text_repel(data=filter(genelist, Species=="AGM"), 
                           aes(label=Gene), size=3) +
  scale_size_manual(values=size.reg, na.value=0.25) +
  scale_color_manual(values=cols.reg, na.value="lightgrey") +
  scale_x_continuous("Fold change (log2)", breaks=c(-10, -5, 0, 5, 10)) +
  scale_y_continuous("FDR-adjusted p-value (-log10)",
                     breaks=c(0, 10, 20, 30, 40, 50)) +
  facet_grid(Strain ~ Daterange) +
  labs(x="Fold change (log2)",
       title="AGM mRNA changes over time, by strain") +
  guides(size="none", col="none")
ggsave("rnaseq/volcano-study-agm.png", units="in", width=7.5, height=4)
rm(genelist)

# cell types - significant changes in populations over time
subtypes <- cibx %>%
            group_by(Celltype, Species, Strain) %>%
            # KW test first
            rstatix::kruskal_test(Score ~ Daterange) %>% 
            filter(p < 0.05) %>%
            select(Celltype, Strain, Species) %>%
            # add back data and run dunn post hoc. only care about baseline difs
            left_join(cibx, by=c("Celltype", "Strain", "Species")) %>%
            group_by(Celltype, Species, Strain) %>%
            rstatix::dunn_test(Score ~ Daterange, ref.group="Baseline") %>%
            filter(p.adj < 0.05) %>%
            rstatix::add_xy_position(scales="free_y")
# agm, malale
subtypes %>%
  filter(Species=="AGM", Strain=="Malale") %>%
  select(Celltype, Species, Strain) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(col="black", fill="lightgrey") +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, Species=="AGM", 
                                         Strain=="Malale"), 
                             color="#e41a1c") +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="DPI", y="Celltype score", title="AGM, Malale")
# agm, 9530537
subtypes %>%
  filter(Species=="AGM", Strain=="9530537") %>%
  select(Celltype, Species, Strain) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(col="black", fill="lightgrey") +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, Species=="AGM", 
                                         Strain=="9530537"), 
                             color="#e41a1c") +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="DPI", y="Celltype score", title="AGM, 9530537")
# cyno, malale
subtypes %>%
  filter(Species=="Cyno", Strain=="Malale") %>%
  select(Celltype, Species, Strain) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(col="black", fill="lightgrey") +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, Species=="Cyno", 
                                         Strain=="Malale"), 
                             color="#e41a1c") +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="DPI", y="Celltype score", title="Cyno, Malale")
# cyno, 9530537
subtypes %>%
  filter(Species=="Cyno", Strain=="9530537") %>%
  select(Celltype, Species, Strain) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(col="black", fill="lightgrey") +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, Species=="Cyno", 
                                         Strain=="9530537"), 
                             color="#e41a1c") +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="DPI", y="Celltype score", title="Cyno, 9530537")
rm(subtypes)

# selected marker genes
genelist <- c("ISG15", "CXCL11", "CCL8", "ICAM4", "OLFM4")
# add baseline for plotting
d0 <- expand.grid(Gene=genelist,
            Species=c("AGM", "Cyno"),
            Strain=c("9530537", "Malale"),
            DPI=0,
            lfc=0)
rmat %>%
  # set daterange to numeric to get an accurate curve
  mutate(DPI=factor(Daterange, levels=names(dpi.dates), labels=dpi.dates),
         DPI=as.integer(as.character(DPI))) %>%
  select(Gene, DPI, lfc, Species, Strain) %>%
  rbind(d0) %>%
  filter(Gene %in% genelist) %>% 
  mutate(Study=paste0(Species, ", ", Strain)) %>%
  ggplot(aes(DPI, lfc)) +
  geom_hline(yintercept=0, linetype=3, col="lightgrey") +
  geom_line(aes(col=Gene, group=Gene), lineend="round") +
  scale_color_brewer(NULL, palette="Dark2") +
  scale_x_continuous(NULL, breaks=dpi.dates) +
  facet_wrap(~Study, nrow=1) +
  labs(y="Fold change (log2)", title="Conserved mRNAs") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("rnaseq/lfc-study.png", units="in", width=7.5, height=2.5)
rm(genelist, d0)

# clean up
rm(rmat)

## DE #2: differences between strains ------------------------------------------
# compare strains at each time-point after infection, controlling for 
# differences from baseline; analyze species separately
# using example from https://support.bioconductor.org/p/65676/

# AGM
de <- filter(meta, Species=="AGM")
de <- DESeqDataSetFromMatrix(cmat[, rownames(de)], de, 
                             ~Strain+Daterange+Strain:Daterange) %>%
      DESeq(test="LRT", reduced=~Strain+Daterange)
rmat <- resultsNames(de)[7:10] %>%
        lapply(get.results, de.obj=de) %>%
        do.call(rbind, .) %>%
        mutate(Species="AGM")

# cyno -- remove D10 since n=2
de <- filter(meta, Species=="Cyno", Daterange!="10")
de <- DESeqDataSetFromMatrix(cmat[, rownames(de)], de, 
                             ~Strain+Daterange+Strain:Daterange) %>%
      DESeq(test="LRT", reduced=~Strain+Daterange)
rmat <- resultsNames(de)[6:8] %>%
        lapply(get.results, de.obj=de) %>%
        do.call(rbind, .) %>%
        mutate(Species="Cyno") %>%
        # add in AGM data
        rbind(rmat)

# extract daterange and re-format regulation
rmat <- rmat %>%
        mutate(Daterange=factor(Group,
                                levels=make.names(paste0("StrainMalale.Daterange", 
                                                         names(cols.dates))), 
                                labels=names(cols.dates)),
               Regulation=factor(Regulation, 
                                 levels=c("Up", "None", "Down"),
                                 labels=c("Malale", "None", "9530537"))) %>%
        select(-Group)
# save results and clean up
rm(de)
write.csv(rmat, "rnaseq/results-strain.csv", row.names=FALSE)

# save for IPA
rmat %>%
  select(Gene, Daterange, Species, lfc, padj) %>%
  reshape2::melt(id.vars=c("Gene", "Daterange", "Species")) %>%
  mutate(Group=interaction(Species, Daterange, variable)) %>%
  reshape2::dcast(Gene ~ Group, value.var="value") %>%
  write.csv("rnaseq/ipa-input-strain.csv", row.names=FALSE)

# ALT OPTION: load pre-saved results matrix
#rmat <- read.csv("rnaseq/results-strain.csv") %>%
#        mutate(Daterange=factor(Daterange, levels=names(cols.dates)))

# strain-based total differences
x <- expand.grid(Species=names(shapes.species),
                 Regulation=names(shapes.strain),
                 Daterange=names(cols.dates)[-1]) %>%
     filter(!(Species=="Cyno" & Daterange=="10"))
rmat %>%
  filter(Significant) %>%
  group_by(Species, Regulation, Daterange) %>%
  summarize(Genes=n(),
            .groups="drop") %>%
  # add zeros
  right_join(x, by=c("Species", "Regulation", "Daterange")) %>%
  replace_na(list(Genes=0)) %>%
  # plot it
  ggplot(aes(Daterange, Genes)) +
  geom_col(aes(fill=Regulation), col="black", position="dodge") +
  geom_text(aes(label=Genes, y=Genes+8, 
                group=interaction(Daterange, Regulation)),
            position=position_dodge(width=0.8)) +
  scale_fill_manual(values=cols.strain) +
  scale_y_continuous("Significantly DE genes", limits=c(NA, 70)) +
  facet_wrap(~Species, ncol=1) +
  labs(x="Days postinfection",
       fill="Higher expr.",
       title="Strain differences, by species") +
  theme(legend.position=c(0.2, 0.82))
ggsave("rnaseq/degenes-strain.png", units="in", width=3.75, height=3.5)
rm(x)

# volcano plots
genelist <- rmat %>%
            filter(Significant) %>%
            group_by(Species, Daterange, Regulation) %>%
            top_n(n=10, wt=abs(lfc)) %>%
            ungroup()
rmat %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation), alpha=0.5) +
  ggrepel::geom_text_repel(data=genelist, aes(label=Gene), size=3) +
  scale_size_manual(values=size.strain, na.value=0.25) +
  scale_color_manual("Higher expr.", values=cols.strain, na.value="lightgrey") +
  facet_grid(Species ~ Daterange) +
  labs(x="Fold change (log2)", y="FDR-adjusted p-value (-log10)",
       title="Strain differences, by species") +
  guides(size="none", col=guide_legend(override.aes=list(size=5))) +
  theme(legend.position=c(0.65, 0.3))
ggsave("rnaseq/volcano-strain.png", units="in", width=7.5, height=4)
# plot just D14/term for main figure
rmat %>%
  filter(Daterange=="≥14/Term") %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation), alpha=0.5) +
  ggrepel::geom_text_repel(data=filter(genelist, Daterange=="≥14/Term"), 
                           aes(label=Gene), size=3) +
  scale_size_manual(values=size.strain, na.value=0.25) +
  scale_color_manual("Higher expr.", values=cols.strain, na.value="lightgrey") +
  facet_wrap(~Species) +
  labs(x="Fold change (log2)", y="FDR-adjusted p-value (-log10)",
       title="Strain differences, ≥14/terminal") +
  guides(size="none", col=guide_legend(override.aes=list(size=5))) +
  theme(legend.position=c(0.5, 0.8))
ggsave("rnaseq/volcano-strain-14term.png", units="in", width=3.75, height=3)
rm(genelist)

# cell types - significant changes in populations over time
subtypes <- cibx %>%
            filter(Daterange != "Baseline") %>%
            group_by(Celltype, Species, Daterange) %>%
            rstatix::wilcox_test(Score ~ Strain) %>%
            filter(p < 0.05) %>%
            rstatix::add_significance() %>%
            rstatix::add_xy_position(x="Daterange", 
                                     group="Strain", 
                                     scales="free_y")
# agm
subtypes %>%
  filter(Species=="AGM") %>%
  select(Celltype, Species) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(aes(group=interaction(Daterange, Strain), fill=Strain), 
               col="black", outlier.shape=21) +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, Species=="AGM")) +
  scale_fill_manual(values=cols.strain) +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="Days postinfection", y="Celltype score", 
       title="MACV strain differences (AGM)")
# cyno
subtypes %>%
  filter(Species=="Cyno") %>%
  select(Celltype, Species) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(aes(group=interaction(Daterange, Strain), fill=Strain), 
               col="black", outlier.shape=21) +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, Species=="Cyno")) +
  scale_fill_manual(values=cols.strain) +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="Days postinfection", y="Celltype score", 
       title="MACV strain differences (CM)")
# plot AGM plasma cell
subtypes %>%
  filter(Species=="AGM", Celltype=="Plasma cells") %>%
  select(Celltype, Species) %>%
  distinct() %>%
  left_join(cibx) %>%
  ggplot(aes(Daterange, Score)) +
  geom_boxplot(aes(group=interaction(Daterange, Strain), fill=Strain), 
               col="black", outlier.shape=21) +
  ggpubr::stat_pvalue_manual(data=filter(subtypes, 
                                         Species=="AGM", 
                                         Celltype=="Plasma cells")) +
  scale_fill_manual(values=cols.strain) +
  facet_wrap(~Celltype, scales="free_y") +
  labs(x="Days postinfection", y="Celltype score", 
       title="Strain celltype differences (AGM)") +
  theme(legend.position=c(0.2, 0.8))
ggsave("rnaseq/celltypes-strain-agm-plasma.png", 
       units="in", width=3.75, height=3.5)

# pathways
pathways <- readxl::read_excel("rnaseq/ipa-output.xlsx", "strain") %>%
            reshape2::melt(id.vars=c("Type", "Pathway"),
                           value.name="zscore") %>%
            na.omit() %>%
            # extract groups
            mutate(Species=str_extract(variable, "^[A-z]+"),
                   Daterange=str_extract(variable, "(?<=\\.)[^\\.]+$"),
                   Daterange=factor(Daterange, levels=names(cols.dates))) %>%
            select(-variable, -Type)


# clean up
rm(rmat, subtypes, pathways)

## DE #3: differences between species ------------------------------------------
# we are specifically interested in how species respond to disease differently
# we are NOT interested in baseline differences
# species "batch" could not be fully regressed. examining lfc correlations

# load study results matrix again
rmat <- read.csv("rnaseq/results-study.csv") %>%
        mutate(Daterange=factor(Daterange, levels=names(cols.dates)))

# correlate between species
# first, get the list of genes that are DE in one or both species
genelist <- rmat %>%
            filter(Significant, !(Strain=="9530537" & Daterange==10)) %>%
            select(Strain, Daterange, Gene) %>%
            distinct()
# second, identify whether these genes were DE in cyno, AGM, or both
genelist <- rmat %>%
            right_join(genelist) %>%
            reshape2::dcast(Gene + Strain + Daterange ~ Species, 
                            value.var="Significant") %>% 
            mutate(Group=factor(Cyno+2*AGM, levels=1:3, 
                                labels=c("Cyno", "AGM", "Both")),
                   Group=as.character(Group)) %>%
            select(Strain, Daterange, Gene, Group)
# third, find the genes with outlier differences (by strain)
# using rosner's test since we have very high n (~13.7k genes)
difs <- rmat %>%
        right_join(genelist) %>%
        reshape2::dcast(Gene + Strain + Daterange + Group ~ Species, 
                        value.var="lfc") %>%
        na.omit() %>%
        mutate(Difference=abs(Cyno-AGM))
# for 9530537
m95 <- difs[difs$Strain=="9530537", ]
x <- EnvStats::rosnerTest(m95$Difference, k=10)
x <- x$all.stats %>%
     filter(Outlier==TRUE) %>%
     select(Obs.Num) %>%
     unlist()
m95 <- m95[x, ]
# for malale
mal <- difs[difs$Strain=="Malale", ]
x <- EnvStats::rosnerTest(mal$Difference, k=10)
x <- x$all.stats %>%
     filter(Outlier==TRUE) %>%
     select(Obs.Num) %>%
     unlist()
mal <- mal[x, ]
outliers <- rbind(m95, mal)
rm(m95, mal, x)
# plot all for supplemental
rmat %>%
  right_join(genelist) %>%
  reshape2::dcast(Gene + Strain + Daterange + Group ~ Species, 
                  value.var="lfc") %>%
  na.omit() %>%
  ggplot(aes(Cyno, AGM)) +
  geom_abline(slope=1, intercept=0, col="darkgrey", linetype=3) +
  geom_point(aes(col=Group, size=Group)) +
  ggrepel::geom_text_repel(data=outliers, aes(label=Gene), size=3) +
  scale_size_manual(values=size.species, na.value=0.25) +
  scale_color_manual("Only DE in", values=cols.species, na.value="lightgrey") +
  ggpubr::stat_cor(aes(label=after_stat(rr.label)), label.y=-5, label.x=5) +
  facet_grid(Strain ~ Daterange) +
  guides(size="none", color=guide_legend(override.aes=list(size=5))) +
  labs(x="Cyno fold change (log2)", y="AGM fold change (log2)",
       title="mRNA correlation between species, by strain") +
  theme(legend.position=c(0.57, 0.85))
ggsave("rnaseq/lfc-species-corr.png",  units="in", width=7.5, height=4)
# plot just D14/terminal
rmat %>%
  right_join(genelist) %>%
  reshape2::dcast(Gene + Strain + Daterange + Group ~ Species, 
                  value.var="lfc") %>%
  na.omit() %>%
  filter(Daterange=="≥14/Term") %>%
  ggplot(aes(Cyno, AGM)) +
  geom_abline(slope=1, intercept=0, col="darkgrey", linetype=3) +
  geom_point(aes(col=Group, size=Group)) +
  ggrepel::geom_text_repel(data=outliers, aes(label=Gene), size=3) +
  scale_size_manual(values=size.species, na.value=0.25) +
  scale_color_manual("Only DE in", values=cols.species, na.value="lightgrey") +
  ggpubr::stat_cor(aes(label=after_stat(rr.label)), label.y=-5, label.x=5) +
  facet_wrap(~Strain) +
  guides(size="none", color=guide_legend(override.aes=list(size=5))) +
  labs(x="Cyno fold change (log2)", y="AGM fold change (log2)",
       title="Species correlation, ≥14/terminal") +
  theme(legend.position=c(0.87, 0.8))
ggsave("rnaseq/lfc-species-corr.png",  units="in", width=3.75, height=3)

# what % of DE genes per strain per timepoint are both/agm only/cyno only?
difs %>%
  group_by(Strain, Daterange, Group) %>%
  summarise(Genes=n(),
            .groups="drop") %>%
  group_by(Strain, Daterange) %>%
  mutate(Proportion=Genes/sum(Genes),
         Group=factor(Group, levels=c("AGM", "Cyno", "Both"))) %>%
  ungroup() %>%
  ggplot(aes(Daterange, Proportion)) +
  geom_hline(yintercept=1/3, linetype=3, col="lightgrey") +
  geom_line(aes(col=Group, group=interaction(Group, Strain)), linewidth=1) +
  geom_point(aes(col=Group), size=3) +
  scale_color_manual("DE in", values=c(cols.species, Both="darkgrey")) +
  scale_y_continuous("Proportion of DE genes", limits=c(0, 0.5)) +
  facet_wrap(~Strain) +
  labs(x="Days postinfection", 
       title="DE gene overlap between species, by strain")
ggsave("rnaseq/degenes-proportions.png", units="in", width=7.5, height=3)

## fin -------------------------------------------------------------------------
sessionInfo()
