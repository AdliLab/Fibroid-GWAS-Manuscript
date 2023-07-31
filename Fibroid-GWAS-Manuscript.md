Leiomyoma_RNAseq_GWAS_scRNAseq
================
Alexander J. Duval
2023-07-17

## Load environment

``` r
silence <- suppressPackageStartupMessages
silence(library(ggplot2))
silence(library(ggrepel))
silence(library(grid))
silence(library(dplyr))
silence(library(DESeq2))
```

    ## Warning: package 'DESeq2' was built under R version 4.0.3

    ## Warning: package 'S4Vectors' was built under R version 4.0.3

    ## Warning: package 'BiocGenerics' was built under R version 4.0.5

    ## Warning: package 'IRanges' was built under R version 4.0.3

    ## Warning: package 'GenomicRanges' was built under R version 4.0.3

    ## Warning: package 'GenomeInfoDb' was built under R version 4.0.5

    ## Warning: package 'SummarizedExperiment' was built under R version 4.0.3

    ## Warning: package 'MatrixGenerics' was built under R version 4.0.3

    ## Warning: package 'Biobase' was built under R version 4.0.3

``` r
silence(library(pheatmap))
silence(library(ggVennDiagram))
silence(library(UpSetR))
silence(library(Seurat))
silence(library(gtools))
silence(library(biomaRt))
```

    ## Warning: package 'biomaRt' was built under R version 4.0.3

    ## Possible Ensembl SSL connectivity problems detected.
    ## Please see the 'Connection Troubleshooting' section of the biomaRt vignette
    ## vignette('accessing_ensembl', package = 'biomaRt')Error in curl::curl_fetch_memory(url, handle = handle) : 
    ##   SSL peer certificate or SSH remote key was not OK: [uswest.ensembl.org] SSL certificate problem: certificate has expired

``` r
silence(library(ComplexHeatmap))
```

    ## Warning: package 'ComplexHeatmap' was built under R version 4.0.3

## GSE169255 RNA-seq analysis

``` r
# load count matrix
gse169255 <- read.table("GSE169255_sample_id_ReadsperGene.txt", sep = "\t", header = T)

row.names(gse169255) <- gse169255$gene_id

gse169255 <- gse169255[, -1]

# remove unpaired samples
gse169255 <- gse169255[, -c(1:6)]

# create coldata
gse169255_coldata <- data.frame(cbind(c(rep("myometrium", 12), rep("fibroid", 6)), c(rep("unpaired", 6), rep("paired", 12))))
gse169255_coldata <- filter(gse169255_coldata, X2 != "unpaired")
row.names(gse169255_coldata) <- colnames(gse169255)
colnames(gse169255_coldata) <- c("tissue", "paired_status")

# differential expression analysis
dds_gse169255 <- DESeqDataSetFromMatrix(countData = gse169255,
                              colData = gse169255_coldata,
                              design = ~ tissue)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
keep_gse169255 <- rowSums(counts(dds_gse169255)) >= 10
dds_gse169255 <- dds_gse169255[keep_gse169255,]


dds_gse169255$tissue <- factor(dds_gse169255$tissue, levels = c("myometrium","fibroid"))

dds_gse169255 <- DESeq(dds_gse169255)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
res_gse169255 <- results(dds_gse169255, pAdjustMethod = 'fdr')
vsd_gse169255 <- vst(dds_gse169255, blind = FALSE)

res_gse169255_df <- data.frame(res_gse169255)


# pick DEGs
gse169255_sig_up <- row.names(filter(res_gse169255_df, log2FoldChange > 0 & padj < 0.05))

gse169255_sig_down <- row.names(filter(res_gse169255_df, log2FoldChange < 0 & padj < 0.05))

gse169255_deg <- c(gse169255_sig_up, gse169255_sig_down)

myColor <- colorRampPalette(c('blue','white', 'red'))(100)

pheatmap(assay(vsd_gse169255)[gse169255_deg,],
                    cluster_rows = TRUE, 
                    show_rownames = FALSE,
                    cluster_cols = TRUE, 
                    annotation_col = gse169255_coldata, 
                    scale = 'row',
                    border_color = NA,
                    color = myColor)
```

    ## `use_raster` is automatically set to TRUE for a matrix with more than
    ## 2000 rows. You can control `use_raster` argument by explicitly setting
    ## TRUE/FALSE to it.
    ## 
    ## Set `ht_opt$message = FALSE` to turn off this message.

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Analyze Debu’s RNA-seq data

``` r
# load data
debu_rna_counts <- read.table("Debu_rna_seq_leio_myo_counts.txt", header = T)

row.names(debu_rna_counts) <- debu_rna_counts[, 1]
debu_rna_counts <- debu_rna_counts[, -1, FALSE]

# create coldata
debu_coldata <- data.frame(c(rep("myometrium", 15), rep("fibroid", 15)))
row.names(debu_coldata) <- colnames(debu_rna_counts)
colnames(debu_coldata) <- "tissue"

# DEseq
dds_debu <- DESeqDataSetFromMatrix(countData = debu_rna_counts,
                              colData = debu_coldata,
                              design = ~ tissue)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
keep <- rowSums(counts(dds_debu)) >= 10
dds_debu <- dds_debu[keep,]


dds_debu$tissue <- factor(dds_debu$tissue, levels = c("myometrium","fibroid"))

dds_debu <- DESeq(dds_debu)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 378 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
res_debu <- results(dds_debu, pAdjustMethod = 'fdr')
vsd_debu <- vst(dds_debu, blind = FALSE)

res_debu_df <- data.frame(res_debu)


# pick DEGs
debu_sig_up <- row.names(filter(res_debu_df, log2FoldChange > 0 & padj < 0.05))

debu_sig_down <- row.names(filter(res_debu_df, log2FoldChange < 0 & padj < 0.05))

debu_deg <- c(debu_sig_up, debu_sig_down)

pheatmap(assay(vsd_debu)[debu_deg,],
                    cluster_rows = TRUE, 
                    show_rownames = FALSE,
                    cluster_cols = TRUE, 
                    annotation_col = debu_coldata, 
                    scale = 'row',
                    border_color = NA,
                    color = myColor)
```

    ## `use_raster` is automatically set to TRUE for a matrix with more than
    ## 2000 rows. You can control `use_raster` argument by explicitly setting
    ## TRUE/FALSE to it.
    ## 
    ## Set `ht_opt$message = FALSE` to turn off this message.

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Read in GWAS candidate gene data

``` r
gwas <- read.table("Kadir_FUMA_job89348_050823/genes.txt", header = T)

gwas_genes <- unique(gwas$symbol)
gwas_ensembl <- unique(gwas$ensg)
```

## Find commonly up and down genes across RNA datasets and compare to GWAS genes

``` r
## all DEGs between both studies
all_degs <- unique(c(gse169255_sig_up, gse169255_sig_down, debu_sig_up, debu_sig_down))

## common genes between all DEGs and top GWAS genes (of 394)
all_degs_gwas <- unique(gwas[gwas$ensg %in% all_degs,]$symbol)

all_degs_gwas_ensg <- unique(gwas[gwas$ensg %in% all_degs,]$ensg)

## all up DEGs between both studies
all_up_degs <- unique(c(gse169255_sig_up, debu_sig_up))

## common genes between all up DEGs and top GWAS genes
all_up_degs_gwas <- unique(gwas[gwas$ensg %in% all_up_degs, ]$symbol)


## all down DEGs between both studies
all_down_degs <- unique(c(gse169255_sig_down, debu_sig_down))

## common genes between all up DEGs and top GWAS genes
all_down_degs_gwas <- unique(gwas[gwas$ensg %in% all_down_degs, ]$symbol)
```

## Venn diagram showing overlap of DEGs from both studies as well as GWAS genes

``` r
ggVennDiagram(list(A = c(debu_sig_up, debu_sig_down), B = c(gse169255_sig_up, gse169255_sig_down), C = gwas$ensg))
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Find pQTL-SNPs within risk loci or within 40 kb (+/- 20kb) of lead snps

``` r
# pQTLs within 20kb window of lead snps
pQTL_summary_stats <- read.csv("pQTL/Ferkingstad_41588_2021_978_MOESM4_ESM.txt", header = T, sep = "\t")

gallagher_genomic_loci_lead_snps <- read.table("Kadir_FUMA_job89348_050823/GenomicRiskLoci.txt", header = T)

gallagher_lead_snps_window_pqtl <- data.frame()
for (i in 1:nrow(gallagher_genomic_loci_lead_snps)){
  chr <- paste0("chr", gallagher_genomic_loci_lead_snps$chr[i])
  pos <- gallagher_genomic_loci_lead_snps$pos[i]
  lead_snp <- gallagher_genomic_loci_lead_snps$rsID[i]
  pqtl_data <- pQTL_summary_stats[pQTL_summary_stats$chr_var == chr, ]
  pqtl_data <- pqtl_data[abs(pos - pqtl_data$pos_var) < 20000, ]
  if (dim(pqtl_data)[1] == 0){
    next()
  } else {
    pqtl <- pqtl_data$gene_prot
    window_variant <- pqtl_data$variant
    beta <- pqtl_data$beta_adj
    lead_snp <- rep(lead_snp, length(pqtl))
    data <- data.frame(cbind(lead_snp, window_variant, pqtl, beta))
    gallagher_lead_snps_window_pqtl <- rbind(gallagher_lead_snps_window_pqtl, data)
  }
}

write.table(gallagher_lead_snps_window_pqtl, file = "Gallagher_Lead_SNP_20kb_window_pQTL.txt", sep = "\t", quote = F, row.names = F)

# draw graph showing number of pQTL snps within window of each lead snp
pqtl_snp_stats <- data.frame()
for (i in 1:nrow(gallagher_genomic_loci_lead_snps)){
  snp <- gallagher_genomic_loci_lead_snps$rsID[i]
  num <- nrow(gallagher_lead_snps_window_pqtl[gallagher_lead_snps_window_pqtl$lead_snp == snp, ])
  pqtl_snp_stats <- rbind(pqtl_snp_stats, c(snp, num))
}

colnames(pqtl_snp_stats) <- c("lead_snp", "pqtl_nums")
pqtl_snp_stats$lead_snp <- factor(pqtl_snp_stats$lead_snp, levels = rev(pqtl_snp_stats$lead_snp))
pqtl_snp_stats$pqtl_nums <- as.numeric(pqtl_snp_stats$pqtl_nums)

# Plot
ggplot(pqtl_snp_stats) +
  geom_col(aes(lead_snp, pqtl_nums)) + 
  scale_y_continuous(breaks=seq(0,12,by=2)) +
  coord_flip()
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Heatmap showing GWAS/RNA-seq common genes and their presence in certain datasets

``` r
# 0 means not present in data, 1 means present, 2 means upregulated, 3 means downregulated
genomic_position <- ifelse(gwas$posMapSNPs == 0, 0, 1)
eqtl <- ifelse(gwas$eqtlMapSNPs == 0, 0, 1)
chromatin_interaction <- ifelse(gwas$ciMap == "No", 0, 1)

debu_deg <- ifelse(gwas$ensg %in% debu_sig_down, 3, ifelse(gwas$ensg %in% debu_sig_up, 2, 0))
gse169255_deg <- ifelse(gwas$ensg %in% gse169255_sig_down, 3, ifelse(gwas$ensg %in% gse169255_sig_up, 2, 0))

gwas_matrix <- cbind(genomic_position, eqtl, chromatin_interaction, debu_deg, gse169255_deg)
gwas_matrix <- t(apply(gwas_matrix, MARGIN = 2, FUN = as.numeric))



colnames(gwas_matrix) <- gwas$symbol
gwas_matrix <- gwas_matrix[,duplicated(colnames(gwas_matrix)) == F]

gwas_matrix <- gwas_matrix[, colnames(gwas_matrix) %in% all_degs_gwas]

gwas_matrix_up <- gwas_matrix[, gwas_matrix[4,] == 2 | gwas_matrix[5,] == 2]

gwas_matrix_down <- gwas_matrix[, gwas_matrix[4,] == 3 | gwas_matrix[5,] == 3]


ComplexHeatmap::Heatmap(gwas_matrix_up, col = c("white", "green4", "red"), cluster_rows = F, rect_gp = gpar(col = "black"), show_heatmap_legend = F)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ComplexHeatmap::Heatmap(gwas_matrix_down, col = c("white", "green4", "blue"), cluster_rows = F, rect_gp = gpar(col = "black"), show_heatmap_legend = F)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

## Calculate significance of overlap between DEGs and FUMA-identified genes

``` r
# gse169255 number of candidate DEGS and non-DEGs
gse169255_candidate_DEG <- length(intersect(c(gse169255_sig_up, gse169255_sig_down), gwas$ensg))
gse169255_candidate_nonDEG <- length(gwas$ensg) - gse169255_candidate_DEG

# gse169255 number of non-candidate DEGS and non-DEGS
gse169255_nonCandidate_DEG <- length(c(gse169255_sig_up, gse169255_sig_down)) - gse169255_candidate_DEG
gse169255_nonCandidate_nonDEG <- length(row.names(gse169255)[-c(which(row.names(gse169255) %in% c(gse169255_sig_up, gse169255_sig_down, gwas$ensg)))])

gse169255_fuma_fisher <- data.frame(
  "Candidate" = c(gse169255_candidate_DEG, gse169255_candidate_nonDEG),
  "Non-Candidate" = c(gse169255_nonCandidate_DEG, gse169255_nonCandidate_nonDEG),
  row.names = c("DEG", "Unchanged"),
  stringsAsFactors = F
)

fisher.test(gse169255_fuma_fisher)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  gse169255_fuma_fisher
    ## p-value = 6.235e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  2.170934 3.445425
    ## sample estimates:
    ## odds ratio 
    ##   2.743401

``` r
# gse169255 number of candidate DEGS and non-DEGs
debu_candidate_DEG <- length(intersect(c(debu_sig_up, debu_sig_down), gwas$ensg))
debu_candidate_nonDEG <- length(gwas$ensg) - debu_candidate_DEG

# debu number of non-candidate DEGS and non-DEGS
debu_nonCandidate_DEG <- length(c(debu_sig_up, debu_sig_down)) - debu_candidate_DEG
debu_nonCandidate_nonDEG <- length(row.names(debu_rna_counts)[-c(which(row.names(debu_rna_counts) %in% c(debu_sig_up, debu_sig_down, gwas$ensg)))])

debu_fuma_fisher <- data.frame(
  "Candidate" = c(debu_candidate_DEG, debu_candidate_nonDEG),
  "Non-Candidate" = c(debu_nonCandidate_DEG, debu_nonCandidate_nonDEG),
  row.names = c("DEG", "Unchanged"),
  stringsAsFactors = F
)

fisher.test(debu_fuma_fisher)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  debu_fuma_fisher
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  3.019632 4.604204
    ## sample estimates:
    ## odds ratio 
    ##   3.734141

## Compare Gallagher, UK Biobank, and Japan Biobank studies

``` r
# look at risk loci overlaps
gallagher_risk_loci <- read.table("Kadir_FUMA_job89348_050823/GenomicRiskLoci.txt", header = T)
gallagher_risk_loci$dataset <- rep("gallagher", n = nrow(gallagher_risk_loci))

japan_biobank_risk_loci <- read.table("Japan_Biobank/FUMA_job259732/GenomicRiskLoci.txt", header = T)
japan_biobank_risk_loci$dataset <- rep("japan", n = nrow(japan_biobank_risk_loci))

uk_biobank_risk_loci <- read.table("UKB/ukbio_round2_FUMA_job253818/GenomicRiskLoci.txt", header = T)
uk_biobank_risk_loci$dataset <- rep("uk", n = nrow(uk_biobank_risk_loci))

# the variable window is defined as half the size that we want
window <- 10000

# japan/gallagher
japan_gallagher_overlap_risk_loci <- data.frame()
for (i in 1:23){
  chr <- i
  japan_data <- japan_biobank_risk_loci[japan_biobank_risk_loci$chr == chr,]
  gallagher_data <- gallagher_risk_loci[gallagher_risk_loci$chr == chr,]
  if (is.null(dim(japan_data)) | is.null(dim(gallagher_data))){
    next()
  }
  for (j in 1:nrow(japan_data)){
    locus <- japan_data[j,]
    locus_start <- locus$start
    locus_end <- locus$end
    start_abs <- abs(locus_start - gallagher_data$start)
    end_abs <- abs(locus_end - gallagher_data$end)
    gallagher_overlap_start <- gallagher_data[which(start_abs < window | end_abs < window),]$start
    gallagher_overlap_end <- gallagher_data[which(start_abs < window | end_abs < window),]$end
    if (length(gallagher_overlap_start) == 0 | length(gallagher_overlap_end) == 0){
      next()
    }
    japan_gallagher_overlap_risk_loci <- rbind(japan_gallagher_overlap_risk_loci, c(paste0("chr", chr), locus$start, locus$end, gallagher_overlap_start, gallagher_overlap_end))
  }
}

colnames(japan_gallagher_overlap_risk_loci) <- c("chr", "japan_biobank_start", "japan_biobank_end", "gallagher_start", "gallagher_end")

write.table(japan_gallagher_overlap_risk_loci, file = "japan_biobank_gallagher_gwas_risk_loci_overlap.txt", row.names = F, quote = F)


# uk/gallagher
uk_gallagher_overlap_risk_loci <- data.frame()
for (i in 1:23){
  chr <- i
  uk_data <- uk_biobank_risk_loci[uk_biobank_risk_loci$chr == chr,]
  gallagher_data <- gallagher_risk_loci[gallagher_risk_loci$chr == chr,]
  if (is.null(dim(uk_data)) | is.null(dim(gallagher_data))){
    next()
  }
  for (j in 1:nrow(uk_data)){
    locus <- uk_data[j,]
    locus_start <- locus$start
    locus_end <- locus$end
    start_abs <- abs(locus_start - gallagher_data$start)
    end_abs <- abs(locus_end - gallagher_data$end)
    gallagher_overlap_start <- gallagher_data[which(start_abs < window | end_abs < window),]$start
    gallagher_overlap_end <- gallagher_data[which(start_abs < window | end_abs < window),]$end
    if (length(gallagher_overlap_start) == 0 | length(gallagher_overlap_end) == 0){
      next()
    }
    uk_gallagher_overlap_risk_loci <- rbind(uk_gallagher_overlap_risk_loci, c(paste0("chr", chr), locus$start, locus$end, gallagher_overlap_start, gallagher_overlap_end))
  }
}

colnames(uk_gallagher_overlap_risk_loci) <- c("chr", "uk_biobank_start", "uk_biobank_end", "gallagher_start", "gallagher_end")

write.table(uk_gallagher_overlap_risk_loci, file = "uk_biobank_gallagher_gwas_risk_loci_overlap.txt", row.names = F, quote = F)


#uk/japan
uk_japan_overlap_risk_loci <- data.frame()
for (i in 1:23){
  chr <- i
  uk_data <- uk_biobank_risk_loci[uk_biobank_risk_loci$chr == chr,]
  japan_data <- japan_biobank_risk_loci[japan_biobank_risk_loci$chr == chr,]
  if (is.null(dim(uk_data)) | is.null(dim(japan_data))){
    next()
  }
  for (j in 1:nrow(uk_data)){
    locus <- uk_data[j,]
    locus_start <- locus$start
    locus_end <- locus$end
    start_abs <- abs(locus_start - japan_data$start)
    end_abs <- abs(locus_end - japan_data$end)
    japan_overlap_start <- japan_data[which(start_abs < window | end_abs < window),]$start
    japan_overlap_end <- japan_data[which(start_abs < window | end_abs < window),]$end
    if (length(japan_overlap_start) == 0 | length(japan_overlap_end) == 0){
      next()
    }
    uk_japan_overlap_risk_loci <- rbind(uk_japan_overlap_risk_loci, c(paste0("chr", chr), locus$start, locus$end, japan_overlap_start, japan_overlap_end))
  }
}

colnames(uk_japan_overlap_risk_loci) <- c("chr", "uk_biobank_start", "uk_biobank_end", "japan_start", "japan_end")

write.table(uk_japan_overlap_risk_loci, file = "uk_biobank_japan_biobank_gwas_risk_loci_overlap.txt", row.names = F, quote = F)


# only one genomic loci overlaps all 3 datasets (chr13:40723944), the rest only overlap 2 datasets, remove all overlapping data (making sure not to remove twice)
japan_only_loci <- paste0("japan_only_", seq_along(1:(nrow(japan_biobank_risk_loci) - nrow(japan_gallagher_overlap_risk_loci) - nrow(uk_japan_overlap_risk_loci) + 1)))

gallagher_only_loci <- paste0("gallagher_only_", seq_along(1:(nrow(gallagher_risk_loci) - nrow(japan_gallagher_overlap_risk_loci) - nrow(uk_gallagher_overlap_risk_loci) + 1)))

uk_only_loci <-  paste0("uk_only_", seq_along((1:(nrow(uk_biobank_risk_loci) - nrow(uk_japan_overlap_risk_loci) - nrow(uk_gallagher_overlap_risk_loci) + 1))))

japan_gallagher_loci <- paste0("japan_gallagher_overlap_", seq_along(1:(nrow(japan_gallagher_overlap_risk_loci)-1)))

uk_gallagher_loci <- paste0("uk_gallagher_overlap_", seq_along(1:(nrow(uk_gallagher_overlap_risk_loci)-1)))

japan_uk_gallagher_loci <- "uk_gallagher_japan_overlap_1"

all_genomic_risk_loci_vector <- c(japan_only_loci, gallagher_only_loci, uk_only_loci, japan_gallagher_loci, uk_gallagher_loci, japan_uk_gallagher_loci)

# create binary matrix for upset plot
all_genomic_risk_loci_df <- data.frame(matrix(ncol = 3, nrow = length(all_genomic_risk_loci_vector)))
for (i in 1:length(all_genomic_risk_loci_vector)){
  variant <- all_genomic_risk_loci_vector[i]
  data <- ifelse(rep(grepl(pattern = "gallagher_only", x = variant), 3), c(1, 0, 0),
                 ifelse(rep(grepl(pattern = "uk_only", x = variant), 3), c(0, 1, 0),
                        ifelse(rep(grepl(pattern = "japan_only", x = variant), 3), c(0, 0, 1),
                               ifelse(rep(grepl(pattern = "japan_gallagher_overlap", x = variant), 3), c(1, 0, 1),
                                      ifelse(rep(grepl(pattern = "uk_gallagher_overlap", x = variant), 3), c(1, 1, 0),
                                             c(1, 1, 1))))))
  all_genomic_risk_loci_df[i,] <- as.numeric(data)
}

colnames(all_genomic_risk_loci_df) <- c("gallagher", "uk_biobank", "japan_biobank")


# Upset plot!
upset(data.frame(all_genomic_risk_loci_df), sets = c("japan_biobank", "uk_biobank", "gallagher"), number.angles = 30, point.size = 3.5, line.size = 2, 
    mainbar.y.label = "Genomic Risk Loci Intersections", sets.x.label = "Genomic Risk Loci Per Dataset",
    text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), keep.order = T)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# look at candidate SNP overlaps
japan_biobank_candidate_snps <- read.table("Japan_Biobank/FUMA_job259732/snps.txt", header = T)

gallagher_candidate_snps <- read.table("Kadir_FUMA_job89348_050823/snps.txt", header = T)

japan_biobank_gallagher_snp_overlap <- merge(japan_biobank_candidate_snps, gallagher_candidate_snps, by = "rsID")

write.table(japan_biobank_gallagher_snp_overlap, file = "japan_biobank_gallagher_gwas_candidate_snp_overlap.txt", row.names = F, quote = F)
```

# Single cell RNA-seq data analysis

## Load single cell RNA-seq data

``` r
dirs <- list.dirs(path = "GSE162122_RAW")[-1]

data_list <- lapply(dirs, FUN = function(x){
  Read10X(x)
})

names(data_list) <- sapply(dirs, FUN = function(x){
  sub('.*/', '', x)
})

# create seurat object lists by tissue type
seurat_myo_list <- lapply(data_list[grepl(pattern = "Myo", names(data_list)) == T], FUN = function(x){
  CreateSeuratObject(counts = x, project = "Myo_Lyo", min.cells = 3, min.features = 200)
})
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
seurat_fib_list <- lapply(data_list[grepl(pattern = "Fibroid", names(data_list)) == T], FUN = function(x){
  CreateSeuratObject(counts = x, project = "Myo_Lyo", min.cells = 3, min.features = 200)
})
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
# add patient ids
for (i in seq_along(seurat_myo_list)){
  seurat_myo_list[[i]][["patient.id"]] <- gsub('_.*', '', names(seurat_myo_list[i]))
}

for (i in seq_along(seurat_fib_list)){
  seurat_fib_list[[i]][["patient.id"]] <- gsub('_.*', '', names(seurat_fib_list[i]))
}
```

## Perform QC on scRNA-seq data

## Identify percent MT, percent hemoglobin, and percent EPCAM/KRT18

``` r
# calculate percent MT
for (i in seq_along(seurat_myo_list)){
  seurat_myo_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_myo_list[[i]], pattern = "^MT-")
}

for (i in seq_along(seurat_fib_list)){
  seurat_fib_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_fib_list[[i]], pattern = "^MT-")
}

# calculate percent hemoglobin
for (i in seq_along(seurat_myo_list)){
  seurat_myo_list[[i]][["percent.hemo"]] <- tryCatch({PercentageFeatureSet(seurat_myo_list[[i]], features = c("HBA1", "HBA2", "HBB"))},
           error = function(e){
             return(0)
           })
}

for (i in seq_along(seurat_fib_list)){
  seurat_fib_list[[i]][["percent.hemo"]] <- tryCatch({PercentageFeatureSet(seurat_fib_list[[i]], features = c("HBA1", "HBA2", "HBB"))},
           error = function(e){
             return(0)
           })
}

# calculate percent EPCAM/KRT18
for (i in seq_along(seurat_myo_list)){
  seurat_myo_list[[i]][["percent.endometrioid"]] <- tryCatch({PercentageFeatureSet(seurat_myo_list[[i]], features = c("EPCAM", "KRT18"))},
           error = function(e){
             return(0)
           })
}

for (i in seq_along(seurat_fib_list)){
  seurat_fib_list[[i]][["percent.endometrioid"]] <- tryCatch({PercentageFeatureSet(seurat_fib_list[[i]], features = c("EPCAM", "KRT18"))},
           error = function(e){
             return(0)
           })
}

lapply(seurat_myo_list, FUN = function(x){
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hemo", "percent.endometrioid"), ncol = 3)
})
```

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent.endometrioid.

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent.endometrioid.

    ## $`GSM4942396_Myometrium-55_12640`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## 
    ## $`GSM4942397_Myometrium-55_11564`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

    ## 
    ## $`GSM4942398_Myometrium-55.11911`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

    ## 
    ## $GSM4942399_061919Myomterium

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

    ## 
    ## $`GSM5023319_Myometrium-55_12745`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

``` r
lapply(seurat_fib_list, FUN = function(x){
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hemo", "percent.endometrioid"), ncol = 3)
})
```

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent.endometrioid.

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent.endometrioid.

    ## $`GSM5023320_Fibroid55-12906-1`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->

    ## 
    ## $`GSM5023320_Fibroid55-12906-2`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-7.png)<!-- -->

    ## 
    ## $`GSM5023321_Fibroid-55_12382`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-8.png)<!-- -->

    ## 
    ## $`GSM5023323_Fibroid-55_12640`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-9.png)<!-- -->

    ## 
    ## $`GSM5023324_Fibroid-55_12843_1_GEXL`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-10.png)<!-- -->

    ## 
    ## $`GSM5023324_Fibroid-55_12843_2_GEXL`

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-11.png)<!-- -->

    ## 
    ## $GSM6509142_Fibroid061919

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-12-12.png)<!-- -->

## Subset based on QC values

``` r
# QC subsetting
seurat_myo_list <- lapply(seurat_myo_list, FUN = function(x){
  subset(x, subset = nFeature_RNA < 2500 & percent.mt < 7 & percent.hemo < 1 & percent.endometrioid == 0)
})

seurat_fib_list <- lapply(seurat_fib_list, FUN = function(x){
  subset(x, subset = nFeature_RNA < 2500 & percent.mt < 7 & percent.hemo < 1 & percent.endometrioid == 0)
})
```

## Normalize data and select integration features

``` r
seurat_myo_list <- suppressMessages(lapply(seurat_myo_list, function(x){
  SCTransform(x)
}))
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  21%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  79%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  21%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  79%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  21%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  79%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  21%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  79%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

``` r
seurat_fib_list <- suppressMessages(lapply(seurat_fib_list, function(x){
  SCTransform(x)
}))
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |=======                                                               |   9%  |                                                                              |=========                                                             |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  19%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  59%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  91%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |=======                                                               |   9%  |                                                                              |=========                                                             |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  19%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  59%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  91%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   9%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  79%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  91%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   9%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  79%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  91%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |============================================================          |  85%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |============================================================          |  85%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  38%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  58%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |=================================================================     |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  38%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  58%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |=================================================================     |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  81%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  81%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

``` r
# select features that are repeatedly variable across datasets for integration
features_myo <- SelectIntegrationFeatures(object.list = seurat_myo_list)

features_fib <- SelectIntegrationFeatures(object.list = seurat_fib_list)
```

## Perform integration

``` r
anchors_myo <- suppressMessages(FindIntegrationAnchors(object.list = seurat_myo_list, anchor.features = features_myo))
```

    ## Warning in CheckDuplicateCellNames(object.list = object.list): Some cell names
    ## are duplicated across objects provided. Renaming to enforce unique cell names.

``` r
anchors_fib <- suppressMessages(FindIntegrationAnchors(object.list = seurat_fib_list, anchor.features = features_fib))
```

    ## Warning in CheckDuplicateCellNames(object.list = object.list): Some cell names
    ## are duplicated across objects provided. Renaming to enforce unique cell names.

``` r
# this command creates an 'integrated' data assay
seurat_combined_myo <- suppressMessages(IntegrateData(anchorset = anchors_myo))

seurat_combined_fib <- suppressMessages(IntegrateData(anchorset = anchors_fib))
```

    ## Warning: Attempting to merge an SCTAssay with another Assay type 
    ## Converting all to standard Assay objects.

    ## Warning: Attempting to merge an SCTAssay with another Assay type 
    ## Converting all to standard Assay objects.

## Run PCA

``` r
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat_combined_myo) <- "integrated"

DefaultAssay(seurat_combined_fib) <- "integrated"


# SCTransform integrated data and run PCA
seurat_combined_myo <- SCTransform(seurat_combined_myo)
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 18117 by 6394

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 88 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 18117 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  19%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  41%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  59%  |                                                                              |============================================                          |  62%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 18117 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  19%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  41%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  59%  |                                                                              |============================================                          |  62%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 47.37542 secs

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Centering data matrix

    ## Set default assay to SCT

``` r
seurat_combined_myo <- RunPCA(seurat_combined_myo, npcs = 30, verbose = FALSE)

seurat_combined_fib <- SCTransform(seurat_combined_fib)
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 19737 by 41948

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 92 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 19737 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 19737 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 2.911297 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Centering data matrix

    ## Set default assay to SCT

``` r
seurat_combined_fib <- RunPCA(seurat_combined_fib, npcs = 30, verbose = FALSE)
```

## Create elbow plots to determine number of PCs

``` r
ElbowPlot(seurat_combined_myo)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ElbowPlot(seurat_combined_fib)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

## Perform UMAP and clustering

``` r
seurat_combined_myo <- RunUMAP(seurat_combined_myo, reduction = "pca", dims = 1:20)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 14:18:03 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:18:03 Read 6394 rows and found 20 numeric columns

    ## 14:18:03 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:18:03 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:18:04 Writing NN index file to temp file /var/folders/q_/mgl3_x0s1rn09zjqyc77n8zc0000gp/T//RtmpKyig2p/file44681eda04f9
    ## 14:18:04 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:18:05 Annoy recall = 100%
    ## 14:18:06 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:18:07 Initializing from normalized Laplacian + noise (using irlba)
    ## 14:18:07 Commencing optimization for 500 epochs, with 258620 positive edges
    ## 14:18:13 Optimization finished

``` r
seurat_combined_myo <- FindNeighbors(seurat_combined_myo, reduction = "pca", dims = 1:20)
```

    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
seurat_combined_myo <- FindClusters(seurat_combined_myo, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 6394
    ## Number of edges: 206658
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9393
    ## Number of communities: 19
    ## Elapsed time: 0 seconds

``` r
seurat_combined_fib <- RunUMAP(seurat_combined_fib, reduction = "pca", dims = 1:20)
```

    ## 14:18:15 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 14:18:15 Read 41948 rows and found 20 numeric columns
    ## 14:18:15 Using Annoy for neighbor search, n_neighbors = 30
    ## 14:18:15 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 14:18:19 Writing NN index file to temp file /var/folders/q_/mgl3_x0s1rn09zjqyc77n8zc0000gp/T//RtmpKyig2p/file44687b7c8741
    ## 14:18:19 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:18:28 Annoy recall = 100%
    ## 14:18:29 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:18:30 Initializing from normalized Laplacian + noise (using irlba)
    ## 14:18:32 Commencing optimization for 200 epochs, with 1924796 positive edges
    ## 14:18:50 Optimization finished

``` r
seurat_combined_fib <- FindNeighbors(seurat_combined_fib, reduction = "pca", dims = 1:20)
```

    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
seurat_combined_fib <- FindClusters(seurat_combined_fib, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 41948
    ## Number of edges: 1333287
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9040
    ## Number of communities: 18
    ## Elapsed time: 7 seconds

# Visualization

``` r
DimPlot(seurat_combined_myo, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Normal Myometrium scRNA-seq")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
DimPlot(seurat_combined_fib, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Uterine Fibroid scRNA-seq")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
DimPlot(seurat_combined_myo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "patient.id") + ggtitle("Normal Myometrium scRNA-seq")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

``` r
DimPlot(seurat_combined_fib, reduction = "umap", label = TRUE, repel = TRUE, group.by = "patient.id") + ggtitle("Uterine Fibroid scRNA-seq")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->

## Find cluster markers

``` r
seurat_combined_myo <- PrepSCTFindMarkers(seurat_combined_myo, assay = "SCT")
```

    ## Only one SCT model is stored - skipping recalculating corrected counts

``` r
myo_de_genes <- FindAllMarkers(seurat_combined_myo, assay = "SCT", only.pos = T)
```

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

    ## Calculating cluster 9

    ## Calculating cluster 10

    ## Calculating cluster 11

    ## Calculating cluster 12

    ## Calculating cluster 13

    ## Calculating cluster 14

    ## Calculating cluster 15

    ## Calculating cluster 16

    ## Calculating cluster 17

    ## Calculating cluster 18

``` r
write.table(myo_de_genes, file = "myo_de_genes_scrnaseq.txt", quote = F, sep = "\t")

seurat_combined_fib <- PrepSCTFindMarkers(seurat_combined_fib, assay = "SCT")
```

    ## Only one SCT model is stored - skipping recalculating corrected counts

``` r
fib_de_genes <- FindAllMarkers(seurat_combined_fib, assay = "SCT", only.pos = T)
```

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

    ## Calculating cluster 9

    ## Calculating cluster 10

    ## Calculating cluster 11

    ## Calculating cluster 12

    ## Calculating cluster 13

    ## Calculating cluster 14

    ## Calculating cluster 15

    ## Calculating cluster 16

    ## Calculating cluster 17

``` r
write.table(fib_de_genes, file = "fib_de_genes_scrnaseq.txt", quote = F, sep = "\t")
```

## Cluster 17 in the Leiomyoma sample represents leftover endometrioid contamination. Remove it and re-cluster the sample

``` r
seurat_combined_fib <- subset(seurat_combined_fib, subset = seurat_clusters != "17")

# re-analyze
seurat_combined_fib <- RunUMAP(seurat_combined_fib, reduction = "pca", dims = 1:20)
```

    ## 14:21:23 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:21:23 Read 41864 rows and found 20 numeric columns

    ## 14:21:23 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:21:23 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:21:27 Writing NN index file to temp file /var/folders/q_/mgl3_x0s1rn09zjqyc77n8zc0000gp/T//RtmpKyig2p/file44682c5d53c7
    ## 14:21:27 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:21:37 Annoy recall = 100%
    ## 14:21:38 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:21:40 Initializing from normalized Laplacian + noise (using irlba)
    ## 14:21:42 Commencing optimization for 200 epochs, with 1921714 positive edges
    ## 14:21:59 Optimization finished

``` r
seurat_combined_fib <- FindNeighbors(seurat_combined_fib, reduction = "pca", dims = 1:20)
```

    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
seurat_combined_fib <- FindClusters(seurat_combined_fib, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 41864
    ## Number of edges: 1331467
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9042
    ## Number of communities: 18
    ## Elapsed time: 8 seconds

``` r
# visualize
DimPlot(seurat_combined_fib, label = T)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# re-perform DE analysis
fib_de_genes <- FindAllMarkers(seurat_combined_fib, assay = "SCT", only.pos = T)
```

    ## Calculating cluster 0
    ## Calculating cluster 1
    ## Calculating cluster 2
    ## Calculating cluster 3
    ## Calculating cluster 4
    ## Calculating cluster 5
    ## Calculating cluster 6
    ## Calculating cluster 7
    ## Calculating cluster 8
    ## Calculating cluster 9
    ## Calculating cluster 10
    ## Calculating cluster 11
    ## Calculating cluster 12
    ## Calculating cluster 13
    ## Calculating cluster 14
    ## Calculating cluster 15
    ## Calculating cluster 16
    ## Calculating cluster 17

``` r
write.table(fib_de_genes, file = "fib_de_genes_scrnaseq.txt", quote = F, sep = "\t")
```

## Plot marker genes using FeaturePlot

``` r
FeaturePlot(seurat_combined_myo, features = c("MYH11", "DCN", "PECAM1", "CD3D", "CD14"))
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
FeaturePlot(seurat_combined_fib, features = c("MYH11", "DCN", "PECAM1", "CD3D", "CD14"))
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

## Identify cell cluster types in individual seurat objects and plot cell frequencies

``` r
seurat_combined_myo[["cell.type"]] <- sapply(seurat_combined_myo$seurat_clusters, FUN = function(x){
  cell.type <- ifelse(x %in% c(0), "smc",
                      ifelse(x %in% c(7), "fibroblast",
                             ifelse(x %in% c(2, 4, 6, 9, 11, 12, 14, 15, 16), "endothelial",
                                    ifelse(x %in% c(1, 3, 8, 10, 18), "lymphoid",
                                           ifelse(x %in% c(5, 13, 17), "myeloid",
                                                  NA)))))
  return(cell.type)
})

seurat_combined_myo[["cell.type"]] <- factor(seurat_combined_myo$cell.type, levels = c("smc", "fibroblast", "endothelial", "lymphoid", "myeloid"))

seurat_combined_fib[["cell.type"]] <- sapply(seurat_combined_fib$seurat_clusters, FUN = function(x){
  cell.type <- ifelse(x %in% c(0, 1, 2, 4, 5, 13, 14), "smc",
                      ifelse(x %in% c(17), "fibroblast",
                             ifelse(x %in% c(7, 9, 10, 11), "endothelial",
                                    ifelse(x %in% c(3, 6, 8, 15), "lymphoid",
                                           ifelse(x %in% c(12, 16), "myeloid",
                                                  NA)))))
  return(cell.type)
})

seurat_combined_fib[["cell.type"]] <- factor(seurat_combined_fib$cell.type, levels = c("smc", "fibroblast", "endothelial", "lymphoid", "myeloid"))

# Plot UMAPs
DimPlot(seurat_combined_myo, group.by = "cell.type")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
DimPlot(seurat_combined_fib, group.by = "cell.type")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

``` r
# make bar charts showing distribution of each cell type
bar_data_myo <- as.data.frame(t(sapply(levels(factor(seurat_combined_myo$cell.type)), FUN = function(x){
  num <- length(seurat_combined_myo$cell.type[which(seurat_combined_myo$cell.type == x)])
  return(c(x,num))
})))

colnames(bar_data_myo) <- c("cell.type", "count")
bar_data_myo$count <- as.numeric(bar_data_myo$count)
bar_data_myo <- mutate(bar_data_myo, percent = (count / sum(bar_data_myo$count)) * 100)
bar_data_myo <- arrange(bar_data_myo, desc(percent))
bar_data_myo$cell.type <- factor(bar_data_myo$cell.type, levels = bar_data_myo$cell.type)

bar_data_fib <- as.data.frame(t(sapply(levels(factor(seurat_combined_fib$cell.type)), FUN = function(x){
  num <- length(seurat_combined_fib$cell.type[which(seurat_combined_fib$cell.type == x)])
  return(c(x,num))
})))

colnames(bar_data_fib) <- c("cell.type", "count")
bar_data_fib$count <- as.numeric(bar_data_fib$count)
bar_data_fib <- mutate(bar_data_fib, percent = (count / sum(bar_data_fib$count)) * 100)
bar_data_fib <- arrange(bar_data_fib, desc(percent))
bar_data_fib$cell.type <- factor(bar_data_fib$cell.type, levels = bar_data_fib$cell.type)


ggplot(bar_data_myo, aes(x = cell.type, y = percent, fill = cell.type)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 45))
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->

``` r
ggplot(bar_data_fib, aes(x = cell.type, y = percent, fill = cell.type)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 45))
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->

``` r
# combine bar data for grouped bar graph
bar_data_myo$tissue <- "myometrium"
bar_data_fib$tissue <- "leiomyoma"

bar_data_all <- rbind(bar_data_myo, bar_data_fib)


bar_data_all$tissue <- factor(bar_data_all$tissue, levels = c("myometrium", "leiomyoma"))
bar_data_all$percent <- as.numeric(bar_data_all$percent)
bar_data_all$cell.type <- factor(bar_data_all$cell.type, levels = c("smc", "endothelial", "lymphoid", "myeloid", "fibroblast"))

ggplot(bar_data_all, aes(x = cell.type, y = percent, fill = tissue)) +
  geom_bar(position = "dodge", stat = "identity")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->

## Feature Plots in Leiomyoma

``` r
FeaturePlot(seurat_combined_fib, features = c("KANK1", "VCAN", "APOBEC3C", "C1QC"), order = T)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

## Save seurat objects

``` r
saveRDS(seurat_combined_myo, file = "GSE162122_scRNA_Myometrium_Combined.RDS")
saveRDS(seurat_combined_fib, file = "GSE162122_scRNA_Fibroid_Combined.RDS")

seurat_combined_myo <- readRDS("GSE162122_scRNA_Myometrium_Combined.RDS")
seurat_combined_fib <- readRDS("GSE162122_scRNA_Fibroid_Combined.RDS")
```

## Plot heatmaps of pseudobulk gene expression of genes identified as differentially expressed between leiomyoma and myometrium

``` r
# Leiomyoma
fib_up_degs_gwas <- FetchData(seurat_combined_fib, vars = c("cell.type", all_up_degs_gwas))
```

    ## Warning: Could not find KCNF1 in the default search locations, found in RNA
    ## assay instead

    ## Warning in FetchData.Seurat(seurat_combined_fib, vars = c("cell.type",
    ## all_up_degs_gwas)): The following requested variables were not found: UGT2B11,
    ## UGT2A1, C11orf87, C17orf61-PLSCR3

``` r
# create pseudobulk data
fib_up_degs_gwas_pseudo <- data.frame()
for (i in 1:length(levels(factor(fib_up_degs_gwas$cell.type)))){
  cell_type <- levels(factor(fib_up_degs_gwas$cell.type))[i]
  data <- fib_up_degs_gwas[fib_up_degs_gwas$cell.type == cell_type, -1]
  data <- colMeans(data)
  fib_up_degs_gwas_pseudo <- rbind(fib_up_degs_gwas_pseudo, data)
}

row.names(fib_up_degs_gwas_pseudo) <- levels(factor(fib_up_degs_gwas$cell.type))
fib_up_degs_gwas_pseudo <- apply(fib_up_degs_gwas_pseudo, MARGIN = 1, as.numeric)
row.names(fib_up_degs_gwas_pseudo) <- colnames(fib_up_degs_gwas[, -1])

# plot heatmaps
gene_order <- c("VMP1", "SH3PXD2A", "ELMOD1", "PDCD11", "NEURL4", "HSPG2", "SLC38A4", "CNTROB", "GPR162", "SORCS3", "PKP3", "KANK1", "B4GALNT4", "NKD1", "KNTC1", "PTDSS2", "MGAT3", "MLXIP", "KCNAB3", "CHD3", "USP46", "KIF5C", "DVL2", "KDM1A", "SALL1", "CDCA7", "CACNA1E", "RILPL1", "NPTXR", "ATP6V0A2", "WT1", "DACT1", "RAP1GAP", "ADAM23", "ERVMER34-1", "RASL11B", "VCAN", "DEPDC7", "WNT4", "EPHB2", "ELOVL4", "DNAH10", "FAM101A", "HIST1H4H", "SULT1E1", "AEN", "APOBEC3F", "DNAH2", "LRRC56", "SHMT2", "rna_KCNF1", "TP53", "MPPED2", "RGS17", "PDIA6", "SLC25A15", "TMEM256-PLSCR3", "RBX1", "SEPT4", "CDH2", "TMEM256", "SNRNP35", "TRIP13", "NBPF3", "HIST1H2BK", "APOBEC3C", "GLIPR2", "SLC38A1", "WRAP53", "HIST1H2BD", "DDX55", "AURKB", "F13A1", "C1QC", "C1QA", "E2F2", "SLC24A4", "GPR160", "TOP2A", "TMEM102", "PFAS", "KDM2B")

cell_type_order <- c("smc", "fibroblast", "endothelial", "lymphoid", "myeloid")

fib_up_degs_gwas_pseudo <- fib_up_degs_gwas_pseudo[gene_order, cell_type_order]

pheatmap(fib_up_degs_gwas_pseudo, scale = "row", cluster_rows = F, cluster_cols = F)
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

## Plot subset of FUMA GWAS targets along with percentage of cells expressing

``` r
fib_up_degs_gwas_pseudo_subset <- data.frame()
for (i in 1:length(levels(factor(fib_up_degs_gwas$cell.type)))){
  cell_type <- levels(factor(fib_up_degs_gwas$cell.type))[i]
  data <- fib_up_degs_gwas[fib_up_degs_gwas$cell.type == cell_type, -1]
  means <- colMeans(data)
  percents <- apply(data, 2, function(x){
    (sum(x > 0)/length(x)) * 100
  })
  all_data <- data.frame(cbind(rep(cell_type, n = length(means)), means, percents))
  colnames(all_data) <- c("cell_type", "mean", "percent")
  all_data$gene <- row.names(all_data)
  row.names(all_data) <- NULL
  fib_up_degs_gwas_pseudo_subset <- rbind(fib_up_degs_gwas_pseudo_subset, all_data)
}

fib_up_degs_gwas_pseudo_subset$cell_type <- factor(fib_up_degs_gwas_pseudo_subset$cell_type, levels = c("smc", "fibroblast", "endothelial", "lymphoid", "myeloid"))

fib_up_degs_gwas_pseudo_subset$mean <- as.numeric(fib_up_degs_gwas_pseudo_subset$mean)

fib_up_degs_gwas_pseudo_subset$percent <- as.numeric(fib_up_degs_gwas_pseudo_subset$percent)

# filter for best marker genes
fib_up_degs_gwas_pseudo_subset <- fib_up_degs_gwas_pseudo_subset[fib_up_degs_gwas_pseudo_subset$gene %in% c("HSPG2", "VMP1", "CHD3", "CNTROB", "KANK1", "VCAN", "RBX1", "PDIA6", "C1QA", "C1QC", "F13A1"), ]

# calculate z-score
fib_up_degs_gwas_pseudo_subset_zscores <- data.frame()
for(i in 1:length(levels(factor(fib_up_degs_gwas_pseudo_subset$gene)))){
  gene <- levels(factor(fib_up_degs_gwas_pseudo_subset$gene))[i]
  data <- fib_up_degs_gwas_pseudo_subset[fib_up_degs_gwas_pseudo_subset$gene == gene, ]
  data$zscore <- sapply(data$mean, function(x){
    (x-mean(data$mean))/sd(data$mean)
  })
  fib_up_degs_gwas_pseudo_subset_zscores <- rbind(fib_up_degs_gwas_pseudo_subset_zscores, data)
}

fib_up_degs_gwas_pseudo_subset_zscores$zscore <- as.numeric(fib_up_degs_gwas_pseudo_subset_zscores$zscore)
fib_up_degs_gwas_pseudo_subset_zscores$gene <- factor(fib_up_degs_gwas_pseudo_subset_zscores$gene, levels = rev(c("HSPG2", "VMP1", "CHD3", "CNTROB", "KANK1", "VCAN", "RBX1", "PDIA6", "C1QA", "C1QC", "F13A1")))

# visualize
ggplot(fib_up_degs_gwas_pseudo_subset_zscores, aes(x = cell_type, y = gene, color = zscore, size = percent)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
```

![](Fibroid-GWAS-Manuscript_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->
