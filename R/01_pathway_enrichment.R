R.version.string        # prints the installed R version
# BiocManager is the installer for Bioconductor packages (clusterProfiler is one of them)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# clusterProfiler: enrichment analysis; org.Hs.eg.db: human gene annotation; enrichplot: plots
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), ask = FALSE, update = FALSE)
# readr reads CSV files, dplyr filters and groups tables, ggplot2 draws plots
install.packages(c("readr", "dplyr", "ggplot2"))
library(clusterProfiler); library(org.Hs.eg.db); library(readr); library(dplyr); library(ggplot2)

project_dir  <- "C:/Users/ADMIN/Documents/human-blastoid-scrna"          # forward slashes in R paths on Windows
markers_file <- file.path(project_dir, "Result", "markers_leiden_r05.csv")
stopifnot(file.exists(markers_file))

markers <- read_csv(markers_file)      # one row per gene per cluster
glimpse(markers)                       # shows the columns so you can check their names

padj_cut <- 0.05; lfc_cut <- 0.5; pct_cut <- 0.25; n_top <- 200     # selection thresholds (my choices)
markers <- markers %>% mutate(group = factor(group, levels = 0:13)) # cluster number as an ordered factor
sig <- markers %>%
  filter(pvals_adj < padj_cut, logfoldchanges > lfc_cut, pct_nz_group > pct_cut) %>%
  group_by(group) %>% slice_max(scores, n = n_top) %>% ungroup()
print(table(sig$group))                # number of genes selected for each cluster

gene_list <- split(sig$names, sig$group)   # named list: one vector of gene symbols per cluster
universe  <- unique(markers$names)         # background: every gene tested in the marker analysis
mapped <- bitr(unique(sig$names), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("mapped:", nrow(mapped), "of", length(unique(sig$names)), "symbols\n")
ego <- compareCluster(
  geneClusters  = gene_list,           # list of gene vectors (one per cluster)
  fun           = "enrichGO",          # Gene Ontology over-representation test
  OrgDb         = org.Hs.eg.db,        # human annotation
  keyType       = "SYMBOL",            # our genes are gene symbols
  ont           = "BP",                # Biological Process
  universe      = universe,            # background gene set
  pAdjustMethod = "BH",                # Benjamini-Hochberg multiple-testing correction
  pvalueCutoff  = 0.05, qvalueCutoff = 0.2
)
res <- as.data.frame(ego)              # results as a regular table
write_csv(res, file.path(project_dir, "Result", "go_bp_by_cluster.csv"))
print(table(res$Cluster))              # number of enriched terms per cluster

p <- dotplot(ego, showCategory = 3, font.size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(project_dir, "Figure", "go_bp_dotplot.png"), p, width = 11, height = 10, dpi = 300)
p
top_terms <- res %>%
  group_by(Cluster) %>% slice_min(p.adjust, n = 3, with_ties = FALSE) %>% ungroup() %>%
  dplyr::select(Cluster, Description, GeneRatio, p.adjust)     # dplyr:: avoids a name clash with another package
print(top_terms, n = 60)
writeLines(capture.output(sessionInfo()), file.path(project_dir, "Result", "R_sessionInfo.txt"))

# (a) How many enriched terms does each cluster have? A cluster missing from the dotplot has none.
table(res$Cluster)

# (b) Genes that drive one term in one cluster; here: oxidative phosphorylation in cluster 3
genes_of <- function(cluster, term) {
  res %>%
    filter(Cluster == cluster, Description == term) %>%
    pull(geneID) %>%                     # genes are stored as one text string separated by "/"
    strsplit("/") %>% unlist()
}
ox <- genes_of("3", "oxidative phosphorylation")
ox                                       # the genes behind the term
sum(startsWith(ox, "MT-"))               # how many are mitochondrially encoded

# (c) The same for ribosome-related terms in cluster 11 (change the term if the name differs in your table)
head(genes_of("11", "cytoplasmic translation"), 30)

# Sensitivity analysis: remove cytoplasmic ribosomal protein genes (RPS/RPL) from the marker lists and repeat the test
sig_nr <- sig %>% filter(!grepl("^RP[SL][0-9]", names))    # drop symbols like RPS3A, RPL5, RPL23A (mitochondrial MRPS/MRPL are kept)
cat("genes removed:", nrow(sig) - nrow(sig_nr), "\n")
gene_list_nr <- split(sig_nr$names, sig_nr$group)           # one vector of gene symbols per cluster

ego_nr <- compareCluster(
  geneClusters  = gene_list_nr,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  universe      = universe,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
res_nr <- as.data.frame(ego_nr)
write_csv(res_nr, file.path(project_dir, "Result", "go_bp_by_cluster_no_ribosomal.csv"))
table(res_nr$Cluster)                                        # compare with the original counts, especially clusters 0, 2 and 11