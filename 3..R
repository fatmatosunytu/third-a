#========================== Install and load required packages for analysis ===================================

# ---- INPUT PATHS ----
root_dir  <- getwd()
meta_dir  <- normalizePath("C:/Users/User/Documents/metadata", winslash = "/")

cel_files <- list.files(meta_dir, pattern = "\\.CEL$", full.names = TRUE)
pheno_csv <- file.path(meta_dir, "pheno.csv")

# (opsiyonel) pheno dosyası varsa oku
if (file.exists(pheno_csv)) {
  pheno <- read.csv(pheno_csv, stringsAsFactors = FALSE)
}

# Affymetrix okumaları
library(affy)  # veya oligo platformuna göre
raw <- ReadAffy(filenames = cel_files)
# ... RMA, QC, batch correction, limma vs. burada devam ...

packages <- c(
  "pkgbuild", "AnnotationDbi", "Biobase", "DOSE", "GEOquery", "GOSemSim", "R.utils",
  "affy", "annotate", "clusterProfiler", "data.table", "dplyr", "enrichplot",
  "genefilter", "ggplot2", "ggrepel", "grid", "gridExtra", "hgu133plus2.db",
  "jsonlite", "limma", "msigdbr", "multiMiR", "org.Hs.eg.db", "pheatmap",
  "sva", "umap", "STRINGdb", "igraph", "ggraph", "tidygraph"
)

#------------------------   Load required packages; assumes they are already installed   ----------------------

invisible(lapply(packages, library, character.only = TRUE))

align_samples <- function(cel_files, pheno) {
  sample_names <- basename(cel_files)
  pheno_df <- as.data.frame(pheno, stringsAsFactors = FALSE)
  pheno_keys <- toupper(basename(trimws(pheno_df$filename)))
  sample_keys <- toupper(sample_names)
  idx <- match(sample_keys, pheno_keys)
  if (any(is.na(idx))) {
    stop("pheno$filename içinde bulunamayan örnek(ler): ",
         paste(sample_names[is.na(idx)], collapse = ", "))
  }
  aligned <- pheno_df[idx, , drop = FALSE]
  rownames(aligned) <- sample_names
  stopifnot(nrow(aligned) == length(sample_names))
  stopifnot(identical(rownames(aligned), sample_names))
  aligned
}

pkgbuild::has_build_tools(debug = TRUE)

args <- commandArgs(trailingOnly = TRUE)
FDR_THRESH   <- if (length(args) >= 1) as.numeric(args[1]) else 0.05
LOGFC_THRESH <- if (length(args) >= 2) as.numeric(args[2]) else 1

# ---- utils: NULL-coalescing ve güvenli multiMiR çekici ----
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_multimir_table <- function(targets, table = c("validated","predicted"),
                                chunk_size = 500, timeout_sec = 120) {
  table <- match.arg(table)
  if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
  n <- length(targets); if (n == 0) return(NULL)
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    tg <- unique(targets[idx[[i]]])
    message(sprintf("multiMiR '%s' chunk %d/%d (n=%d)", table, i, length(idx), length(tg)))
    flush.console()
    mm <- try(
      R.utils::withTimeout(
        multiMiR::get_multimir(org = "hsa", target = tg, table = table, summary = FALSE),
        timeout = timeout_sec, onTimeout = "silent"
      ),
      silent = TRUE
    )
    if (inherits(mm, "try-error") || is.null(mm)) next
    if (!methods::is(mm, "multiMiR") || nrow(mm@data) == 0) next
    out[[i]] <- mm@data
  }
  if (all(vapply(out, is.null, logical(1)))) return(NULL)
  do.call(rbind, out)
}

#-----------------------------   Load CEL files and normalize with RMA    ---------------------------------------------

root_dir <- normalizePath(Sys.getenv("POMPE_ROOT", getwd()), winslash = "/")
meta_dir <- normalizePath("C:/Users/User/Documents/metadata", winslash = "/")

stopifnot(dir.exists(meta_dir))

dir.create(file.path(root_dir, "results", "qc"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "results", "deg"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "results", "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "results", "miRNA_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "logs"),             recursive = TRUE, showWarnings = FALSE)

writeLines(c(capture.output(sessionInfo())), file.path(root_dir, "logs", "session_info.txt"))

set.seed(20240724)

pheno <- data.table::fread(file.path(meta_dir, "pheno.csv"),
                           na.strings = c("", "NA", "NaN"))

pheno$filename <- basename(trimws(pheno$filename))
cel_dir   <- meta_dir
cel_files <- file.path(cel_dir, pheno$filename)
missing   <- pheno$filename[!file.exists(cel_files)]
if (length(missing)) stop("Missing CEL files in metadata/: ", paste(missing, collapse = ", "))

stopifnot(all(c("sample","group","filename") %in% names(pheno)))
g <- tolower(trimws(as.character(pheno$group)))
g[g %in% c("control","healthy","normal")] <- "control"
g[grepl("pompe", g)] <- "pompe"
pheno$group <- factor(ifelse(g %in% c("control","pompe"),
                             ifelse(g == "control","Control","Pompe"), NA_character_),
                      levels = c("Control","Pompe"))
stopifnot(!any(is.na(pheno$group)))

if (!"is_melas" %in% names(pheno)) pheno$is_melas <- FALSE
pheno$is_melas <- tolower(trimws(as.character(pheno$is_melas))) %in%
  c("1","true","t","yes","y")

melas_idx <- which(pheno$is_melas)
if (length(melas_idx) > 1) {
  warning("More than one MELAS sample detected; using the first flagged sample.")
  melas_idx <- melas_idx[1]
}

if (!"batch" %in% names(pheno)) pheno$batch <- NA

# ----------------   ExpressionSet ID eşitleme (assayData / phenoData / protocolData)   -----------------------

pheno_aligned <- align_samples(cel_files, pheno)
if (inherits(pheno_aligned, "data.table")) data.table::setDF(pheno_aligned)

if (anyDuplicated(rownames(pheno_aligned))) {
  rownames(pheno_aligned) <- make.unique(rownames(pheno_aligned))
}

if (!exists("eset")) {
  stopifnot(exists("cel_files"))
  raw  <- ReadAffy(filenames = cel_files)
  eset <- rma(raw)
  expr <- exprs(eset)
}

Biobase::sampleNames(eset) <- rownames(pheno_aligned)

expr <- Biobase::exprs(eset)                 
colnames(expr) <- Biobase::sampleNames(eset) 

ord <- match(Biobase::sampleNames(eset), rownames(pheno_aligned))
stopifnot(!any(is.na(ord)))                
pheno_aligned <- pheno_aligned[ord, , drop = FALSE]

Biobase::pData(eset) <- pheno_aligned
Biobase::protocolData(eset) <- Biobase::AnnotatedDataFrame(
  data.frame(row.names = Biobase::sampleNames(eset))
)

methods::validObject(eset)  # TRUE dönmeli
expr_matrix <- Biobase::exprs(eset)
colnames(expr_matrix) <- Biobase::sampleNames(eset)

#========================== Conduct quality control and filter low-expression probes ===================================

iqr <- apply(expr, 1, IQR)
thr <- quantile(iqr, 0.20, na.rm = TRUE)   # alt %20 IQR dışarı
keep <- iqr > thr
expr_f <- expr[keep, ]
message("Low-expression removed: ", sum(!keep), " probes (", round(mean(!keep)*100,2), "%)")
write.csv(data.frame(probe=rownames(expr)[!keep], IQR=iqr[!keep]),
          file.path(root_dir,"results","qc","low_expression_removed.csv"), row.names = FALSE)

pca <- prcomp(t(expr_f), scale. = TRUE)
pca_df <- data.frame(pca$x[,1:2], group = pheno_aligned$group,
                     batch = pheno_aligned$batch, sample = pheno_aligned$sample)
gg <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = factor(batch), label = sample)) +
  geom_point(size=3) + ggrepel::geom_text_repel(size=2.8) +
  theme_bw(14) + labs(title="PCA – RMA + IQR filtre", shape = "batch")
ggplot2::ggsave(file.path(root_dir,"results","qc","pca_groups_batches.png"), gg, width=9, height=7)

#========================== Identify differentially expressed genes with limma ===================================

expr_cb <- if (exists("expr_f")) expr_f else expr
stopifnot(is.matrix(expr_cb))
stopifnot(exists("pheno_aligned"))
groups <- factor(pheno_aligned$group, levels = c("Control","Pompe"))
design <- model.matrix(~0 + groups)
colnames(design) <- c("Control","Pompe")
stopifnot(ncol(expr_cb) == nrow(design))

print(dim(expr_cb)); print(dim(design)); head(design)

  # 5) limma

  fit  <- limma::lmFit(expr_bc, design_cov)
  lc_p <- grep("^groupPompe$", colnames(design_cov), value = TRUE)
  lc_c <- grep("^groupControl$", colnames(design_cov), value = TRUE)
  stopifnot(length(lc_p)==1, length(lc_c)==1)
  cont <- limma::makeContrasts(PompeVsControl = !!as.name(lc_p) - !!as.name(lc_c),
                               levels = design_cov)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  tt_all <- limma::topTable(fit2, coef="PompeVsControl", number=Inf, adjust.method="BH")
  
  # 6) anotasyon + collapse + DEG filtre
  aa <- .annotate_and_collapse(tt_all)
  deg <- subset(aa$tt_gene, adj.P.Val < FDR_THRESH & abs(logFC) >= LOGFC_THRESH)
  
  # 7) yazımlar
  out_dir <- file.path(root_dir, "results", "deg", out_tag)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(tt_all,             file.path(out_dir, "all_probes_BH.csv"))
  write.csv(aa$tt_gene,         file.path(out_dir, "collapsed_by_gene.csv"), row.names=TRUE)
  write.csv(deg,                file.path(out_dir, sprintf("DEG_FDR_lt_%s_logFC_ge_%s.csv", FDR_THRESH, LOGFC_THRESH)),
            row.names = TRUE)
  
  # küçük özet & dönüş
  up_n   <- sum(deg$logFC > 0, na.rm = TRUE)
  down_n <- sum(deg$logFC < 0, na.rm = TRUE)
  message(sprintf("[%s] samples=%d, DEGs=%d (up=%d, down=%d)",
                  out_tag, ncol(expr_use), nrow(deg), up_n, down_n))
  
  return(list(deg=deg, tt_all=tt_all, tt_gene=aa$tt_gene,
              n_up=up_n, n_down=down_n,
              out_dir=out_dir, n_samples=ncol(expr_use), tag=out_tag))

# ===================== PPI analysis via STRINGdb =====================

run_string_ppi_smart <- function(deg_tbl, root_dir,
                                 score_min   = 700,
                                 top_nodes   = 200,
                                 label_top_n = 20) {
  stopifnot(nrow(deg_tbl) > 1)
  
  ppi_dir <- file.path(root_dir, "results", "ppi")
  dir.create(ppi_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1) STRING'e map
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
  map_col   <- if ("ENTREZID" %in% colnames(deg_tbl)) "ENTREZID" else "SYMBOL"
  deg_map   <- string_db$map(deg_tbl, map_col, removeUnmappedRows = TRUE)
  if (nrow(deg_map) < 2) { message("STRING map yetersiz."); return(invisible(NULL)) }
  
  # 2) Etkileşimleri çek
  hits      <- deg_map$STRING_id
  ppi_edges <- string_db$get_interactions(hits)
  if (is.null(ppi_edges) || !nrow(ppi_edges)) { message("PPI kenarı yok."); return(invisible(NULL)) }
  
  # 3) Temizlik ve skor filtresi
  if ("combined_score" %in% names(ppi_edges)) {
    ppi_edges <- subset(ppi_edges, combined_score >= score_min)
  }
  # self-loop ve duplikatlar
  ppi_edges <- subset(ppi_edges, protein1 != protein2)
  ppi_edges$u <- pmin(ppi_edges$protein1, ppi_edges$protein2)
  ppi_edges$v <- pmax(ppi_edges$protein1, ppi_edges$protein2)
  ppi_edges <- unique(ppi_edges[, c("u","v", setdiff(names(ppi_edges), c("protein1","protein2","u","v")))])
  names(ppi_edges)[1:2] <- c("protein1","protein2")
  
  # 4) ID → SYMBOL eşleme (etiketler için)
  id2sym <- unique(deg_map[, c("STRING_id","SYMBOL","logFC")])
  ppi_edges$SYMBOL1 <- id2sym$SYMBOL[match(ppi_edges$protein1, id2sym$STRING_id)]
  ppi_edges$SYMBOL2 <- id2sym$SYMBOL[match(ppi_edges$protein2, id2sym$STRING_id)]
  
  # 5) Grafik için graph objesi
  g <- igraph::graph_from_data_frame(ppi_edges[, c("protein1","protein2")], directed = FALSE)
  if (igraph::ecount(g) == 0) { message("Graf boş."); return(invisible(NULL)) }
  
  # En büyük bağlı komponent
  comps <- igraph::components(g)
  g <- igraph::induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  
  # Düğüm meta: SYMBOL ve logFC
  V(g)$SYMBOL <- id2sym$SYMBOL[match(names(V(g)), id2sym$STRING_id)]
  V(g)$logFC  <- id2sym$logFC [match(names(V(g)), id2sym$STRING_id)]
  V(g)$dir    <- ifelse(is.na(V(g)$logFC), "NA",
                        ifelse(V(g)$logFC > 0, "Up", "Down"))
  
  # En yoğun top_nodes'a indirgeme
  degv <- igraph::degree(g)
  keep <- names(sort(degv, decreasing = TRUE))[seq_len(min(top_nodes, length(degv)))]
  g2   <- igraph::induced_subgraph(g, keep)
  
  # Hub etiketleri
  degv2 <- igraph::degree(g2)
  hubs  <- names(sort(degv2, decreasing = TRUE))[seq_len(min(label_top_n, length(degv2)))]
  V(g2)$label <- ifelse(names(V(g2)) %in% hubs, V(g2)$SYMBOL, NA)
  
  # Kenarlara skor (kalınlık) ekle
  e_df <- igraph::as_data_frame(g2, what = "edges")
  e_key <- paste(ppi_edges$protein1, ppi_edges$protein2)
  w     <- ppi_edges$combined_score[match(paste(e_df$from, e_df$to), e_key)]
  igraph::E(g2)$w <- scales::rescale(w %||% 700, to = c(0.2, 1.5))
  
  # 6) Görselleştirme
  set.seed(42)
  lay <- igraph::layout_with_fr(g2)
  plt <- ggraph::ggraph(g2, layout = lay) +
    ggraph::geom_edge_link(aes(width = ..index..), alpha = 0.3, show.legend = FALSE) +
    ggraph::geom_node_point(aes(shape = dir), size = 3) +
    ggraph::geom_node_text(aes(label = label), repel = TRUE, size = 3) +
    scale_shape_manual(values = c(Up = 16, Down = 1, NA = 4)) +
    ggraph::theme_void()
  
  # 7) Çıktılar
  dir.create(ppi_dir, showWarnings = FALSE, recursive = TRUE)
  # tam tablo
  write.csv(ppi_edges, file.path(ppi_dir, "ppi_edges_full.csv"), row.names = FALSE)
  # filtreli düğüm/kenar tabloları
  nodes_out <- data.frame(STRING_id = names(V(g2)),
                          SYMBOL    = V(g2)$SYMBOL,
                          logFC     = V(g2)$logFC,
                          degree    = degv2,
                          dir       = V(g2)$dir,
                          stringsAsFactors = FALSE)
  edges_out <- igraph::as_data_frame(g2, what = "edges")
  write.csv(nodes_out, file.path(ppi_dir, "ppi_nodes_filtered.csv"), row.names = FALSE)
  write.csv(edges_out, file.path(ppi_dir, "ppi_edges_filtered.csv"), row.names = FALSE)
  ggsave(file.path(ppi_dir, "ppi_network_filtered.png"), plt, width = 9, height = 7, dpi = 300)
  
  message("STRING mapped: ", nrow(deg_map), "/", nrow(deg_tbl),
          " | edges(full): ", nrow(ppi_edges),
          " | nodes(filtered): ", igraph::vcount(g2),
          " | edges(filtered): ", igraph::ecount(g2))
  invisible(list(graph = g2, nodes = nodes_out, edges = edges_out))
}

# ÇAĞRI:
deg_for_ppi <- if (exists("melas_res") && !is.null(melas_res$deg_noM_sig)) melas_res$deg_noM_sig else deg
run_string_ppi_smart(deg_for_ppi, root_dir, score_min = 700, top_nodes = 200, label_top_n = 20)

nrow(tt_all)
table("SYMBOL_is_NA" = is.na(tt_all$SYMBOL))
mean(!is.na(tt_all$SYMBOL))  # kapsama oranı
probes <- rownames(expr_cb)
all(rownames(tt_all) %in% probes)

head(AnnotationDbi::select(hgu133plus2.db,
                           keys=head(probes, 5),
                           columns=c("SYMBOL","GENENAME"),
                           keytype="PROBEID"))


####### OPSİYONEL ###########

# tt_all senin tüm gen tablon
#all_ids <- tt_all$SYMBOL[!is.na(tt_all$SYMBOL)]
#mapped   <- string_db$get_aliases(all_ids)  # mevcut eşleşmeler
#mapped_ids <- unique(mapped$preferred_name)

#not_mapped <- setdiff(all_ids, mapped_ids)
#length(not_mapped)        # kaç tanesi eşleşmedi?
#head(not_mapped, 20)      # ilk 20’sini göster

# toplam diferansiyel gen
nrow(deg)

# yukarı ve aşağı düzenlenen ayrı ayrı
up_n   <- sum(deg$logFC > 0, na.rm=TRUE)
down_n <- sum(deg$logFC < 0, na.rm=TRUE)

cat("\n=== ÖZET ===\n")
cat("Design boyutu: ", nrow(design), " örnek x ", ncol(design), " değişken\n", sep = "")
cat("Toplam gen (probe-level): ", nrow(tt_all), "\n", sep = "")
cat("Tekilleştirilmiş gen sayısı: ", nrow(tt_gene), "\n", sep = "")
write.csv(deg,     file.path(root_dir,"results","deg", sprintf("DEG_FDR_lt_%s.csv", FDR_THRESH)), row.names = TRUE)
cat("İlk 5 DEG:\n"); print(utils::head(deg[, c("SYMBOL","logFC","adj.P.Val","GENENAME")], 5))
cat("Toplam DEG:", nrow(deg), "\n")
cat("Up-regulated:", up_n, "\n")
cat("Down-regulated:", down_n, "\n")


######## OPSİYONEL #######
# Kaç DEG STRING’de eşleşti?
mapped_n <- nrow(deg_mapped); total_n <- nrow(deg_for_ppi)
message("STRING mapped: ", mapped_n, "/", total_n)

# Kaç kenar geldi, min skor?
if (exists("ppi_edges")) {
  message("PPI edges: ", nrow(ppi_edges))
  if ("combined_score" %in% names(ppi_edges))
    message("score range: ", range(ppi_edges$combined_score, na.rm=TRUE))
}

# PNG/CSV üretildi mi?
print(file.exists(file.path(ppi_dir,"ppi_edges.csv")))
print(file.exists(file.path(ppi_dir,"ppi_network.png")))

de_sig <- deg             # FDR < FDR_THRESH gene-level
expr_use <- expr_cb       # batch/SVA sonrası matris (grafikler için)

if (!exists("deg")) {
  p_thresh   <- get0("p_thresh",   ifnotfound = get0("FDR_THRESH",   ifnotfound = 0.05))
  lfc_thresh <- get0("lfc_thresh", ifnotfound = get0("LOGFC_THRESH", ifnotfound = 1))
  stopifnot(exists("tt_gene"))
  deg <- subset(tt_gene, adj.P.Val < p_thresh & abs(logFC) >= lfc_thresh)
}
total_deg <- nrow(deg)
up_deg    <- sum(deg$logFC > 0, na.rm = TRUE)
down_deg  <- sum(deg$logFC < 0, na.rm = TRUE)
cat("Toplam DEG:", total_deg, "| Up:", up_deg, "| Down:", down_deg, "\n")

to_scalar_logical <- function(x) {
  if (is.logical(x)) return(isTRUE(any(x, na.rm = TRUE)))
  if (is.numeric(x)) return(isTRUE(any(x != 0, na.rm = TRUE)))
  if (is.character(x)) {
    y <- tolower(trimws(x))
    return(isTRUE(any(y %in% c("1","true","t","yes","y"), na.rm = TRUE)))
  }
  FALSE
}

prepare_melas <- function(pheno, expr, results_dir = file.path(getwd(), "results")) {
  ids <- colnames(expr)
  if (!"filename" %in% names(pheno))
    stop("pheno_aligned içinde 'filename' kolonu yok (GSM...CEL).")
  
  idx <- match(ids, pheno$filename)
  stopifnot(!any(is.na(idx)))
  pheno <- pheno[idx, , drop = FALSE]
  rownames(pheno) <- ids
  
  run_res <- run_deg_analysis(expr, pheno, run_melas = TRUE)
  deg_cov    <- run_res$deg_cov
  deg_cov_sig <- run_res$deg_cov_sig
  deg_noM    <- run_res$deg_noM
  deg_noM_sig <- run_res$deg_noM_sig
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(deg_cov,  file.path(results_dir, "DE_all_covariate_MELAS_included.csv"))
  write.csv(deg_noM,  file.path(results_dir, "DE_all_noMELAS.csv"))
  write.csv(deg_cov_sig, file.path(results_dir,
                                   sprintf("DE_sig_covariate_MELAS_included_FDR%s_LFC%s.csv",
                                           FDR_THRESH, LOGFC_THRESH)))
  write.csv(deg_noM_sig, file.path(results_dir,
                                   sprintf("DE_sig_noMELAS_FDR%s_LFC%s.csv",
                                           FDR_THRESH, LOGFC_THRESH)))
  
  list(pheno_aligned = pheno,
       deg_cov = deg_cov, deg_cov_sig = deg_cov_sig,
       deg_noM = deg_noM, deg_noM_sig = deg_noM_sig)
}

melas_flag <-
  ( "is_melas" %in% names(pheno_aligned) && to_scalar_logical(pheno_aligned$is_melas) ) ||
  ( "specialcase" %in% names(pheno_aligned) && to_scalar_logical(pheno_aligned$specialcase) )

pheno_aligned$is_melas <- factor(ifelse(melas_flag, "Yes","No"), levels = c("No","Yes"))

pheno_aligned$group <- factor(pheno_aligned$group, levels = c("Control","Pompe"))
if ("age_at_baseline_mounth" %in% names(pheno_aligned))
  pheno_aligned$age_at_baseline_mounth <- suppressWarnings(as.numeric(pheno_aligned$age_at_baseline_mounth))
if ("sex" %in% names(pheno_aligned))
  pheno_aligned$sex <- factor(pheno_aligned$sex)  # "M"/"F"
if ("rin" %in% names(pheno_aligned))
  pheno_aligned$rin <- suppressWarnings(as.numeric(pheno_aligned$rin))

stopifnot(nrow(pheno_aligned) == ncol(expr_f))
stopifnot(identical(rownames(pheno_aligned), colnames(expr_f)))

run_deg_analysis <- function(expr, pheno, run_melas = FALSE,
                             lfc_threshold = LOGFC_THRESH, fdr_threshold = FDR_THRESH) {
  design_cov <- model.matrix(~0 + group + is_melas + age_at_baseline_mounth + sex + rin,
                             data = pheno)
  colnames(design_cov) <- make.names(colnames(design_cov))
  fit_cov  <- limma::lmFit(expr, design_cov)
  cont_cov <- limma::makeContrasts(Pompe - Control, levels = design_cov)
  fit2_cov <- limma::eBayes(limma::contrasts.fit(fit_cov, cont_cov))
  cov_res  <- limma::topTable(fit2_cov, coef = 1, number = Inf, adjust.method = "BH")
  out <- list(deg_cov = cov_res,
              deg_cov_sig = subset(cov_res, adj.P.Val < fdr_threshold & abs(logFC) >= lfc_threshold))
  if (run_melas) {
    keep_idx <- which(pheno$is_melas == "No")
    expr_noM <- expr[, keep_idx, drop = FALSE]
    ph_noM   <- droplevels(pheno[keep_idx, ])
    ph_noM$group <- factor(ph_noM$group, levels = c("Control","Pompe"))
    design_noM <- model.matrix(~0 + group + age_at_baseline_mounth + sex + rin, data = ph_noM)
    colnames(design_noM) <- make.names(colnames(design_noM))
    fit_noM  <- limma::lmFit(expr_noM, design_noM)
    cont_noM <- limma::makeContrasts(Pompe - Control, levels = design_noM)
    fit2_noM <- limma::eBayes(limma::contrasts.fit(fit_noM, cont_noM))
    deg_noM    <- limma::topTable(fit2_noM, coef = 1, number = Inf, adjust.method = "BH")
    out$deg_noM <- deg_noM
    out$deg_noM_sig <- subset(deg_noM, adj.P.Val < fdr_threshold & abs(logFC) >= lfc_threshold)
  }
  out
}

plot_dir <- file.path(root_dir, "results", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f  # batch/SVA sonrası
stopifnot(ncol(expr_plot) == nrow(pheno_aligned))

stopifnot(exists("tt_all"), exists("tt_gene"), exists("deg"))

volc_main <- within(tt_all, {
  negLog10FDR <- -log10(adj.P.Val)
  sig <- adj.P.Val < FDR_THRESH
  dir <- ifelse(sig & logFC > 0, "Up",
                ifelse(sig & logFC < 0, "Down", "NS"))
})

lab_main <- head(volc_main[order(volc_main$adj.P.Val, -abs(volc_main$logFC)), ], 15)
lab_main$SYMBOL[is.na(lab_main$SYMBOL) | lab_main$SYMBOL == ""] <- rownames(lab_main)[is.na(lab_main$SYMBOL) | lab_main$SYMBOL == ""]

volc <- ggplot(volc_main, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(shape = dir), alpha = 0.6) +
  geom_hline(yintercept = -log10(FDR_THRESH), linetype = 2) +
  geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = 3) +
  ggrepel::geom_text_repel(
    data = lab_main,
    aes(label = SYMBOL),
    size = 3, max.overlaps = 100
  ) +
  theme_bw(14) +
  labs(title = "Volcano (Pompe vs Control) — Main",
       x = "log2 Fold Change", y = "-log10(FDR)")
ggsave(file.path(plot_dir, "volcano_main.png"), volc, width = 8, height = 6, dpi = 300)

if (exists("tt_all_s") && exists("deg_sens")) {
  volc_s <- within(tt_all_s, {
    negLog10FDR <- -log10(adj.P.Val)
    sig <- adj.P.Val < FDR_THRESH
    dir <- ifelse(sig & logFC > 0, "Up",
                  ifelse(sig & logFC < 0, "Down", "NS"))
  })
  lab_s <- head(volc_s[order(volc_s$adj.P.Val, -abs(volc_s$logFC)), ], 15)
  lab_s$SYMBOL[is.na(lab_s$SYMBOL) | lab_s$SYMBOL == ""] <- rownames(lab_s)[is.na(lab_s$SYMBOL) | lab_s$SYMBOL == ""]
  
  p_s <- ggplot(volc_s, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(shape = dir), alpha = 0.6) +
    geom_hline(yintercept = -log10(FDR_THRESH), linetype = 2) +
    geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = 3) +
    ggrepel::geom_text_repel(
      data = lab_s,
      aes(label = SYMBOL),
      size = 3, max.overlaps = 100
    ) +
    theme_bw(14) +
    labs(title = "Volcano (Pompe vs Control) — MELAS excluded",
         x = "log2 Fold Change", y = "-log10(FDR)")
  ggsave(file.path(plot_dir, "volcano_noMELAS_main.png"), p_s, width = 8, height = 6, dpi = 300)
}

stopifnot(exists("tt_gene"), exists("deg"))
stopifnot(nrow(deg) > 0)

TOP_N <- min(50, nrow(deg))
stopifnot(TOP_N >= 1)

tt_gene$SYMBOL <- as.character(tt_gene$SYMBOL)
deg$SYMBOL     <- as.character(deg$SYMBOL)

top_syms <- head(deg$SYMBOL[order(deg$adj.P.Val)], TOP_N)

reps <- rownames(tt_gene)[ match(top_syms, tt_gene$SYMBOL) ]
reps <- unique(reps[!is.na(reps)])

sub_expr <- expr_plot[ intersect(rownames(expr_plot), reps), , drop = FALSE ]
stopifnot(nrow(sub_expr) > 0)   # mantıksal test

rownames(sub_expr) <- tt_gene$SYMBOL[ match(rownames(sub_expr), rownames(tt_gene)) ]

ann_col <- data.frame(Group = pheno_aligned$group,
                      row.names = rownames(pheno_aligned))
if ("batch" %in% names(pheno_aligned))
  ann_col$Batch <- pheno_aligned$batch

if ("is_melas" %in% names(pheno_aligned)) {
  melas_raw <- pheno_aligned$is_melas
  melas_log <- if (is.logical(melas_raw)) {
    melas_raw
  } else {
    v <- tolower(trimws(as.character(melas_raw)))
    v %in% c("1","true","t","yes","y")
  }
  melas_log[is.na(melas_log)] <- FALSE  # boşları "No" say
  ann_col$MELAS <- ifelse(melas_log, "Yes", "No")
}

ord <- order(ann_col$Group, ann_col$Batch, na.last = TRUE)
ann_col <- ann_col[ord, , drop = FALSE]
sub_expr <- sub_expr[, rownames(ann_col), drop = FALSE]

if (TOP_N >= 1) {  
  sub_expr_z <- t(scale(t(sub_expr)))  # satır bazlı z-score
  
  if (nrow(sub_expr_z) >= 2) {
    ph <- pheatmap::pheatmap(sub_expr_z,
                             annotation_col = ann_col,
                             show_rownames = TRUE, show_colnames = TRUE,
                             clustering_distance_rows = "correlation",
                             clustering_distance_cols = "correlation",
                             clustering_method = "average",
                             main = sprintf("Top %d DEG (z-score by gene) – Pompe vs Control", nrow(sub_expr_z)),
                             filename = file.path(plot_dir, sprintf("heatmap_top_FDR%s_main.png", FDR_THRESH))
    )
    
  } else if (nrow(sub_expr_z) == 1) {
    ph <- pheatmap::pheatmap(sub_expr_z,
                             annotation_col = ann_col,
                             show_rownames = TRUE, show_colnames = TRUE,
                             cluster_rows = FALSE, cluster_cols = TRUE,
                             main = "Only 1 gene passed filter (no row clustering)",
                             filename = file.path(plot_dir, sprintf("heatmap_top_FDR%s_main.png", FDR_THRESH))
    )
    
  } else {
    message("Heatmap atlandı: seçilen genlerden hiçbiri ifade matrisinde bulunamadı (probe vs gene adı?).")
  }
  
} else {
  message(sprintf("Heatmap atlandı: FDR<%s ile yeterli gen yok.", FDR_THRESH))
}


if (exists("tt_gene_s") && exists("expr_sens") && exists("ph_sens")) {
  top_gene_tbl_s <- head(tt_gene_s[order(tt_gene_s$adj.P.Val, -abs(tt_gene_s$logFC)), ], TOP_N)
  top_probes_s <- rownames(top_gene_tbl_s)
  sub_expr_s <- expr_sens[intersect(rownames(expr_sens), top_probes_s), , drop = FALSE]
  sub_expr_s_z <- t(scale(t(sub_expr_s)))
  
  ann_col_s <- data.frame(
    Group = ph_sens$group
  )
  if ("batch" %in% names(ph_sens)) ann_col_s$Batch <- ph_sens$batch
  rownames(ann_col_s) <- rownames(ph_sens)
  
  ord_s <- order(ann_col_s$Group, ann_col_s$Batch, decreasing = FALSE, na.last = TRUE)
  sub_expr_s_z <- sub_expr_s_z[, ord_s, drop = FALSE]
  ann_col_s    <- ann_col_s[ord_s, , drop = FALSE]
  
  png(file.path(plot_dir, sprintf("heatmap_top%d_genes_noMELAS.png", TOP_N)),
      width = 1100, height = 1400, res = 150)
  pheatmap(sub_expr_s_z,
           annotation_col = ann_col_s,
           show_rownames = TRUE, show_colnames = TRUE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           main = sprintf("Top %d genes (z-score) — Pompe vs Control (no MELAS)", TOP_N))
  dev.off()
}

message("Plots written to: ", plot_dir)

gsea_dir <- file.path(root_dir, "results", "deg")
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)

stopifnot(exists("tt_gene"))
ranks <- tt_gene$t
names(ranks) <- tt_gene$SYMBOL
ranks <- ranks[!is.na(ranks) & nzchar(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

if (length(ranks) < 50) {
  warning("GSEA için yeterli sayıda gen yok gibi (rank < 50). Sonuçlar güvenilmez olabilir.")
}

msig <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
  error = function(e) {
    msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  }
)
term2gene <- unique(msig[, c("gs_name","gene_symbol")])
gsea_go <- GSEA(geneList = ranks, TERM2GENE = term2gene,
                pvalueCutoff = FDR_THRESH, verbose = FALSE)

write.csv(as.data.frame(gsea_go),
          file.path(gsea_dir, "GSEA_GO_BP_msigdb.csv"), row.names = FALSE)

png(file.path(gsea_dir, "GSEA_GO_BP_dotplot.png"), width=1000, height=800, res=130)
print(dotplot(gsea_go, showCategory = 20))
dev.off()

if (nrow(as.data.frame(gsea_go)) > 0) {
  top_terms <- head(gsea_go@result$ID, 2)
  for (id in top_terms) {
    png(file.path(gsea_dir, paste0("GSEA_GO_BP_runningplot_", gsub("[^A-Za-z0-9]+","_", id), ".png")),
        width=1100, height=700, res=130)
    print(gseaplot2(gsea_go, geneSetID = id))
    dev.off()
  }
}

sym2ent <- bitr(tt_gene$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
sym2ent <- sym2ent[!is.na(sym2ent$ENTREZID) & nzchar(sym2ent$SYMBOL), ]
sym2ent$abs_rank <- abs(ranks[sym2ent$SYMBOL])
sym2ent <- sym2ent[order(sym2ent$abs_rank, decreasing = TRUE), ]
sym2ent_1to1 <- sym2ent[!duplicated(sym2ent$SYMBOL), c("SYMBOL","ENTREZID")]

ranks_kegg <- ranks[sym2ent_1to1$SYMBOL]
names(ranks_kegg) <- sym2ent_1to1$ENTREZID
ranks_kegg <- ranks_kegg[!duplicated(names(ranks_kegg))]
ranks_kegg <- sort(ranks_kegg, decreasing = TRUE)

if (length(ranks_kegg) >= 50) {
  gsea_kegg <- gseKEGG(geneList = ranks_kegg,
                       organism = "hsa",
                       pvalueCutoff = 0.1,
                       verbose = FALSE)
  write.csv(as.data.frame(gsea_kegg),
            file.path(gsea_dir, "GSEA_KEGG.csv"), row.names = FALSE)
  
  png(file.path(gsea_dir, "GSEA_KEGG_dotplot.png"), width=1000, height=800, res=130)
  print(dotplot(gsea_kegg, showCategory = 20))
  dev.off()
  
  top_kegg <- gsea_kegg@result$ID[1]
  if (!is.null(top_kegg) && !is.na(top_kegg)) {
    pathview::pathview(gene.data = ranks_kegg, pathway.id = top_kegg,
                       species = "hsa", out.suffix = "gsea",
                       kegg.native = TRUE, out.dir = gsea_dir)
  }
} else {
  message("KEGG GSEA atlandı: ENTREZ eşleşmesi az ( < 50 ).")
}

# Reactome pathway enrichment
entrez_deg <- sym2ent_1to1$ENTREZID[sym2ent_1to1$SYMBOL %in% deg$SYMBOL]
reactome <- ReactomePA::enrichPathway(gene = entrez_deg,
                                      pvalueCutoff = FDR_THRESH,
                                      readable = TRUE)
write.csv(as.data.frame(reactome),
          file.path(gsea_dir, "Reactome_enrichment.csv"), row.names = FALSE)
if (nrow(as.data.frame(reactome)) > 0) {
  png(file.path(gsea_dir, "Reactome_dotplot.png"), width=1000, height=800, res=130)
  print(dotplot(reactome, showCategory = 20))
  dev.off()
}

# GSVA/ssGSEA sample-level scores
gsva_dir <- file.path(root_dir, "results", "gsva")
dir.create(gsva_dir, showWarnings = FALSE, recursive = TRUE)
gs_list <- split(msig$gene_symbol, msig$gs_name)
probe2sym <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr_cb),
                                   columns = "SYMBOL", keytype = "PROBEID")
probe2sym <- probe2sym[!is.na(probe2sym$SYMBOL), ]
expr_gene <- rowsum(expr_cb[probe2sym$PROBEID, , drop = FALSE], probe2sym$SYMBOL)
gsva_scores <- GSVA::gsva(expr_gene, gs_list, method = "ssgsea", verbose = FALSE)
write.csv(gsva_scores, file.path(gsva_dir, "ssgsea_scores.csv"))

if (exists("tt_gene_s")) {
  ranks_s <- tt_gene_s$t
  names(ranks_s) <- tt_gene_s$SYMBOL
  ranks_s <- ranks_s[!is.na(ranks_s) & nzchar(names(ranks_s))]
  ranks_s <- sort(ranks_s, decreasing = TRUE)
  
  if (length(ranks_s) >= 50) {
    gsea_go_s <- GSEA(geneList = ranks_s, TERM2GENE = term2gene,
                      pvalueCutoff = FDR_THRESH, verbose = FALSE)
    write.csv(as.data.frame(gsea_go_s),
              file.path(gsea_dir, "GSEA_GO_BP_msigdb_noMELAS.csv"), row.names = FALSE)
    
    png(file.path(gsea_dir, "GSEA_GO_BP_dotplot_noMELAS.png"), width=1000, height=800, res=130)
    print(dotplot(gsea_go_s, showCategory = 20))
    dev.off()
  } else {
    message("no-MELAS GSEA atlandı: yeterli rank yok.")
  }
}

if (exists("ph")) {
  pdf(file.path(root_dir,"results","plots","volcano_heatmap_Fig1.pdf"), width=16, height=8, family="Times")
  grid.arrange(
    arrangeGrob(volc, top = textGrob("A", gp=gpar(fontsize=22, fontface="bold"), x=unit(0,"npc"), hjust=0)),
    arrangeGrob(ph$gtable, top = textGrob("B", gp=gpar(fontsize=22, fontface="bold"), x=unit(0,"npc"), hjust=0)),
    ncol = 2
  )
  dev.off()
}

# =======================    ENRICHMENT: GO / KEGG / GSEA + MELAS       =======================

enrich_dir <- file.path(root_dir, "results", "enrichment"); dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir   <- file.path(root_dir, "results", "plots");       dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

stopifnot(exists("tt_gene"), exists("deg"))
if (nrow(deg) == 0) stop("No DEGs detected; adjust thresholds.")
sig_symbols      <- unique(na.omit(as.character(deg$SYMBOL)))
universe_symbols <- unique(na.omit(as.character(tt_gene$SYMBOL)))

writeLines(sig_symbols,      file.path(enrich_dir, "_sig_gene_symbols.txt"))
writeLines(universe_symbols, file.path(enrich_dir, "_universe_gene_symbols.txt"))

sym2entrez      <- bitr(universe_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_universe <- unique(na.omit(sym2entrez$ENTREZID))
sig2entrez      <- bitr(sig_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_sig      <- unique(na.omit(sig2entrez$ENTREZID))

meta <- list(time=as.character(Sys.time()),
             n_sig_genes=length(sig_symbols), n_universe_genes=length(universe_symbols),
             n_sig_entrez=length(entrez_sig), n_universe_entrez=length(entrez_universe),
             adjust_method="BH", pvalueCutoff=FDR_THRESH, qvalueCutoff=FDR_THRESH)
writeLines(jsonlite::toJSON(meta, pretty=TRUE), file.path(enrich_dir, "_params_enrichment_main.json"))

go_bp <- enrichGO(gene = sig_symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", universe = universe_symbols,
                  pAdjustMethod = "BH", pvalueCutoff = FDR_THRESH, qvalueCutoff = FDR_THRESH,
                  readable = TRUE)
write.csv(as.data.frame(go_bp), file.path(enrich_dir,"GO_BP_main.csv"), row.names = FALSE)

ek_args <- list(gene = entrez_sig, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = FDR_THRESH)
if ("universe" %in% names(formals(clusterProfiler::enrichKEGG))) ek_args$universe <- entrez_universe
kegg <- do.call(clusterProfiler::enrichKEGG, ek_args)
write.csv(as.data.frame(kegg), file.path(enrich_dir,"KEGG_main.csv"), row.names = FALSE)

gl_df <- tt_gene[, c("SYMBOL","t","logFC")]
gl_df$t[is.na(gl_df$t)] <- gl_df$logFC[is.na(gl_df$t)]
gl_df <- merge(gl_df, sym2entrez, by.x="SYMBOL", by.y="SYMBOL", all.x=TRUE)
gl_df <- gl_df[!is.na(gl_df$ENTREZID), c("ENTREZID","t")]
geneList <- gl_df$t; names(geneList) <- gl_df$ENTREZID; geneList <- sort(geneList, decreasing = TRUE)

gsea_go   <- NULL
gsea_kegg <- NULL
if (length(geneList) >= 20) {
  gsea_go   <- try(gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                         ont = "BP", pAdjustMethod = "BH", verbose = FALSE), silent = TRUE)
  gsea_kegg <- try(gseKEGG(geneList = geneList, organism = "hsa", pAdjustMethod = "BH", verbose = FALSE), silent = TRUE)
  if (!inherits(gsea_go, "try-error"))   write.csv(as.data.frame(gsea_go),   file.path(enrich_dir,"GSEA_GO_BP_main.csv"),   row.names = FALSE)
  if (!inherits(gsea_kegg, "try-error")) write.csv(as.data.frame(gsea_kegg), file.path(enrich_dir,"GSEA_KEGG_main.csv"),    row.names = FALSE)
}
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  p_go_dot <- dotplot(go_bp, showCategory = 15, title = "GO BP (BH, explicit universe)")
  ggsave(file.path(plot_dir,"GO_BP_dotplot_main.png"), p_go_dot, width=10, height=9, dpi=300)
  p_go_bar <- barplot(go_bp, showCategory = 15, title = "GO BP (BH, explicit universe)")
  ggsave(file.path(plot_dir,"GO_BP_barplot_main.png"), p_go_bar, width=10, height=9, dpi=300)
  
  hsGO <- try(godata('org.Hs.eg.db', ont="BP"), silent=TRUE)
  if (!inherits(hsGO,"try-error")) {
    s <- try(pairwise_termsim(go_bp, semData=hsGO), silent=TRUE)
    if (!inherits(s,"try-error")) {
      p_emap <- try(emapplot(s, showCategory=15), silent=TRUE)
      if (!inherits(p_emap,"try-error")) ggsave(file.path(plot_dir,"GO_BP_emapplot_main.png"), p_emap, width=10, height=10, dpi=300)
    }
  }
}
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  p_kegg <- barplot(kegg, showCategory = 15, title = "KEGG (BH)")
  ggsave(file.path(plot_dir,"KEGG_barplot_main.png"), p_kegg, width=10, height=9, dpi=300)
}
if (!is.null(gsea_go) && !inherits(gsea_go,"try-error") && nrow(as.data.frame(gsea_go))>0) {
  p_gsea1 <- gseaplot2(gsea_go, geneSetID = head(gsea_go@result$ID,1), title = "GSEA GO BP (top set)")
  ggsave(file.path(plot_dir,"GSEA_GO_BP_topset_main.png"), p_gsea1, width=10, height=6, dpi=300)
}
if (!is.null(gsea_kegg) && !inherits(gsea_kegg,"try-error") && nrow(as.data.frame(gsea_kegg))>0) {
  p_gsea2 <- gseaplot2(gsea_kegg, geneSetID = head(gsea_kegg@result$ID,1), title = "GSEA KEGG (top pathway)")
  ggsave(file.path(plot_dir,"GSEA_KEGG_topset_main.png"), p_gsea2, width=10, height=6, dpi=300)
}

enrich_melas <- function(deg_tbl, tag) {
  sig_syms <- unique(na.omit(as.character(deg_tbl$SYMBOL)))
  if (length(sig_syms) < 3) return(invisible(NULL))
  sig_ent  <- unique(na.omit(bitr(sig_syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID))
  
  go <- enrichGO(gene=sig_syms, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                 ont="BP", universe=universe_symbols,
                 pAdjustMethod="BH", pvalueCutoff=FDR_THRESH, qvalueCutoff=FDR_THRESH, readable=TRUE)
  write.csv(as.data.frame(go), file.path(enrich_dir, paste0("GO_BP_",tag,".csv")), row.names=FALSE)
  
  ek_args <- list(gene=sig_ent, organism="hsa", pAdjustMethod="BH", pvalueCutoff=FDR_THRESH)
  if ("universe" %in% names(formals(clusterProfiler::enrichKEGG))) ek_args$universe <- entrez_universe
  kk <- do.call(clusterProfiler::enrichKEGG, ek_args)
  write.csv(as.data.frame(kk), file.path(enrich_dir, paste0("KEGG_",tag,".csv")), row.names=FALSE)
  
  if (!is.null(go) && nrow(as.data.frame(go))>0) {
    ggsave(file.path(plot_dir, paste0("GO_BP_dotplot_",tag,".png")),
           dotplot(go, showCategory=15, title=paste0("GO BP (",tag,")")), width=10, height=9, dpi=300)
  }
  if (!is.null(kk) && nrow(as.data.frame(kk))>0) {
    ggsave(file.path(plot_dir, paste0("KEGG_barplot_",tag,".png")),
           barplot(kk, showCategory=15, title=paste0("KEGG (",tag,")")), width=10, height=9, dpi=300)
  }
  invisible(list(go=go, kegg=kk))
}

if (exists("deg_sens")) {
  invisible(enrich_melas(deg_sens, "noMELAS"))
  
  has_melas_col <- "is_melas" %in% names(pheno_aligned)
  melas_log <- if (has_melas_col) {
    x <- pheno_aligned$is_melas
    if (is.logical(x)) x else {
      v <- tolower(trimws(as.character(x)))
      v %in% c("1","true","t","yes","y")
    }
  } else rep(FALSE, nrow(pheno_aligned))
  
  
} else if (exists("pheno_aligned") && any(melas_log, na.rm = TRUE)) {
  
  expr_base <- if (exists("expr_use")) expr_use else expr_cb
  
  keep <- !melas_log
  expr_s <- expr_base[, keep, drop = FALSE]
  ph_s   <- droplevels(pheno_aligned[keep, , drop = FALSE])
  
  if (nlevels(ph_s$group) >= 2 && all(table(ph_s$group) > 0)) {
    design_s <- model.matrix(~ 0 + group, data = ph_s)
    colnames(design_s) <- gsub("^group", "", colnames(design_s))  # "Control","Pompe"
    
    fit_s  <- limma::lmFit(expr_s, design_s)
    cont_s <- limma::makeContrasts(Pompe - Control, levels = design_s)
    fit2_s <- limma::eBayes(limma::contrasts.fit(fit_s, cont_s))
    
    tt_all_s <- limma::topTable(fit2_s, coef = 1, number = Inf, adjust.method = "BH")
    
    
    expr_s   <- expr_cb[, keep_idx, drop=FALSE]
    ph_s     <- droplevels(pheno_aligned[keep_idx, ])
    design_s <- model.matrix(~ 0 + group, data = ph_s); colnames(design_s) <- levels(ph_s$group)
    fit_s2   <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr_s, design_s),
                                                   limma::makeContrasts(Pompe - Control, levels = design_s)))
    
    tt_all_s$abs_logFC <- abs(tt_all_s$logFC)
    tt_gene_s <- tt_all_s[!is.na(tt_all_s$SYMBOL) & nzchar(tt_all_s$SYMBOL), ]
    tt_gene_s <- tt_gene_s[order(tt_gene_s$SYMBOL, -tt_gene_s$abs_logFC), ]
    tt_gene_s <- tt_gene_s[!duplicated(tt_gene_s$SYMBOL), ]
    deg_sens  <- subset(tt_gene_s, adj.P.Val < FDR_THRESH)
    write.csv(deg_sens,
              file.path(enrich_dir,
                        sprintf("DEG_noMELAS_FDR_lt_%s.csv", FDR_THRESH)),
              row.names=FALSE)
    invisible(enrich_melas(deg_sens, "noMELAS"))
    
  } else {
    message("MELAS hariçte iki grup kalmadı; Pompe-Control karşılaştırması atlandı.")
    tt_all_s <- NULL
    deg_sens <- NULL
  }
  
  sum_lines <- c(
    sprintf("Main GO terms (q<%s): %d", FDR_THRESH, ifelse(is.null(go_bp), 0, nrow(as.data.frame(go_bp)))),
    sprintf("Main KEGG pathways (q<%s): %d", FDR_THRESH, ifelse(is.null(kegg), 0, nrow(as.data.frame(kegg)))),
    sprintf("Universe size: SYMBOL=%d | ENTREZ=%d", length(universe_symbols), length(entrez_universe)),
    "Notes: BH correction everywhere; GO uses explicit universe (all tested gene-level).",
    "Sensitivity: If MELAS present, enrichment repeated without MELAS; GSEA added to reduce threshold dependence."
  )
  writeLines(sum_lines, file.path(enrich_dir,"_ENRICHMENT_SUMMARY.txt"))
  
  # ====================       miRNA Target Enrichment (FAIR, offline)         ====================
  
  out_dir  <- file.path(root_dir, "results", "miRNA_analysis")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(20240724)
  
  stopifnot(exists("tt_gene"), "SYMBOL" %in% names(tt_gene))
  stopifnot(exists("deg"),     "SYMBOL" %in% names(deg))
  if (nrow(deg) == 0) stop("No DEGs detected; adjust thresholds.")
  universe_symbols <- unique(tt_gene$SYMBOL[!is.na(tt_gene$SYMBOL) & nzchar(tt_gene$SYMBOL)])
  sig_genes        <- unique(deg$SYMBOL[!is.na(deg$SYMBOL) & nzchar(deg$SYMBOL)])
  
  meta <- list(
    dataset          = "GSE38680",
    n_control        = if (exists("groups")) sum(groups=="Control") else NA_integer_,
    n_pompe          = if (exists("groups")) sum(groups=="Pompe") else NA_integer_,
    n_universe_genes = length(universe_symbols),
    n_sig_genes      = length(sig_genes),
    fdr_threshold    = FDR_THRESH,
    adjust_method    = "BH",
    note_small_n     = "Small sample size; interpret enrichment cautiously."
  )
  jsonlite::write_json(meta, file.path(out_dir, "_mirna_meta.json"), pretty = TRUE)
  
  
  enrich_mirna <- function(gene_set, universe_set, which_table=c("validated","predicted")) {
    which_table <- match.arg(which_table)
    df <- safe_multimir_table(universe_set, which_table)
    if (is.null(df)) return(NULL)
    
    target_col <- c("target_symbol","target.gene","target")[
      c("target_symbol","target.gene","target") %in% names(df)][1]
    if (is.na(target_col)) return(NULL)
    
    df <- df[!is.na(df[[target_col]]) & nzchar(df[[target_col]]),
             c("mature_mirna_id", target_col, "database")]
    colnames(df) <- c("miRNA","SYMBOL","database")
    df <- df[df$SYMBOL %in% universe_set, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    
    map_list <- split(unique(df$SYMBOL), df$miRNA)
    res <- lapply(names(map_list), function(m) {
      tgt <- map_list[[m]]
      a <- length(intersect(tgt, gene_set))
      b <- length(setdiff(tgt, gene_set))
      c <- length(setdiff(gene_set, tgt))
      d <- length(setdiff(universe_set, union(tgt, gene_set)))
      if ((a+b)==0 || (c+d)==0) return(NULL)
      ft <- suppressWarnings(fisher.test(matrix(c(a,b,c,d), 2, byrow=TRUE), alternative="greater"))
      data.frame(miRNA=m, p.value=ft$p.value, targets_in_deg=a, targets_total=a+b, stringsAsFactors=FALSE)
    })
    res <- dplyr::bind_rows(res); if (is.null(res) || !nrow(res)) return(NULL)
    res$FDR <- p.adjust(res$p.value, "BH"); res$method <- which_table
    res[order(res$FDR, -res$targets_in_deg), ]
  }
  
  tab_val <- enrich_mirna(sig_genes, universe_symbols, "validated")
  tab_pre <- enrich_mirna(sig_genes, universe_symbols, "predicted")
  enr_main <- dplyr::bind_rows(tab_val, tab_pre)
  if (!is.null(enr_main) && nrow(enr_main)) {
    write.csv(enr_main, file.path(out_dir,"miRNA_enrichment_main.csv"), row.names = FALSE)
    top15 <- head(enr_main[order(enr_main$FDR, -enr_main$targets_in_deg), ], 15)
    if (nrow(top15)) {
      p <- ggplot(top15, aes(x=reorder(miRNA, -log10(FDR)), y=-log10(FDR), fill=method)) +
        geom_col() +
        coord_flip() +
        labs(title = "miRNA Target Enrichment (BH, universe = all tested genes)",
             x = "miRNA", y = expression(-log[10]*"FDR")) +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir,"miRNA_enrichment_top15.png"), p, width=8, height=6, dpi=300)
    }
  } else {
    message("miRNA enrichment: no signal (check DEG size / universe).")
  }
  
  if (!is.null(tab_val) && nrow(tab_val)) {
    best <- head(tab_val$miRNA[order(tab_val$FDR)], 5)
    #    mm_net <- tryCatch(
    #      multiMiR::get_multimir(org="hsa", target = sig_genes, table = "validated", summary = FALSE),
    #      error = function(e) NULL
    #    )
    
    
    mm_net_df <- safe_multimir_table(sig_genes, "validated")
    if (!is.null(mm_net_df)) {
      target_col <- c("target_symbol","target.gene","target")[
        c("target_symbol","target.gene","target") %in% names(mm_net_df)][1]
      if (!is.na(target_col)) {
        net <- mm_net_df[mm_net_df$mature_mirna_id %in% best &
                           mm_net_df[[target_col]] %in% sig_genes,
                         c("mature_mirna_id", target_col, "database")]
        colnames(net) <- c("miRNA","SYMBOL","database")
        if (nrow(net)) write.csv(net, file.path(out_dir,"network_top5_miRNA_to_DEG.csv"), row.names = FALSE)
      }
    }
    
    if (!is.null(mm_net) && nrow(mm_net@data)) {
      df <- mm_net@data
      target_col <- c("target_symbol","target.gene","target")[
        c("target_symbol","target.gene","target") %in% names(df)][1]
      if (!is.na(target_col)) {
        net <- df[df$mature_mirna_id %in% best & df[[target_col]] %in% sig_genes,
                  c("mature_mirna_id", target_col, "database")]
        colnames(net) <- c("miRNA","SYMBOL","database")
        if (nrow(net)) write.csv(net, file.path(out_dir,"network_top5_miRNA_to_DEG.csv"), row.names = FALSE)
      }
    }
  }
  
  sig_genes_s <- NULL
  if (exists("pheno_aligned") && (("is_melas" %in% names(pheno_aligned)) || ("specialcase" %in% names(pheno_aligned)))) {
    mel_flag <- if ("is_melas" %in% names(pheno_aligned)) pheno_aligned$is_melas else pheno_aligned$specialcase
    mel_flag <- tolower(as.character(mel_flag)) %in% c("1","true","t","yes","y")
    if (any(mel_flag)) {
      if (exists("deg_sens") && "SYMBOL" %in% names(deg_sens)) {
        sig_genes_s <- unique(deg_sens$SYMBOL)
      } else if (exists("expr_cb")) {
        keep_idx <- which(!mel_flag)
        if (length(keep_idx) >= 4) {
          expr_s <- expr_cb[, keep_idx, drop=FALSE]
          ph_s   <- droplevels(pheno_aligned[keep_idx, ])
          des_s  <- model.matrix(~ 0 + group, data = transform(ph_s, group=factor(group, levels=c("Control","Pompe"))))
          colnames(des_s) <- levels(factor(ph_s$group, levels=c("Control","Pompe")))
          fit_s  <- limma::lmFit(expr_s, des_s)
          cont_s <- limma::makeContrasts(Pompe - Control, levels = des_s)
          fit_s2 <- limma::eBayes(limma::contrasts.fit(fit_s, cont_s))
          tt_s   <- limma::topTable(fit_s2, number=Inf, adjust.method="BH")
          tt_s$SYMBOL <- tt_s$SYMBOL %||% NA_character_
          if (!"SYMBOL" %in% names(tt_s)) {
          }
          if ("SYMBOL" %in% names(tt_s)) {
            tt_s$abs_logFC <- abs(tt_s$logFC)
            tt_sg <- tt_s[!is.na(tt_s$SYMBOL) & nzchar(tt_s$SYMBOL), ]
            tt_sg <- tt_sg[order(tt_sg$SYMBOL, -tt_sg$abs_logFC), ]
            tt_sg <- tt_sg[!duplicated(tt_sg$SYMBOL), ]
            sig_genes_s <- unique(tt_sg$SYMBOL[tt_sg$adj.P.Val < FDR_THRESH])
          }
        }
      }
    }
  }
  
  
  if (!is.null(sig_genes_s) && length(sig_genes_s)) {
    tab_val_s <- enrich_mirna(sig_genes_s, universe_symbols, "validated")
    tab_pre_s <- enrich_mirna(sig_genes_s, universe_symbols, "predicted")
    enr_sens  <- dplyr::bind_rows(tab_val_s, tab_pre_s)
    if (!is.null(enr_sens) && nrow(enr_sens)) {
      write.csv(enr_sens, file.path(out_dir,"miRNA_enrichment_noMELAS.csv"), row.names = FALSE)
      
      top_main <- head(enr_main$miRNA[order(enr_main$FDR)], 10)
      top_sens <- head(enr_sens$miRNA[order(enr_sens$FDR)], 10)
      ovlp     <- intersect(top_main, top_sens)
      writeLines(c(
        sprintf("Top10 overlap (validated+predicted): %d", length(ovlp)),
        paste("Common:", paste(ovlp, collapse=", "))
      ), file.path(out_dir,"_sensitivity_overlap.txt"))
    }
  }
  
  writeLines(c(capture.output(sessionInfo())), file.path(out_dir,"session_info.txt"))
  
  message("miRNA enrichment (FAIR/offline) tamamlandı. Çıktılar: ", out_dir)
  
  # ======================        miRNA VIS (multiMiR; FAIR; ayrı DB'ler)       ======================
  
  mir_dir <- file.path(root_dir, "results", "miRNA_analysis")
  dir.create(mir_dir, recursive = TRUE, showWarnings = FALSE)
  
  writeLines("WARNING: small sample size (9 Pompe vs 10 Control). Interpret with caution.",
             file.path(mir_dir, "_NOTE_small_n.txt"))
  
  sig_syms <- if (exists("deg")) head(unique(deg$SYMBOL), 50) else character(0)
  dfv <- data.frame()
  dfp <- data.frame()
  if (length(sig_syms) < 2) {
    message("[miRNA] Yeterli anlamlı gen yok; multiMiR adımları atlandı.")
  } else {
    
    mm_val <- multiMiR::get_multimir(org = "hsa",
                                     target = sig_syms,
                                     table = "validated",
                                     summary = FALSE)
    
    class(mm_val); isS4(mm_val); slotNames(mm_val)
    
    mm_to_df <- function(x) {
      if (is.null(x)) return(NULL)
      if (isS4(x) && "data" %in% slotNames(x)) as.data.frame(x@data) else
        if (is.data.frame(x)) x else NULL
    }
    
    dfv <- mm_val@data
    pcol_v <- intersect(c("p.value","p_value"), names(dfv))[1]
    if (!is.na(pcol_v)) {
      dfv$FDR <- p.adjust(dfv[[pcol_v]], method = "BH")
    }
    
    pltv <- dfv |>
      dplyr::count(database, mature_mirna_id, name="Target_Count") |>
      dplyr::group_by(database) |>
      dplyr::slice_max(Target_Count, n = 10, with_ties = FALSE) |>
      dplyr::ungroup()
    
    p_val <- ggplot(pltv, aes(x = mature_mirna_id, y = Target_Count, fill = database)) +
      geom_bar(stat="identity", position = position_dodge(width=.7), width=.6) +
      coord_flip() + theme_minimal(base_size = 12) +
      labs(title="Validated miRNA Targets (multiMiR; BH‑FDR applied)",
           x="miRNA", y="Target count", fill="DB")
    ggsave(file.path(mir_dir, "Validated_miRNA_targets.png"), p_val, width=8, height=6, dpi=300)
    write.csv(dfv, file.path(mir_dir, "validated_full_with_FDR.csv"), row.names = FALSE)
  }  
  else {
    message("[miRNA] Validated sonuç yok ya da erişilemedi.")
  }
  
  mm_pred <- tryCatch(
    multiMiR::get_multimir(org="hsa", target=sig_syms, table="predicted", summary=TRUE),
    error = function(e) NULL
  )
  if (!is.null(mm_pred) && nrow(mm_pred@data) > 0) {
    dfp <- mm_pred@data
    dfp$FDR <- p.adjust(dfp$p_value, method = "BH")
    
    pltp <- dfp |>
      dplyr::count(database, mature_mirna_id, name="Target_Count") |>
      dplyr::group_by(database) |>
      dplyr::slice_max(Target_Count, n = 10, with_ties = FALSE) |>
      dplyr::ungroup()
    
    p_pred <- ggplot(pltp, aes(x = mature_mirna_id, y = Target_Count, fill = database)) +
      geom_bar(stat="identity", position=position_dodge(width=.7), width=.6) +
      coord_flip() + theme_minimal(base_size = 12) +
      labs(title="Predicted miRNA Targets (multiMiR; BH‑FDR applied)",
           x="miRNA", y="Target count", fill="DB")
    ggsave(file.path(mir_dir, "Predicted_miRNA_targets.png"), p_pred, width=8, height=6, dpi=300)
    write.csv(dfp, file.path(mir_dir, "predicted_full_with_FDR.csv"), row.names = FALSE)
  } else {
    message("[miRNA] Predicted sonuç yok ya da erişilemedi.")
  }
  
  # ==================== Drug discovery from miRNA targets ====================
  gene_targets <- unique(c(
    if (exists("dfv") && "target_symbol" %in% names(dfv)) dfv$target_symbol else character(),
    if (exists("dfp") && "target_symbol" %in% names(dfp)) dfp$target_symbol else character()
  ))
  if (length(gene_targets) > 0) {
    drug_dir <- file.path(mir_dir, "drug_discovery")
    dir.create(drug_dir, recursive = TRUE, showWarnings = FALSE)
    try({
      dg <- rDGIdb::queryDGIdb(gene_targets)
      drug_tbl <- as.data.frame(dg@matches)
      write.csv(drug_tbl, file.path(drug_dir, "candidate_drugs.csv"), row.names = FALSE)
    }, silent = TRUE)
  }
}

#================= Drug discovery for miRNA targets (rDGIdb) ==================

expr2 <- na.omit(expr_use)
expr2 <- expr2[!duplicated(rownames(expr2)), ]
set.seed(123)
um <- umap::umap(t(expr2), n_neighbors = 8, random_state = 123)
ump_df <- data.frame(Sample = rownames(um$layout),
                     UMAP1 = um$layout[,1], UMAP2 = um$layout[,2],
                     Group = groups)
p_umap <- ggplot(ump_df, aes(UMAP1, UMAP2, color = Group, label = Sample)) +
  geom_point(size=3) + ggrepel::geom_text_repel(size=3, max.overlaps = 15) +
  theme_minimal(base_size = 13) + labs(title="UMAP QC (batch/SVA sonrası ifade)")
ggsave(file.path(mir_dir, "QC_UMAP.png"), p_umap, width=7, height=5, dpi=300)

meta_mir <- list(dataset="GSE38680", n_pompe=sum(groups=="Pompe"),
                 n_control=sum(groups=="Control"), adj="BH",
                 note_small_n=TRUE, date=as.character(Sys.Date()))
jsonlite::write_json(meta_mir, file.path(mir_dir, "_analysis_meta.json"), pretty=TRUE)



## ==== MODEL SAĞLIK KONTROL BLOĞU ====
## Varsayımlar: expr_cc (veya expr_noM), ph_cc (veya ph_noM), design, cont (contrasts) değişkenlerin var.

# 0) Hangi objeleri kullanıyoruz? (adlarını kendi akışına göre ayarla)
X_expr   <- if (exists("expr_cc")) expr_cc else expr_noM
PH       <- if (exists("ph_cc"))   ph_cc   else ph_noM
DESIGN   <- design
CONTRAST <- tryCatch(cont, error = function(e) NULL)

cat(">> Boyutlar: expr = ", nrow(X_expr), "x", ncol(X_expr),
    " | design = ", nrow(DESIGN), "x", ncol(DESIGN), "\n", sep="")

# 1) Hizalama ve boyut eşitliği
stopifnot(identical(colnames(X_expr), rownames(PH)))
stopifnot(nrow(DESIGN) == ncol(X_expr))

# 2) Grup seviyesi kontrolü
if (!"group" %in% names(PH)) stop("PH içinde 'group' yok.")
tbl_grp <- table(PH$group, useNA="ifany")
cat(">> Grup sayımları:\n"); print(tbl_grp)
stopifnot(length(setdiff(names(tbl_grp), NA)) >= 2, all(tbl_grp[names(tbl_grp)!= ""] > 0))

# 3) Design kalitesi
cat(">> Design sütunları:\n"); print(colnames(DESIGN))
if (anyNA(DESIGN) || any(is.infinite(DESIGN))) stop("Design içinde NA/Inf var.")
rk <- qr(DESIGN)$rank
cat(">> Design rank = ", rk, " / ", ncol(DESIGN), "\n", sep="")
if (rk < ncol(DESIGN)) warning("Design tam-rank değil (koliner kovaryat olabilir).")

# 4) Artık serbestlik derecesi (residual d.f.) > 0 mı?
rdf <- ncol(X_expr) - rk
cat(">> Residual d.f. = ", rdf, "\n", sep="")
stopifnot(rdf > 0)

# 5) Kontrast kontrolü
if (is.null(CONTRAST)) {
  stop("CONTRAST nesnesi bulunamadı. 'makeContrasts(...)' üretildi mi?")
} else {
  cat(">> Contrast OK. Boyut: ", nrow(CONTRAST), "x", ncol(CONTRAST), "\n", sep="")
  # contrast'ta geçen sütunlar design içinde var mı?
  need_cols <- colnames(CONTRAST)
  if (!all(need_cols %in% colnames(DESIGN))) {
    stop("Kontrastta geçen sütun(lar) design'da yok: ",
         paste(setdiff(need_cols, colnames(DESIGN)), collapse=", "))
  }
}

# 6) Mini-fit ile hızlı duman testi
suppressMessages({
  fit_test  <- limma::lmFit(X_expr, DESIGN)
  fit_test2 <- limma::eBayes(limma::contrasts.fit(fit_test, CONTRAST))
})
stopifnot(inherits(fit_test2, "MArrayLM"))
cat(">> Mini-fit başarılı. İlk 3 p-değeri:\n")
print(head(fit_test2$p.value[,1], 3))

cat("==> Model doğrulamaları TAMAM. Analize güvenle devam edebilirsin. ✅\n")


# ===================== Drug–gene sensitivity analysis with PharmacoGx =====================
if (requireNamespace("PharmacoGx", quietly = TRUE)) {
  library(PharmacoGx)
  gdsc <- downloadPSet("GDSC")
  cor_res <- PharmacoGx::drugGeneSensitivityCor(
    pSet = gdsc,
    drugs = drugNames(gdsc),
    genes = "GAA",
    mDataType = "rna"
  )
  print(head(cor_res))
}

-------------------------------------------------------------------------------
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("PharmacoGx")

library(PharmacoGx)

pSet <- PharmacoSet(
  name        = "PipelinePsSet",
  molecularProfiles = list(
    rnaseq = rnaseq_matrix   # veya microarray, proteomic, vb.
  ),
  sensitivityInfo     = sensitivity_info_df,      # IC50/GI50 değerleri
  sensitivityRaw      = sensitivity_raw_array,
  sensitivityProfiles = sensitivity_profiles_df,
  drug                = drug_info_df,
  cell                = cell_info_df
)


# Pearson korelasyon örneği
corMat <- summarizeSensitivityProfiles(
  pSet,
  summary.stat = "auc_recomputed"
)

geneDrugCor <- correlateProfiles(
  sensitivity.profile = corMat,
  molecular.profile = rnaseq_matrix
)
-------------------------------------------------------------------------------
  
  
  # ---- Drug–gene sensitivity analysis (PharmacoGx) ----
if (requireNamespace("PharmacoGx", quietly = TRUE)) {
  library(PharmacoGx)
  cat("\n>> PharmacoGx: obtaining GDSC2 pharmacogenomic dataset...\n")
  pset <- tryCatch(PharmacoGx::downloadPSet("GDSC2"), error = function(e) NULL)
  if (!is.null(pset)) {
    sens <- PharmacoGx::summariseSensitivityProfiles(pset, sensitivity.measure = "AAC")
    expr_pset <- t(PharmacoGx::exprs(pset))
    cells <- intersect(rownames(expr_pset), rownames(sens))
    genes <- intersect(rownames(X_expr), colnames(expr_pset))
    if (length(cells) > 1 && length(genes) > 0) {
      expr_sub <- expr_pset[cells, genes, drop = FALSE]
      sens_sub <- sens[cells, , drop = FALSE]
      cor_mat  <- cor(expr_sub, sens_sub, use = "pairwise.complete.obs")
      out_dir <- file.path(root_dir, "results", "pharmacogx")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(cor_mat, file.path(out_dir, "drug_gene_correlations.csv"))
      cat(">> PharmacoGx correlation matrix written to results/pharmacogx\n")
    } else {
      warning("Insufficient overlap between PSet and expression matrix; skipping PharmacoGx analysis.")
    }
  } else {
    warning("downloadPSet('GDSC2') failed; skipping PharmacoGx analysis.")
  }
} else {
  warning("Package 'PharmacoGx' not installed; skipping drug–gene sensitivity step.")
}







# Girdi: 'deg' (SYMBOL sütunu var), 'root_dir'
library(jsonlite); library(httr)

out_dir <- file.path(root_dir, "results", "drug")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- unique(na.omit(deg$SYMBOL))
chunk <- 150L
chunks <- split(genes, ceiling(seq_along(genes)/chunk))

all_hits <- list()
for (i in seq_along(chunks)) {
  q <- paste(chunks[[i]], collapse = ",")
  url <- paste0("https://dgidb.org/api/v2/interactions?genes=", URLencode(q))
  res <- try(GET(url, timeout(60)), silent = TRUE)
  if (inherits(res, "try-error") || http_error(res)) next
  js  <- content(res, as="text", encoding="UTF-8")
  x   <- fromJSON(js, flatten = TRUE)
  if (!length(x$matchedTerms)) next
  
  # matchedTerms$interactions altında drug–gene etkileşimleri var
  for (k in seq_len(nrow(x$matchedTerms))) {
    mt <- x$matchedTerms[k,]
    if (!length(mt$interactions[[1]])) next
    df <- mt$interactions[[1]][, c("drugName","interactionTypes","geneName",
                                   "sources","pmids","score","interactionDirection"), drop=FALSE]
    df$queryGene <- mt$searchTerm
    all_hits[[length(all_hits)+1]] <- df
  }
}
dgi <- if (length(all_hits)) do.call(rbind, all_hits) else data.frame()
if (nrow(dgi)) {
  # Temizlik
  dgi$interactionTypes <- vapply(dgi$interactionTypes, function(z) paste(z, collapse=";"), "")
  dgi$sources          <- vapply(dgi$sources,          function(z) paste(z, collapse=";"), "")
  dgi$pmids            <- vapply(dgi$pmids,            function(z) paste(z, collapse=";"), "")
  # Tekilleştir & skorla
  dgi <- unique(dgi)
  write.csv(dgi, file.path(out_dir, "drug_candidates_DGIdb.csv"), row.names = FALSE)
  message("DGIdb: ", nrow(dgi), " etkileşim yazıldı -> drug_candidates_DGIdb.csv")
} else {
  message("DGIdb'den kayıt dönmedi (gen listesi çok küçük olabilir).")
}


#Not: rDGIdb paketine gerek yok; REST ile hallediyoruz. İstersen sonra dgi tablosunu up/down DEG bilgisiyle birleştirip “up-gen hedefleyen inhibitörler” gibi öncelik kuralları ekleyebilirsin.




library(msigdbr); library(clusterProfiler); library(dplyr); library(stringr)

# 1) Up/Down setleri
up_gen   <- unique(deg$SYMBOL[deg$logFC > 0])
down_gen <- unique(deg$SYMBOL[deg$logFC < 0])

# 2) LINCS L1000 kimyasal imzaları (up/down ayrı gene set’ler)
msig_lincs <- msigdbr(species = "Homo sapiens",
                      category = "C2",
                      subcategory = "CP:LINCS_L1000") %>%
  select(gs_name, gene_symbol)

# Ayrı koleksiyonlara böl
lincs_up   <- msig_lincs %>% filter(str_ends(gs_name, "_UP"))
lincs_down <- msig_lincs %>% filter(str_ends(gs_name, "_DN"))

# 3) Reversal mantığı:
#  - up_gen  ~ LINCS *_DN ile zenginleşsin
#  - down_gen~ LINCS *_UP ile zenginleşsin
en_up_dn   <- enricher(up_gen,   TERM2GENE = lincs_down, pAdjustMethod = "BH")
en_dn_up   <- enricher(down_gen, TERM2GENE = lincs_up,   pAdjustMethod = "BH")

# 4) Aday bileşik adını temizle (gs_name: DRUGNAME_<dose>_DN/UP gibi)
clean_name <- function(x) toupper(gsub("(_UP|_DN)$","",gsub("^LINCS_L1000_","",x)))
get_top <- function(enr, n=30) {
  if (is.null(enr) || nrow(as.data.frame(enr))==0) return(NULL)
  df <- as.data.frame(enr)
  df$drug <- clean_name(df$ID)
  df %>% arrange(p.adjust) %>% group_by(drug) %>% slice_min(order_by = p.adjust, n = 1) %>% ungroup()
}

cand_rev <- bind_rows(
  get_top(en_up_dn, n=50) %>% mutate(direction="UPvsDN"),
  get_top(en_dn_up, n=50) %>% mutate(direction="DNvsUP")
)

dir.create(file.path(root_dir,"results","drug"), showWarnings = FALSE, recursive = TRUE)
if (!is.null(cand_rev) && nrow(cand_rev)) {
  write.csv(cand_rev,
            file.path(root_dir,"results","drug","drug_candidates_LINCS_reversal.csv"),
            row.names = FALSE)
  message("LINCS reversal adayları yazıldı -> drug_candidates_LINCS_reversal.csv")
}


#İpuçları:
#• p.adjust küçük ve GeneRatio büyük olanlar daha umut vericidir.
#• Aynı ilaç birden fazla doz/koşulda çıkabilir; kod en iyi p-değerli satırı tekilleştiriyor.
#• İstersen çıktıdan up/down tutarlılığı iste (ilaç hem UPvsDN’de hem DNvsUP’de anlamlıysa +1 puan gibi).






#Pratikte ikisinin kesişimi en iyi önceliklendirme listesi olur.
#Kesişim/önceliklendirme (opsiyonel)


# dgi (DGIdb) ve cand_rev (LINCS) hazırsa:
if (exists("dgi") && nrow(dgi) && exists("cand_rev") && nrow(cand_rev)) {
  dgi$drug_std  <- toupper(dgi$drugName)
  cand_rev$drug_std <- toupper(cand_rev$drug)
  pri <- inner_join(cand_rev, distinct(dgi, drug_std), by = "drug_std")
  write.csv(pri, file.path(root_dir,"results","drug","drug_candidates_intersection.csv"), row.names = FALSE)
  message("Kesişim listesi yazıldı -> drug_candidates_intersection.csv (daha güçlü adaylar)")
}
