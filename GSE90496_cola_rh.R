setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/GSE90496")
library(cola)

tumor_type = read.table("GSE90496_tumor_type.csv", sep = ";", header = TRUE, comment.char = "")
tumor_col = structure(unique(tumor_type$color), names = unique(tumor_type$tumor_type))
tumor_type = structure(tumor_type$tumor_type, names = tumor_type$subtype)

lt = readRDS("GSE90496_islands_processed.rds")
mat = as.matrix(lt$mat)
anno = lt$anno
anno$subtype2 = gsub("\\s+", " ", anno$subtype2)
anno$subtype2 = gsub("\\s+$", "", anno$subtype2)

anno = data.frame(meth_class = anno$subtype2, tumor_type = tumor_type[anno$subtype2])

anno_col = list(
    tumor_type = tumor_col
)

mat = adjust_matrix(mat)

library(matrixStats)
rh = hierarchical_partition(mat, cores = 4, 
    top_value_method = c("SD", "ATC"), max_k = 8, 
    partition_method = c("kmeans", "skmeans"),
    scale_rows = FALSE, anno = anno, top_n = 1000, 
    subset = 500, group_diff = 0.25, min_n_signatures = 1000,
    filter_fun = function(mat) {
        s = rowSds(mat)
        order(-s)[1:30000]
    })

saveRDS(rh, file = "GSE90496_cola_rh.rds")






cola_report(rh, output = "GSE90496_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'GSE90496'")
