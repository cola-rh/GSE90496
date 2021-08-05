library(data.table)
mat = fread("GSE90496_beta.txt")
mat = as.data.frame(mat)
rownames(mat) = mat[, 1]
mat = mat[, grepl("SAMPLE", colnames(mat))]

colnames(mat) = gsub(" ", "_", colnames(mat))

tb = read.csv("GSE90496_sample.csv")

subtype1 = gsub("^(.*?),.*$", "\\1", tb$Title)
subtype2 = gsub(", sample.*$", "", tb$Title)
sample_id = gsub("^.*, sample\\s(\\d+).*$", "\\1", tb$Title)
sample_id = paste0("SAMPLE_", sample_id)


anno = data.frame(subtype1 = subtype1, subtype2 = subtype2)
rownames(anno) = sample_id
anno = anno[colnames(mat), ]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19",
    package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
probe = IlluminaHumanMethylation450kanno.ilmn12.hg19 # change to a short name

mat = mat[rownames(mat) %in% rownames(getAnnotation(probe, what = "Locations")), ]

l = getAnnotation(probe, what = "Locations")$chr %in% paste0("chr", 1:22) &
    is.na(getAnnotation(probe, what = "SNPs.137CommonSingle")$Probe_rs)
mat = mat[l, ]
cgi_anno = getAnnotation(probe, "Islands.UCSC")$Relation_to_Island[l]

mat = mat[cgi_anno == "Island", ]

saveRDS(list(mat = mat, anno = anno), file = "GSE90496_islands_processed.rds")


