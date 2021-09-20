library(DESeq2)

args <- commandArgs(TRUE)
Type <- args[1]
input_dir <- args[2]
output_dir <- args[3]
meta_path <- args[4]

wd <- getwd()
input_matrix <- paste0('matrix_', Type, '.txt')
input_info <- meta_path
output_file <- paste0('diffexp_result_', Type, '.txt')


info_data <- read.table(input_info, header = TRUE, row.names = 1)
setwd(input_dir)
input_data <- read.table(input_matrix, header = TRUE, row.names = 1)
input_data <- round(input_data, digits = 0)
input_data <- as.matrix(input_data)


dds <- DESeqDataSetFromMatrix(countData = input_data, colData = info_data, design = ~Group)
dds <- tryCatch(
    {
        DESeq(dds)
    },
    error = function(e) {
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds <- nbinomWaldTest(dds)
        dds
    }
)
result <- results(dds, contrast = c(colnames(info_data), unique(info_data[, 1])), alpha = 0.05, pAdjustMethod = 'BH')
summary(result)
result <- result[order(result$padj), ]
result_data <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized = TRUE)), by = 'row.names', sort = FALSE)
names(result_data)[1] <- 'Gene'

setwd(wd)
setwd(output_dir)
write.table(result_data, file = output_file, sep = '\t', quote = F, row.names = F)
