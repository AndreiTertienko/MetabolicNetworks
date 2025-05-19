library(ALDEx2)
library(readr)

df <- read_tsv("C:/Users/A/Documents/MSc_Thesis/output/absolute_abundance_matrix.tsv", col_names = TRUE)
#df <- read.table("C:/Users/A/Documents/MSc_Thesis/output/absolute_abundance_matrix.tsv", row.names=1)
df_df <- as.data.frame(df)
rownames(df_df) <- df_df[,1]
df_df <- df_df[,-1]

samples = colnames(df_df)


sample_cohorts <- vector()
patients <- vector()
day0 <- 0
day84 <- 0
for (sample in samples) {
    patient = unlist(strsplit(sample, "_"))[2]
    patients <- append(patients, patient)
    if (grepl("Visit1", sample, fixed = TRUE)) {
        day0 <-  day0 + 1
        #sample_cohorts <- append(sample_cohorts, "day0")
    } else if (grepl("Visit3", sample, fixed = TRUE)) {
        #sample_cohorts <- append(sample_cohorts, "day84")
        day84 <- day84 + 1
    }
}

i <- 0
k <- 0
while (i < day0 | k < day84) {
    if (i < day0){
        sample_cohorts <- append(sample_cohorts, "day0")
        i <- i + 1
    }
    if (k < day84){
        sample_cohorts <- append(sample_cohorts, "day84")
        k <- k + 1
    }
}


# DESeq2 analysis
library(DESeq2)
sample_info <- data.frame(samples, patients, sample_cohorts)
#DESeq2_cohorts <- data.frame(samples, sample_cohorts)
#sample_patients <- data.frame(samles, patients)
names(sample_info) <- c('samples', 'patient_id', 'cohort')
rownames(sample_info) <- sample_info[,1]
sample_info[1] <- NULL
all(rownames(sample_info) %in% colnames(df_df))
dds <- DESeqDataSetFromMatrix(countData = df_df, colData = sample_info, design = ~ patient_id + cohort)

featureData <- data.frame(gene=rownames(df_df))
mcols(dds) <- DataFrame(mcols(dds), featureData)

#dds <- dds[keep,]

dds <- DESeq(dds, sfType = "poscounts", test = "Wald")
res05 <- results(dds, alpha=0.05)
resOrdered <- res05[order(res05$pvalue),]
summary(res05)
deseq_out <- data.frame(microbe = rownames(resOrdered), resOrdered, row.names = NULL)
write.table(deseq_out, file = "DA_results/DESeq2.tsv", sep="\t", quote = FALSE, row.names = FALSE)


# edgeR analysis
library(edgeR)
dgList <- DGEList(counts=df_df, genes=rownames(df_df), group=sample_cohorts)
z = calcNormFactors(dgList)
design = model.matrix(~ sample_cohorts + patients)
dge <- estimateDisp(z, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit,coef = 2)
edge.res <- topTags(lrt, n = Inf)
edge.res.all <- edge.res$table %>% as.data.frame()
colnames(edge.res.all)[1] = "Microbe"

write.table(edge.res.all, file='DA_results/edgeR_adj.tsv', quote=F, sep="\t", row.names = FALSE)


process_network <- function(network_file, 
                            outname,
                            deseq_file = "DA_results/DESeq2.tsv",
                            edger_file = "DA_results/edgeR.tsv") {
    # Read the network TSV (no header)
    net <- read.csv(network_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(net) <- c("Vertex.A", "Vertex.B")
    
    # Add Type columns: if value starts with "s__", then microbe, else metabolite
    net$Type.A <- ifelse(grepl("^d__", net$Vertex.A), "microbe", "metabolite")
    net$Type.B <- ifelse(grepl("^d__", net$Vertex.B), "microbe", "metabolite")
    
    # Create a 'microbe' column: use Vertex.A if it's a microbe; otherwise, use Vertex.B
    net$microbe <- ifelse(net$Type.A == "microbe", net$Vertex.A,
                          ifelse(net$Type.B == "microbe", net$Vertex.B, NA))
    
    # Read p-value files
    aldex <- read_tsv(aldex_file, col_names = TRUE)
    deseq <- read_tsv(deseq_file, col_names = TRUE)
    edger <- read_tsv(edger_file, col_names = TRUE)
    
    # Rename the aldex p-value columns
    names(aldex)[names(aldex) == "we.ep"] <- "ALDEx2.Welch"
    names(aldex)[names(aldex) == "wi.ep"] <- "ALDEx2.Wilcox"

    # Merge using the microbe column.
    merged_net <- merge(net, aldex, by = "microbe", all.x = TRUE)

    # Rename the edger p-value column
    edger_subset <- edger[, c("genes", "PValue")]
    names(edger_subset)[names(edger_subset) == "PValue"] <- "edgeR.Exact" 
    names(edger_subset)[names(edger_subset) == "genes"] <- "microbe" 
   
    merged_net <- merge(merged_net, edger_subset, by = "microbe", all.x = TRUE)
    
    # Only use the pvalue column from DESeq2 and rename it to DESeq2.Wald.
    deseq_subset <- deseq[, c("microbe", "pvalue")]
    names(deseq_subset)[names(deseq_subset) == "pvalue"] <- "DESeq2.Wald"
    
    merged_net <- merge(merged_net, deseq_subset, by = "microbe", all.x = TRUE)
    
    # Save the merged network to file
    write.table(merged_net, file = outname, sep = "\t", quote = FALSE, row.names = FALSE)
    return(merged_net)
}


basepath <- "C:/Users/A/Documents/MSc_Thesis/"
cohorts <- c("merged_cohorts")
zeros <- c("zeros", "nozeros")
thresholds <- c("0.85", "0.90", "0.95")
filepaths <- vector()
outnames <- vector()

for (cohort in cohorts) {
    for (zero in zeros) {
        for (threshold in thresholds) {
            filename <- paste0("network_", threshold, "_", zero, ".csv")
            filepath <- file.path(basepath, "merged_fva", cohort, filename)
            filepaths <- append(filepaths, filepath)
            outname <- paste0("DA_results/", "DA_", filename)
            outnames <- append(outnames, outname)
}}}

mapply(process_network, filepaths, outnames)
