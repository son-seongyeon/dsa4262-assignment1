suppressPackageStartupMessages(library(bambu))
se <- readRDS('bambu_output/se.rds')

sink('bambu_output/analysis_answers.txt')

cat('=== 4.1 show(se) ===\n')
show(se)

cat('\n\n=== 4.2a novel transcripts ===\n')
print(table(mcols(se)$novelTranscript))

cat('\n\n=== 4.2b top 5 expressed genes (by total CPM across samples) ===\n')
se_gene <- transcriptToGeneExpression(se)
gene_cpm <- assays(se_gene)$CPM
gene_total <- rowSums(gene_cpm)
top5 <- head(sort(gene_total, decreasing=TRUE), 5)
print(top5)

cat('\n\n=== 4.2c transcripts per gene (all annotated+novel) ===\n')
geneid <- rowData(se)$GENEID
tx_per_gene <- table(geneid)
cat('min:', min(tx_per_gene), ' max:', max(tx_per_gene), ' mean:', mean(tx_per_gene), '\n')

cat('\n\n=== 4.2d transcripts per gene (expressed, >=10 reads in at least one sample) ===\n')
counts <- assays(se)$counts
expressed <- rowSums(counts >= 10) > 0
geneid_expr <- geneid[expressed]
tx_per_gene_expr <- table(geneid_expr)
cat('n expressed transcripts:', sum(expressed), '\n')
cat('min:', min(tx_per_gene_expr), ' max:', max(tx_per_gene_expr), ' mean:', mean(tx_per_gene_expr), '\n')

cat('\n\n=== 4.2e single-exon novel transcripts ===\n')
novel <- mcols(se)$novelTranscript
exon_counts <- elementNROWS(rowRanges(se))
cat('total novel transcripts:', sum(novel), '\n')
cat('single-exon novel transcripts:', sum(novel & exon_counts == 1), '\n')

cat('\n\n=== 4.2f highest expressed transcript among unique-read-count-only transcripts ===\n')
uniqueCounts <- assays(se)$uniqueCounts
fullCounts <- assays(se)$fullLengthCounts
totalCounts <- assays(se)$counts
has_ambiguity <- rowSums(totalCounts) != rowSums(uniqueCounts)
unique_only <- !has_ambiguity & rowSums(uniqueCounts) > 0
cat('n transcripts with unique-only counts:', sum(unique_only), '\n')
if (sum(unique_only) > 0) {
  uc_sub <- rowSums(uniqueCounts[unique_only,,drop=FALSE])
  top_unique <- sort(uc_sub, decreasing=TRUE)[1]
  print(top_unique)
}

sink()
cat('saved to bambu_output/analysis_answers.txt\n')
