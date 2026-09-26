suppressPackageStartupMessages(library(bambu))
se <- readRDS('bambu_output/se.rds')

sink('verify_task4_output.txt')

cat('=== dim(se) ===\n')
print(dim(se))

cat('\n=== novel transcript count (mcols novelTranscript) ===\n')
novel <- mcols(se)$novelTranscript
print(table(novel))

cat('\n=== single-exon novel transcript count ===\n')
exon_counts <- elementNROWS(rowRanges(se))
cat('novel & single-exon:', sum(novel & exon_counts==1), '\n')
cat('distribution of exon counts among novel transcripts:\n')
print(table(exon_counts[novel]))

cat('\n=== gene-level CPM top5 (method: transcriptToGeneExpression then sum CPM across samples) ===\n')
se_gene <- transcriptToGeneExpression(se)
gene_cpm <- assays(se_gene)$CPM
gene_total <- rowSums(gene_cpm)
print(head(sort(gene_total, decreasing=TRUE), 8))
cat('per-sample CPM for top gene (sanity check no single sample dominates unfairly):\n')
top1 <- names(sort(gene_total, decreasing=TRUE))[1]
print(gene_cpm[top1,])

cat('\n=== transcripts per gene: all annotated+novel, denominator = all genes with >=1 transcript in se ===\n')
geneid <- rowData(se)$GENEID
tx_per_gene <- table(geneid)
cat('n genes (denominator):', length(tx_per_gene), '\n')
cat('min:', min(tx_per_gene), 'max:', max(tx_per_gene), 'mean:', mean(tx_per_gene), '\n')

cat('\n=== transcripts per gene: expressed only (>=10 reads in >=1 sample), same gene set restricted ===\n')
counts <- assays(se)$counts
expressed <- rowSums(counts >= 10) > 0
cat('n expressed transcripts:', sum(expressed), '\n')
geneid_expr <- geneid[expressed]
tx_per_gene_expr <- table(geneid_expr)
cat('n genes with >=1 expressed transcript (denominator):', length(tx_per_gene_expr), '\n')
cat('min:', min(tx_per_gene_expr), 'max:', max(tx_per_gene_expr), 'mean:', mean(tx_per_gene_expr), '\n')

cat('\n=== ENST00000331825 verification ===\n')
uniqueCounts <- assays(se)$uniqueCounts
totalCounts <- assays(se)$counts
idx <- which(rownames(se) == 'ENST00000331825')
cat('rowname found:', length(idx)>0, '\n')
cat('total counts per sample:', round(totalCounts[idx,],2), '\n')
cat('unique counts per sample:', round(uniqueCounts[idx,],2), '\n')
cat('sum unique:', sum(uniqueCounts[idx,]), '\n')
cat('sum total:', sum(totalCounts[idx,]), '\n')
cat('equal (fully unique, no ambiguity)?', isTRUE(all.equal(sum(uniqueCounts[idx,]), sum(totalCounts[idx,]))), '\n')

cat('\n=== recompute unique-only selection from scratch, verify top pick ===\n')
has_ambiguity <- abs(rowSums(totalCounts) - rowSums(uniqueCounts)) > 1e-6
unique_only <- (!has_ambiguity) & (rowSums(uniqueCounts) > 0)
cat('n transcripts fully-unique (no ambiguity), by strict equality check:', sum(unique_only), '\n')
uc_sub <- rowSums(uniqueCounts[unique_only,,drop=FALSE])
top3 <- sort(uc_sub, decreasing=TRUE)[1:3]
print(top3)

cat('\n=== 3 selected novel transcripts: expression totals cross-check ===\n')
for (txid in c('BambuTx381','BambuTx757','BambuTx369')) {
  i <- which(rownames(se)==txid)
  cat(txid, 'total counts per sample:', round(totalCounts[i,],2), ' SUM:', round(sum(totalCounts[i,]),2), '\n')
}

sink()
cat('done, see verify_task4_output.txt\n')
