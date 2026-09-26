suppressPackageStartupMessages(library(bambu))
se <- readRDS('bambu_output/se.rds')

novel <- mcols(se)$novelTranscript
counts <- assays(se)$counts
total_counts <- rowSums(counts)

novel_idx <- which(novel)
novel_counts <- total_counts[novel_idx]
names(novel_counts) <- rownames(se)[novel_idx]

top_novel <- sort(novel_counts, decreasing=TRUE)[1:10]
cat('=== Top 10 novel transcripts by total read count ===\n')
print(top_novel)

cat('\n=== Details for top 3 ===\n')
for (txid in names(top_novel)[1:3]) {
  idx <- which(rownames(se) == txid)
  gr <- rowRanges(se)[[idx]]
  gene <- rowData(se)$GENEID[idx]
  novelGene <- mcols(se)$novelGene[idx]
  n_exons <- length(gr)
  exon_widths <- width(gr)
  total_width <- sum(exon_widths)
  chr <- as.character(seqnames(gr)[1])
  start_pos <- min(start(gr))
  end_pos <- max(end(gr))
  strand_val <- as.character(strand(gr)[1])
  cat(sprintf('\ntranscript: %s | gene: %s | novelGene: %s\n', txid, gene, novelGene))
  cat(sprintf('  location: chr%s:%d-%d (%s strand)\n', chr, start_pos, end_pos, strand_val))
  cat(sprintf('  n_exons: %d | total_exonic_length: %d bp\n', n_exons, total_width))
  cat(sprintf('  total reads (sum across samples): %.1f\n', total_counts[idx]))
  cat(sprintf('  per-sample counts: %s\n', paste(round(counts[idx,],1), collapse=', ')))
}
