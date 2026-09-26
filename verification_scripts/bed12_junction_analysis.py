"""
Task 4.4 verification: for each of the 3 selected novel transcript loci,
count how many reads in the submitted reads.bb track (decoded to BED12 with
`bigBedToBed`) overlap the locus, exactly match Bambu's predicted intron gap,
or instead span that region as one continuous block (no gap).

Input: reads_decoded.bed, produced with:
    bigBedToBed task4_4_ucsc_track_files/reads.bb reads_decoded.bed

Used to produce the junction-support table in Task 4.4 of Assignment1_Report.md.
"""

# (chrom, gap_start, gap_end) for Bambu's predicted intron at each locus,
# derived from the BED12 blockSizes/blockStarts of novel_transcripts.bb
loci = {
    'BambuTx381': ('chr11', 82689597, 82689645),   # 48 bp predicted intron
    'BambuTx757': ('chr10', 96750365, 96750433),   # 68 bp predicted intron
    'BambuTx369': ('chr10', 120355054, 120355114), # 60 bp predicted intron
}

# full transcript span per locus, used to decide whether a read "overlaps"
spans = {
    'BambuTx381': ('chr11', 82689529, 82689941),
    'BambuTx757': ('chr10', 96750264, 96750965),
    'BambuTx369': ('chr10', 120354640, 120355223),
}

results = {k: {'overlap': 0, 'exact_gap': 0, 'spans_gap_one_block': 0} for k in loci}

with open('reads_decoded.bed') as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        block_sizes = [int(x) for x in fields[10].rstrip(',').split(',')]
        block_starts = [int(x) for x in fields[11].rstrip(',').split(',')]
        blocks = [(start + bs, start + bs + sz) for bs, sz in zip(block_starts, block_sizes)]
        gaps = [(blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1)]

        for name, (lchrom, gstart, gend) in loci.items():
            schrom, sstart, send = spans[name]
            if chrom != schrom or end <= sstart or start >= send:
                continue
            results[name]['overlap'] += 1
            if any(gs == gstart and ge == gend for gs, ge in gaps):
                results[name]['exact_gap'] += 1
            if any(bs <= gstart and be >= gend for bs, be in blocks):
                results[name]['spans_gap_one_block'] += 1

for name, r in results.items():
    print(name, r)
