# --- Download and convert all data to BED format 
curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb #ENCODE cCRE data 
bigBedToBed encodeCcreCombined.bb cCREs_hg38.bed 

curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/fantom5/hg38.cage_peak.bb #FANTOM5 TSS peaks data (robust peaks across cell lines -- indicator of functional activity)
bigBedToBed hg38.cage_peak.bb TSS_peaks_hg38.bed 

curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb # MANE protein-coding gene coordinates 
bigBedtoBed mane.bb mane.bed 

# Trim and sort data 
bedtools sort -i 5utr_variants.bed > 5utr_variants.sorted.bed 
cut -f1-4 TSS_peaks_hg38.bed > TSS_peaks_trimmed.bed 
cut -f1-4,7,8,18-20 mane.bed > mane_trimmed.bed

# ---- Annotate 5'UTR variants with cCREs 
# check column layout and cut irrelevant fields. retaining variant coordinates (1-3), cCRE ID (7), ENCODE label/classification (13), and UCSC label (16)
bedtools intersect -a 5utr_variants.bed -b cCREs_hg38.bed -wo > 5utr_cCRE_overlap.bed
cut -f1-3,7,13,16,19-21 5utr_cCRE_overlap.bed > 5utr_cCRE_trimmed.bed 

# Collapse overlapping cCREs. Grouping on variant (1-3), collapse cCRE annotations
bedtools groupby -i 5utr_cCRE_trimmed.bed -g 1,2,3 -c 4,5,6 -o collapse,collapse,collapse > 5utr_cCRE_annotated.bed

# ----- TSS Annotation (direct overlap) - Does this variant directly disrupt a TSS peak?
# Annotate 5'UTR variants with trimmed FANTOM5 TSS peaks and collapse per-variant
bedtools intersect -a 5utr_variants.bed -b TSS_peaks_trimmed.bed -wo > 5utr_TSS_overlap.bed
bedtools groupby -i 5utr_TSS_overlap_sorted.bed -grp 1-3 -c 7 -o collapse > 5utr_TSS_intersect.bed 

# ----- TSS Annotation (proximal promoter check) - Does the variant disrupt nearby TSS peaks? How close are the peaks? 
# Variants near but not inside peaks that might disrupt regulatory regions just upstream/downstream of TSSs. can adjust -w but kept to default. 
bedtools window -a 5utr_variants.bed -b TSS_peaks_trimmed.bed -w 1000 > 5utr_TSS_window1kb.bed

# Nearest non-overlapping TSS. 
bedtools closest -a 5utr_variants.sorted.bed -b TSS_peaks_trimmed.bed -d -io > 5utr_TSS_closest.bed 

# ---- Gene mapping 
# Identify closest upstream/downstream gene + calculate nt distance (-D)
 bedtools closest -a 5utr_variants.sorted.bed -b mane_trimmed.bed -D a -id > 5utr_upstream_gene.bed 
 bedtools closest -a 5utr_variants.sorted.bed -b mane_trimmed.bed -D a -iu > 5utr_downstream_gene.bed 

