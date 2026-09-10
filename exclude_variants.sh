awk 'NR > 1 && $1 == 11 && $3 >= 61293499 && $3 <= 61854782 {print $2}' \
ldblk_1kg_eur/snpinfo_1kg_hm3 \
> FADS_variants_GRCh37.txt
