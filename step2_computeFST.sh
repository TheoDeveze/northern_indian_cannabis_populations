#!/bin/bash

DIR="FST_outputs"
OUT_LONG="fst_long.tsv"

echo -e "dataset\treplicate\tcomparison\tfst" > "$OUT_LONG"

for f in "$DIR"/*_chrom.windowed.weir.fst
do
    filename=$(basename "$f")

    dataset=$(echo "$filename" | cut -d'_' -f1)

    K=$(echo "$filename" | grep -oE "K[0-9]+" | head -n1)

    replicate=$(echo "$filename" | grep -oE "rep[0-9]+" | sed 's/rep//')

    cluster1=$(echo "$filename" | grep -oE "cluster[0-9]+" | head -n1)
    cluster2=$(echo "$filename" | grep -oE "cluster[0-9]+" | tail -n1)

    threshold=$(echo "$filename" | grep -oE "(pur|majo)" | head -n1)

    if [[ -z "$dataset" || -z "$K" || -z "$cluster1" || -z "$cluster2" || -z "$threshold" || -z "$replicate" ]]; then
        echo "ERROR parsing: $filename" >&2
        continue
    fi

    c1=$(echo "$cluster1" | sed 's/cluster/c/')
    c2=$(echo "$cluster2" | sed 's/cluster/c/')

    comparison="${K}_${c1}_vs_${c2}_${threshold}"

    fst=$(awk '
        BEGIN {sum=0; wsum=0}
        NR>1 && $5 != "nan" {
            sum += $5 * $4
            wsum += $4
        }
        END {
            if (wsum > 0) printf "%.8f", sum/wsum;
            else print "NA";
        }
    ' "$f")

    echo -e "${dataset}\t${replicate}\t${comparison}\t${fst}" >> "$OUT_LONG"

done
