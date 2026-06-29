#!/bin/bash

INPUT_DIR="PI_site_rarefaction"
OUTPUT="pi_site_rarefaction_summary.tsv"

echo -e "dataset\tK\tcluster\tassignment\treplicate\tmean_pi" > ${OUTPUT}

for file in ${INPUT_DIR}/*.sites.pi
do

    base=$(basename ${file})

    # =====================================================
    # parse metadata from filename
    # expected:
    # 3d_K2_cluster1_majo_rep1.sites.pi
    # 3d_K2_cluster1_pur_rep1.sites.pi
    # =====================================================

    dataset=$(echo ${base} | grep -oE "^(3d|litt)")

    K=$(echo ${base} | grep -oE "K[0-9]+" | grep -oE "[0-9]+")
    cluster=$(echo ${base} | grep -oE "cluster[0-9]+" | grep -oE "[0-9]+")
    rep=$(echo ${base} | grep -oE "rep[0-9]+" | grep -oE "[0-9]+")

    # pur vs majo (Q matrix > 0.7 or all samples)
    assignment=$(echo ${base} | grep -oE "(pur|majo)")

    mean_pi=$(awk '
        NR>1 && $NF != "nan" {
            sum += $NF;
            n++
        }
        END {
            if (n > 0) print sum/n;
            else print "NA"
        }
    ' ${file})

    echo -e "${dataset}\t${K}\t${cluster}\t${assignment}\t${rep}\t${mean_pi}" >> ${OUTPUT}

done
