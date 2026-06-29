#!/bin/bash

module load StdEnv/2023
module load bcftools
module load vcftools

# =========================================================
# Parameters
# =========================================================

VCF_3D="renamed_3D_forFST.vcf"
VCF_LITT="renamed_litt_forFST.vcf"

NBOOT=1000
OUTDIR="FST_K7_DAPC_majo"

WINDOW_SMALL=100000
WINDOW_LARGE=1000000000

SEED=420

mkdir -p ${OUTDIR}


# =========================================================
# Reproducible subsampling
# =========================================================

subsample_file() {

    infile=$1
    n=$2
    seed=$3
    outfile=$4

    shuf --random-source=<(yes ${seed}) ${infile} \
        | head -n ${n} > ${outfile}
}


# =========================================================
# Compute FST
# =========================================================

compute_fst() {

    vcf=$1
    file1=$2
    file2=$3
    label=$4
    rep=$5

    prefix=${OUTDIR}/${label}_rep${rep}


    # 100 kb windows
    vcftools \
        --vcf ${vcf} \
        --weir-fst-pop ${file1} \
        --weir-fst-pop ${file2} \
        --fst-window-size ${WINDOW_SMALL} \
        --out ${prefix}_100kb \
        > /dev/null 2>&1


    # Whole chromosome FST
    vcftools \
        --vcf ${vcf} \
        --weir-fst-pop ${file1} \
        --weir-fst-pop ${file2} \
        --fst-window-size ${WINDOW_LARGE} \
        --out ${prefix}_chrom \
        > /dev/null 2>&1
}


# =========================================================
# Pairwise FST with equal sample size
# =========================================================

run_pairwise() {

    vcf=$1
    pop1=$2
    pop2=$3
    label=$4
    nmin=$5


    echo "Running ${label} with n=${nmin}"


    for ((rep=1; rep<=NBOOT; rep++))
    do

        sub1="tmp_${label}_${rep}_1.txt"
        sub2="tmp_${label}_${rep}_2.txt"


        # Reproducible but different seed for each replicate
        seed_rep=$((SEED + rep))


        subsample_file \
            ${pop1} \
            ${nmin} \
            ${seed_rep} \
            ${sub1}


        subsample_file \
            ${pop2} \
            ${nmin} \
            ${seed_rep} \
            ${sub2}


        compute_fst \
            ${vcf} \
            ${sub1} \
            ${sub2} \
            ${label} \
            ${rep}


        rm ${sub1} ${sub2}

    done
}


# =========================================================
# Main analysis
# Only DAPC K7 majority clusters
# =========================================================

VCF=${VCF_LITT}

files=(litt_K7_cluster*_majo_DAPC.txt)


if [ ${#files[@]} -lt 2 ]; then

    echo "Not enough clusters found"
    exit 1

fi


echo "Processing litt K7 DAPC majority clusters"


# =========================================================
# Determine minimum cluster size
# =========================================================

nmin_global=999999


for f in "${files[@]}"
do

    n=$(wc -l < ${f})


    if [ ${n} -lt ${nmin_global} ]; then

        nmin_global=${n}

    fi

done


echo "Global nmin = ${nmin_global}"


# =========================================================
# Pairwise comparisons
# =========================================================

for ((i=0; i<${#files[@]}; i++))
do

    for ((j=i+1; j<${#files[@]}; j++))
    do

        f1=${files[$i]}
        f2=${files[$j]}


        name1=$(basename ${f1} .txt)
        name2=$(basename ${f2} .txt)


        label="${name1}_VS_${name2}"


        run_pairwise \
            ${VCF} \
            ${f1} \
            ${f2} \
            ${label} \
            ${nmin_global}

    done

done
