#!/bin/bash

module load StdEnv/2023
module load bcftools
module load vcftools

# -----------------------------
# Paramètres
# -----------------------------
VCF_3D="renamed_3D_forFST.vcf"
VCF_LITT="renamed_litt_forFST.vcf"

NBOOT=1000
OUTDIR="FST_windowed_v2"

WINDOW_SMALL=100000
WINDOW_LARGE=1000000000

SEED=420

mkdir -p ${OUTDIR}

# -----------------------------
# Subsampling reproductible
# -----------------------------
subsample_file() {
    infile=$1
    n=$2
    seed=$3
    outfile=$4

    shuf --random-source=<(yes ${seed}) ${infile} | head -n ${n} > ${outfile}
}

# -----------------------------
# FST computation
# -----------------------------
compute_fst() {
    vcf=$1
    file1=$2
    file2=$3
    label=$4
    rep=$5

    prefix=${OUTDIR}/${label}_rep${rep}

    vcftools \
        --vcf ${vcf} \
        --weir-fst-pop ${file1} \
        --weir-fst-pop ${file2} \
        --fst-window-size ${WINDOW_SMALL} \
        --out ${prefix}_100kb > /dev/null 2>&1

    vcftools \
        --vcf ${vcf} \
        --weir-fst-pop ${file1} \
        --weir-fst-pop ${file2} \
        --fst-window-size ${WINDOW_LARGE} \
        --out ${prefix}_chrom > /dev/null 2>&1
}

# -----------------------------
# Pairwise avec nmin global
# -----------------------------
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

        # seed différent par réplication pour diversité mais reproductible
        seed_rep=$((SEED + rep))

        subsample_file ${pop1} ${nmin} ${seed_rep} ${sub1}
        subsample_file ${pop2} ${nmin} ${seed_rep} ${sub2}

        compute_fst ${vcf} ${sub1} ${sub2} ${label} ${rep}

        rm ${sub1} ${sub2}
    done
}

# -----------------------------
# MAIN LOOP
# -----------------------------

datasets=("3d" "litt")
thresholds=("pur" "majo")
Ks=("2" "3" "4" "7")

for dataset in "${datasets[@]}"
do
    if [ "$dataset" == "3d" ]; then
        VCF=${VCF_3D}
    else
        VCF=${VCF_LITT}
    fi

    for K in "${Ks[@]}"
    do
        for th in "${thresholds[@]}"
        do
            files=(${dataset}_K${K}_cluster*_${th}.txt)

            if [ ${#files[@]} -lt 2 ]; then
                continue
            fi

            echo "Processing ${dataset} | K${K} | ${th}"

            # -----------------------------
            # Calcul nmin GLOBAL
            # -----------------------------
            nmin_global=999999

            for f in "${files[@]}"
            do
                n=$(wc -l < ${f})
                if [ $n -lt $nmin_global ]; then
                    nmin_global=$n
                fi
            done

            echo "Global nmin for ${dataset} K${K} ${th} = ${nmin_global}"

            # -----------------------------
            # Comparaisons pairwise
            # -----------------------------
            for ((i=0; i<${#files[@]}; i++))
            do
                for ((j=i+1; j<${#files[@]}; j++))
                do
                    f1=${files[$i]}
                    f2=${files[$j]}

                    name1=$(basename ${f1} .txt)
                    name2=$(basename ${f2} .txt)

                    label="${name1}_VS_${name2}"

                    run_pairwise ${VCF} ${f1} ${f2} ${label} ${nmin_global}
                done
            done
        done
    done
done
