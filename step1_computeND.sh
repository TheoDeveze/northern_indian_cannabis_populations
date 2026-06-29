#!/bin/bash

module load StdEnv/2023
module load bcftools
module load vcftools

VCF_3D="renamed_3D_forFST.vcf"
VCF_LITT="renamed_litt_forFST.vcf"
NBOOT=1000
OUTDIR="PI_site_rarefaction"

SEED=420

mkdir -p ${OUTDIR}


# subsampling
subsample_file() {

    infile=$1
    n=$2
    seed=$3
    outfile=$4

    shuf --random-source=<(yes ${seed}) ${infile} \
        | head -n ${n} > ${outfile}
}


# compute per site pi

compute_pi() {

    vcf=$1
    keepfile=$2
    label=$3
    rep=$4

    prefix=${OUTDIR}/${label}_rep${rep}

    vcftools \
        --vcf ${vcf} \
        --keep ${keepfile} \
        --site-pi \
        --out ${prefix} \
        > /dev/null 2>&1
}


# run rarefaction
run_rarefaction() {

    vcf=$1
    popfile=$2
    label=$3
    nmin=$4

    echo "Running ${label} with n=${nmin}"

    for ((rep=1; rep<=NBOOT; rep++))
    do

        sub="tmp_${label}_${rep}.txt"

        # seed
        seed_rep=$((SEED + rep))

        subsample_file ${popfile} ${nmin} ${seed_rep} ${sub}

        compute_pi ${vcf} ${sub} ${label} ${rep}

        rm ${sub}

    done
}

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

            if [ ${#files[@]} -lt 1 ]; then
                continue
            fi

            echo "================================================="
            echo "Processing ${dataset} | K${K} | ${th}"
            echo "================================================="

            # -------------------------------------------------
            # GLOBAL MIN SAMPLE SIZE
            # -------------------------------------------------

            nmin_global=999999

            for f in "${files[@]}"
            do

                n=$(wc -l < ${f})

                if [ $n -lt $nmin_global ]; then
                    nmin_global=$n
                fi

            done

            echo "Global nmin = ${nmin_global}"

            # -------------------------------------------------
            # RUN EACH CLUSTER
            # -------------------------------------------------

            for f in "${files[@]}"
            do

                name=$(basename ${f} .txt)

                run_rarefaction \
                    ${VCF} \
                    ${f} \
                    ${name} \
                    ${nmin_global}

            done

        done
    done
done
