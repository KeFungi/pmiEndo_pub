#!/bin/bash

busco_tsv=$1
outpath=$2

if [[ ! -d $outpath ]]
    then mkdir $outpath
    else rm ${outpath}/*
    fi

spcodel=( $(head -1 $busco_tsv | cut -d$'\t' -f 2- ) )
buscoidl=( $(awk '{print $1}' $busco_tsv | sed '1d' ) )

for buscoid in ${buscoidl[@]}; do
    echo 'start' $buscoid
    geneline=( $(grep ^${buscoid} $busco_tsv | cut -d$'\t' -f 2-) )

    if [[ ${#spcodel[@]} != ${#geneline[@]} ]]; then
        echo 'table has different column number in rows'
        echo 'around' ${geneline}
        exit 1
    fi
    
    ncol=$(( ${#spcodel[@]} - 1 ))
    for i in $(seq -s' ' 0 $ncol); do
        spcode=${spcodel[$i]}
        gene=${geneline[$i]}
        echo '>'${spcode} >> ${outpath}/${buscoid}.fasta
        sed -n "/${gene}/,/>/p" ../CDS/${spcode}_CDS.fasta | sed '1d;$d' >> ${outpath}/${buscoid}.fasta
    done
done

