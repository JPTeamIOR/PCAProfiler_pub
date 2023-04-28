#!/bin/bash

outdir=$(realpath ./export/)
cd ./output

mkdir -p ${outdir}


for name in ./* ; do
    bname=$(basename $name)
    quantSf="${name}/star_salmon/${bname}/quant.sf"
    if [ -f ${outdir}/${bname}/quant.sf ]; then
        gzip ${outdir}/${bname}/quant.sf
    else
        if [ -f ${outdir}/${bname}/quant.sf.gz ] ; then
            echo "${bname} skipped"
        else
            echo $bname 
            mkdir -p ${outdir}/${bname}/
            if [ -f ${quantSf} ]; then
                cp ${quantSf} ${outdir}/${bname}/
                gzip ${outdir}/${bname}/quant.sf
            else
                echo "Missing ${quantSf}" 
            fi
        fi
    fi
done

cd ..
