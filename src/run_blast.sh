#!/bin/bash 

if [ "$#" -ne 4 ]; then
    echo "Illegal number of parameters"
    exit
fi


BLASTDB=$1

sudo docker run --rm \
    -v $BLASTDB:/blast/blastdb:ro \
    -v $2:/blast/queries:ro \
    -v $2:/blast/results:rw \
    ncbi/blast \
    blastn -query /blast/queries/$3 -db $4 \
    -max_target_seqs 3 \
    -out /blast/results/$3_$4.out

