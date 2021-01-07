#! /bin/bash

test $# = 2 || {
    cat <<EOF
Usage: $0 x.fastq y.fastq
 or:   $0 x.fastq.gz y.fastq.gz 

Read 2 fastq files, and write them interleaved
Drop the quality info and 3rd line of fastq. 
Keep just the first word of header lines, and append "/1" and "/2" if they are otherwise
identical.  Assumes 1 fasta per 2 lines, i.e. no line wrapping.
EOF
    exit
}

fastaTab () {
    gzip -cdf "$@" |
        sed -e 's/^@ */\>/' -e 's/ .*//' |
        paste - - - - | cut -f 1,2
}

paste <(fastaTab "$1") <(fastaTab "$2") |
    awk '$1 == $3 {$1 = $1 "/1"; $3 = $3 "/2"} $1 = $1' OFS="\n"

