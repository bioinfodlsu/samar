#! /bin/bash

test $# = 2 || {
    cat <<EOF
Usage: $0 x.fasta y.fasta
or:    $0 x.fasta.gz y.fasta.gz

Read 2 fasta files, and write them interleaved.  Keep just the first
word of header lines, and append "/1" and "/2" if they are otherwise
identical.  Assumes 1 fasta per 2 lines, i.e. no line wrapping.
EOF
    exit
}

fastaTab () {
    gzip -cdf "$@" |
        sed -e 's/^>  */>/' -e 's/ .*//' |
        paste - -
}

paste <(fastaTab "$1") <(fastaTab "$2") |
    awk '$1 == $3 {$1 = $1 "/1"; $3 = $3 "/2"} $1 = $1' OFS="\n"

