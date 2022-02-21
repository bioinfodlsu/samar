#! /bin/bash

test $# = 1 || {
    cat <<EOF
Usage: $0 x.fastq
 or:   $0 x.fastq.gz

Read a fasts file, and output fasta file.
Drop 3rd and 4th lines of fastq.
Keep just the first word of header lines.
Assumes 1 fasta per 2 lines, i.e. no line wrapping.
EOF
    exit
}

gzip -cdf "$1" |
   awk '(NR%4 == 1)||(NR%4==2) ' |
   sed -e 's/^@ */\>/' -e 's/ .*//'