#Working to extract the sequences from the full legnth L1s (sequences >= 5000bp) from the input .bed file.

while (<>) {
@F=split/\s+/;

print if $F[3]=~m/^L1/ and $F[2]-$F[1] >= 5000

}
