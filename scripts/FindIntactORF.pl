#Locating the intact open reading frames (ORFs) in each of the extracted L1 sequences by finding the longest alignment to ORF1 and ORF2 among all found ORFs using the Blastp file as input.
while (<>) {
@F=split/\s+/;
print if $F[8] eq 1 and $F[9] eq 338 and $F[23] eq 1 and $F[24] eq 1275
}
