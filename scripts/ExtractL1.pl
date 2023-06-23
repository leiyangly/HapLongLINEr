while (<>) {
@F=split/\s+/;

print if $F[3]=~m/^L1/ and $F[2]-$F[1] >= 5000

}
