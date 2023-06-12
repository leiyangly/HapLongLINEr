while (<>) {
@F=split/\s+/;

print if ($F[3] eq "L1HS" or $F[3] eq "L1PA2" or $F[3] eq "L1PA3") and $F[2]-$F[1] >= 5000

}
