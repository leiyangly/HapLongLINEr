while (<>) {
@F=split/\s+/;
print if $F[8] eq 1 and $F[9] eq 338 and $F[23] eq 1 and $F[24] eq 1275
}
