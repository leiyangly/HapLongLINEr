while (<>) {
@F=split/\s+/;

if (m/^\>/) {
$F[1]=~s/^\[//;
$F[3]=~s/\]$//;
@G=split/\_/,$F[0];
$G[0]=~s/^\>//;
$len=abs($F[3]-$F[1]);
print "$G[0]\_$G[1]\_$G[2]\t$F[1]\t$F[3]\t+\t$len\t$G[3]\n" if $F[3] >= $F[1];
print "$G[0]\_$G[1]\_$G[2]\t$F[3]\t$F[1]\t-\t$len\t$G[3]\n" if $F[3] < $F[1];
}
}
