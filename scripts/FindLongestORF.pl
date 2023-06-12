while (<>) {
@F=split/\s+/;

$name=$F[0];
$name=~s/\_[0-9]+$//;
$number=$F[0];
$number=~s/^.+\_//;
if ($F[3] >= $len{$name}{$F[1]}) {
$len{$name}{$F[1]}=$F[3];
$info{$name}{$F[1]}=$_;
}
$index{$name}=1;
END {
foreach $key (keys %index) {
@A=split/\t/,$info{$key}{L1rpORF1p};
@B=split/\t/,$info{$key}{L1rpORF2p};
@C=split/\_/,$A[0];
@D=split/\_/,$B[0];
$orf1len=$C[2]-$C[1];
$orf2len=$D[2]-$D[1];
#print "$orf1len\t$orf2len";
print "$info{$key}{L1rpORF1p}\t$info{$key}{L1rpORF2p}" if exists $info{$key}{L1rpORF1p} and exists $info{$key}{L1rpORF2p};
}}
}
