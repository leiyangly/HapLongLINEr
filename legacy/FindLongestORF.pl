#Finding the longest Open Reading Frame (ORF) for each of the extracted L1 sequences 
while (<>) {
@F=split/\s+/;
chomp;

#Extracting the name of the L1
$name=$F[0];
$name=~s/\_[0-9]+$//;

#Extracting the number of the L1 sequence
$number=$F[0];
$number=~s/^.+\_//;

#Checking if the current ORF is longer than the previous longest ORF for the given L1
if ($F[3] >= $len{$name}{$F[1]}) {
$len{$name}{$F[1]}=$F[3];
$info{$name}{$F[1]}=$_;
}

#Setting the value of the current L1 sequence as 1 in the %index hash
$index{$name}=1;

END {
#Iterate over each L1 sequence stored in the %index hash
foreach $key (keys %index) {
#Split the information of the longest ORF1 and ORF2 for the current L1 sequence
@A=split/\t/,$info{$key}{L1rpORF1p};
@B=split/\t/,$info{$key}{L1rpORF2p};
@C=split/\_/,$A[0];
@D=split/\_/,$B[0];

#Calculate the lengths of ORF1 and ORF2
$orf1len=$C[2]-$C[1];
$orf2len=$D[2]-$D[1];

#print "$orf1len\t$orf2len";
print "$info{$key}{L1rpORF1p}\t$info{$key}{L1rpORF2p}\n" if exists $info{$key}{L1rpORF1p} and exists $info{$key}{L1rpORF2p};
}}
}
