#Processing each line of input 
while (<>) {
@F=split/\s+/;

#Checking to see if each line starts with '>'
if (m/^\>/) {         
$F[1]=~s/^\[//;              # Remove '[' from the start of the second field
$F[3]=~s/\]$//;              # Remove ']' from the end of the fourth field
@G=split/\_/,$F[0];          # Split the first field on '_' into an array @G
$G[0]=~s/^\>//;              # Remove '>' from the start of the first element of @G
$len=abs($F[3]-$F[1]);       # Calculate the absolute length difference between the third and first fields

# Print the information of the sequence in forward orientation if the third field is greater than or equal to the first field
print "$G[0]\_$G[1]\_$G[2]\t$F[1]\t$F[3]\t+\t$len\t$G[3]\n" if $F[3] >= $F[1];

# Print the information of the sequence in reverse orientation if the third field is less than the first field
print "$G[0]\_$G[1]\_$G[2]\t$F[3]\t$F[1]\t-\t$len\t$G[3]\n" if $F[3] < $F[1];
}
}
