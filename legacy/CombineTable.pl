#Integrating the ORF finding and liftover information into a single file output.
while (<>) {
@F=split/\s+/;

#Combine lines into %PLUS hash based on conditions
if (scalar(@ARGV) eq 3) {$PLUS{$F[0]}=$_ if $F[3]-$F[2] >= 200 and (!$PLUS{$F[0]} or $PLUS{$F[0]}=~m/\_/g)}

#Store lines in the %MINUS hash if specific conditions are met
if (scalar(@ARGV) eq 2) {$MINUS{$F[0]}=$_ if $F[3]-$F[2] >= 200 and (!$MINUS{$F[0]} or $MINUS{$F[0]}=~m/\_/g)}

#Taking the lines of input and splitting it into the G-array and storing the new strings in the %INTACT hash 
if (scalar(@ARGV) eq 1) {$"="\_"; @G=split/\_/, $F[0]; $len=scalar(@G)-4; $name="@G[0..$len]";$g=$name."\_".$G[-3]."\_".$G[-2];$INTACT{$g}=$_;}

#Performing multiple calculations and generating final output
if (scalar(@ARGV) eq 0) {
$m=$F[1]-1999;
$p=$F[2]+2000;
$px=$F[2]+1;
$mx=$F[1]+1;
$mkey=$F[0]."\:".$m."\-".$F[1];
$pkey=$F[0]."\:".$px."\-".$p;
$ikey=$F[0]."\_".$mx."\_".$F[2];

#Split the relevant values into the specified arrays for further processing
@M=split/\t/,$MINUS{$mkey};
@P=split/\t/,$PLUS{$pkey};
@I=split/\t/,$INTACT{$ikey};

#Determining additional attributes for the final output 
$intact="present";
$intact="intact" if $INTACT{$ikey};
$chr="NA";
$chr=$M[5] if $M[5] eq $P[5];
$strand="NA";
$strand="\+" if $M[4] eq $P[4] and $M[4] eq $F[5];
$strand="\-" if $M[4] eq $P[4] and $M[4] ne $F[5];
$start="NA";
$start=$M[8] if $M[7] <= $P[8] and $M[5] eq $P[5];
$start=$P[8] if $P[7] <= $M[8] and $M[5] eq $P[5];
$end="NA";
$end=$P[7] if $M[7] <= $P[8] and $M[5] eq $P[5];
$end=$M[7] if $P[7] <= $M[8] and $M[5] eq $P[5];
if ($end < $start) {$swap=$end; $end=$start; $start=$swap}

#Generate the final output line using extracted data and print it as the script's output
print "$F[0]\_$F[1]\_$F[2]\_$F[5]\_$F[4]\_$F[3]\_$intact\t$chr\_$start\_$end\_$strand\n";
}
}
