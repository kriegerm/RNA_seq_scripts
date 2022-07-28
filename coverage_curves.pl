#!/usr/bin/perl
​
unless(@ARGV) { die "usage:\tcoverage_curves.pl input.sam  > output.txt\n\n"; }
​
$file = $ARGV[0];
​
open(IN,"<$file");
while($l =<IN>){
        if($l=~/^HWI-/) { ### MIGHT HAVE TO CHANGE TO MATCH READ NAME CONVENTION
                # NOTE THAT STRAND ASSIGNMENT IS REVERSED TO COMPENSATE FOR PROBLEM WITH ILLUMINA READS
                @cols=split(/\t/,$l);
                if (
                        $cols[1] == 0 || 
                        $cols[1] == 99 ||
                        $cols[1] == 147 ||
                        $cols[1] == 97 ||
                        $cols[1] == 145
                        ) { 
                                $strand = "R"; 
                        }
                elsif (
                        $cols[1] == 16 ||
                        $cols[1] == 83 ||
                        $cols[1] == 163 ||
                        $cols[1] == 81 ||
                        $cols[1] == 161
                        ) {
                                $strand = "F";
                        } 
                # These values indicate strand is ambiguous or poor quality; no strand assigned
                elsif (
                        $cols[1] == 113 ||
                        $cols[1] == 177 ||
                        $cols[1] == 65 ||
                        $cols[1] == 129 ||
                        $cols[1] == 4 ||
                        $cols[1] == 141 ||
                        $cols[1] == 77
                        ) {
                        } 
                else { print "ERROR:[$cols[1]]\n"; } 
​
                $cols[5]=~/(\d+)M/;
                $length = $1;
                for (1..$length) {
                        $COVERAGE{$strand}{$cols[3]}++;
                        $cols[3]++;
                }
        }elsif($l=~/\s+LN\:(\d+)\s/) {
                $gsize = $1;
                #print "genome size=[$gsize]\n";
        }
        
}
​
for $i (1..$gsize){
#for $i (1..500){
​
        # print "$i\t";## uncomment if you want the position listed
        if(exists($COVERAGE{F}{$i})){
                print "$COVERAGE{F}{$i}\t";
        } else {
                print "0\t";
        }
        if(exists($COVERAGE{R}{$i})){
                print "$COVERAGE{R}{$i}\n";
        } else {
                print "0\n";
        }       
}
