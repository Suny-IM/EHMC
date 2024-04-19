use strict;
use warnings;
use Data::Dumper;

#check if all cluster from mix

my $file1=$ARGV[0];
#genomeInformation.csv


open(FILE,"$file1") or die "can't open $file1!";

<FILE>;

my %hash;

while(my $line=<FILE>){
    chomp($line);
    my @data=split/,/,$line;
    $hash{$data[0]}{'completeness'}=$data[1];
    $hash{$data[0]}{'contamination'}=$data[2];
    $hash{$data[0]}{'length'}=$data[3];
    $hash{$data[0]}{'N50'}=$data[4];
}
close(FILE);

my $file2=$ARGV[1];

open(FILE,"$file2") or die "can't open $file2!";
#combine_S_C_*.csv


while(my $line=<FILE>){
    chomp($line);
    my @data=split/,/,$line;
    if($hash{$data[0]}){
        print $line."\t".$hash{$data[0]}{'completeness'}."\t".$hash{$data[0]}{'contamination'}."\t".$hash{$data[0]}{'length'}."\t".$hash{$data[0]}{'N50'}."\t";
        if(($hash{$data[0]}{'completeness'}>=95)&&($hash{$data[0]}{'contamination'}<=5)){
            print "A\n";
        }elsif(($hash{$data[0]}{'completeness'}>=90)&&($hash{$data[0]}{'completeness'}<95)&&($hash{$data[0]}{'contamination'}<=5)){
            print "B\n";
        }elsif(($hash{$data[0]}{'completeness'}<90)&&($hash{$data[0]}{'contamination'}<=5)){
            print "C\n";
        }else{
            print "D\n";
        }
    }
}

close(FILE);