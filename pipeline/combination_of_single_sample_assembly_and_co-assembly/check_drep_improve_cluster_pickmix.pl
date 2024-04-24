use strict;
use warnings;
use Data::Dumper;

#check if all cluster from mix

my $file1=$ARGV[0];
#cluster_*_single_score.txt


open(FILE,"$file1") or die "can't open $file1!";


my %hash_s;

while(my $line=<FILE>){
    chomp($line);
    my @data=split/\t/,$line;
    $hash_s{$data[0]}=$data[1];
}
close(FILE);



my $file2=$ARGV[1];
#combine_S_C_*cluster.csv

open(FILE2,"$file2") or die "can't open $file2!";

while(my $line=<FILE2>){
    chomp($line);
    my @data=split/,/,$line;
    if($data[6] > $hash_s{$data[1]}){
        print $line."\n";
    }
}
close(FILE2);


