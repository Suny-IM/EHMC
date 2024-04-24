use strict;
use warnings;
use Data::Dumper;

#check if all cluster from mix

my $file1=$ARGV[0];
#combine_S_C.csv

open(FILE,"$file1") or die "can't open $file1!";

<FILE>;

my %hash;
my %hash_2;

while(my $line=<FILE>){
    chomp($line);
    my @data=split/,/,$line;
    $hash{$data[1]}++;
    if($data[0]=~/Cluster/){
        $hash_2{$data[1]}++;
    }
}
close(FILE);

for my $key (sort keys%hash){
    if($hash_2{$key}){
        if($hash{$key} eq $hash_2{$key}){
            print $key."\n";
        }
    }
}
