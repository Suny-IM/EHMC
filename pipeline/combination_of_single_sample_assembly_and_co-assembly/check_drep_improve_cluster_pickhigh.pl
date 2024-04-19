use strict;
use warnings;
use Data::Dumper;

#check if all cluster from mix

my $file1=$ARGV[0];
#combine_S_C_387cluster.csv


open(FILE,"$file1") or die "can't open $file1!";


my %hash_s;

while(my $line=<FILE>){
    chomp($line);
    my @data=split/,/,$line;
    if($hash_s{$data[1]}){
        if($data[6] > $hash_s{$data[1]}){
            $hash_s{$data[1]}=$data[6];
        }
    }else{
        $hash_s{$data[1]}=$data[6];
    }
}
close(FILE);

for my $key (sort keys%hash_s){
    print $key."\t".$hash_s{$key}."\n";
}
