use strict;
use warnings;
use Data::Dumper;

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

my %contigs;

open(FILE,"$file1") or die "can't open!";

while(my $line=<FILE>){
    chomp($line);
    $contigs{$line}++;
}

close(FILE);

my $signal=0;

open(FILE2,"$file2") or die "can't open!";


while(my $line=<FILE2>){
    if($line=~/^>/){
        my @data=split/ /,$line;
        if($contigs{$data[0]}){
            print $line;
            $signal=1;
        }
        next;
    }
    if($signal==1){
       print $line;
       $signal=0;
    }
}

close(FILE2);