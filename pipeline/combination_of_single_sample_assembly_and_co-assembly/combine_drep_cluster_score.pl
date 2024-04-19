use strict;
use warnings;
use Data::Dumper;

my $file1=$ARGV[0];
#Sdb.csv


open(FILE,"$file1") or die "can't open $file1!";

<FILE>;

my %hash;

while(my $line=<FILE>){
    chomp($line);
    my @data=split/,/,$line;
    $hash{$data[0]}=$data[1];
}
close(FILE);


my $file2=$ARGV[1];
#Cdb.csv


open(FILE2,"$file2") or die "can't open $file2!";

my $header=<FILE2>;
chomp($header);
print $header.",score\n";

while(my $line=<FILE2>){
    chomp($line);
    my @data=split/,/,$line;
    if($hash{$data[0]}){
        print $line.",$hash{$data[0]}\n";
    }
}

close(FILE2);