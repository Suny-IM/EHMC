use strict;
use warnings;
use Data::Dumper;


my $file = $ARGV[0];


open(FILE,"$file") or die "can't open $file!";


while(my $line=<FILE>){
    chomp($line);
    my @data=split/\t/,$line;
    my @subdata=split/\|/,$data[0];
    if($subdata[-1]=~/g__/){
        if($data[1]>0){
            print $line."\n";
        }
        
    }
}

close(FILE);