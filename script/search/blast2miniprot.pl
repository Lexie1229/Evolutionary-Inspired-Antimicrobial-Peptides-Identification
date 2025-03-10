use strict;
use warnings;

while (my $line = <>) {
    chomp $line;
    my @fields = split /\t/, $line;
    
    if ($fields[2] eq "-") {
        ($fields[6], $fields[5]) = ($fields[5], $fields[6]);
    }
    
    print join("\t", @fields) . "\n";
}