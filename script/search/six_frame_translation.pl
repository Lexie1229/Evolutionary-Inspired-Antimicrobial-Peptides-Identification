use strict;
use warnings;
use autodie;

# six-frame translation
# perl six_frame_translation.pl in.fa out.fa

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

# Read fasta file: names and sequences
my ( $dnaseq, $names ) = read_fasta($infile);

# Six-frame translation
my $six_orf_seq = {};
foreach my $name ( @{$names} ) {

    my @orfs;
    push @orfs, ['+1', $dnaseq->{$name}];
    push @orfs, ['+2', substr $dnaseq->{$name}, 1]; 
    push @orfs, ['+3', substr $dnaseq->{$name}, 2];
    my $rc = reverse $dnaseq->{$name}; # Reverse
    $rc =~ tr/ATCGatcg/TAGCtagc/; # Complementary
    push @orfs, ['-1', $rc];
    push @orfs, ['-2', (substr $rc, 1)];
    push @orfs, ['-3', (substr $rc, 2)];

    my @pros;
    for my $i ( 0 .. $#orfs ) {
        my $count = $orfs[$i][0];
        my $orf = $orfs[$i][1];

        for ( my $j = 0; $j < ( length($orf) - 2 ); $j += 3 ) {
            my $codon = substr( $orf, $j, 3 );
            $pros[$i][0] = $count;
            $pros[$i][1] .= codon2aa($codon);
        }
    }

    $six_orf_seq -> {$name} = \@pros;

}

# Write to fasta file or stdout
my $fh_out;
if ( $outfile ) {
    open $fh_out, ">$outfile";
} else {
    $fh_out = \*STDOUT; 
}
for my $name (@{$names}) {
    for my $l (@{$six_orf_seq->{$name}}) {
        my ($c, $s) = @$l;
        $name =~ s/\s+//g;
        $c =~ s/\s+//g;
        print $fh_out ">$name" . ":" . "$c" . "\n";
        print $fh_out "$s\n";
    }
}
close $fh_out if $outfile;

exit;

# Read fasta file: names and sequences
sub read_fasta {
    my $fasta_file = shift;
    my ( $dnaseq, $names ) = ( {}, [] );
    my $cur_name;

    if ($fasta_file){
        open my $fh_in, '<', $fasta_file;
        while ( my $line = <$fh_in> ) {
            chomp $line;
            if ( $line =~ /^>/ ) {
                # Extract sequence name
                $cur_name = substr $line, 1;
                push @{$names}, $cur_name;
            }
            elsif ( $line =~ /^#/ ) {
                # Skip comment line
            }
            elsif ( $line =~ /^\w+/ ) {
                $line = uc $line; 
                next unless $cur_name;
                $dnaseq->{$cur_name} .= $line;
            }
        }
        close $fh_in;
    } else {
        # Read from stdin
        while ( my $line = <> ) {
            chomp $line;
            if ( $line =~ /^>/ ) {
                # Extract sequence name
                $cur_name = substr $line, 1;
                push @{$names}, $cur_name;
            }
            elsif ( $line =~ /^#/ ) {
                # Skip comment line
            }
            elsif ( $line =~ /^\w+/ ) {
                $line = uc $line;
                next unless $cur_name;
                $dnaseq->{$cur_name} .= $line;
            }
        }
    }

    return ( $dnaseq, $names );
}

# codon2aa: Codon comparison table
sub codon2aa {
    my ($codon) = @_;
    $codon = uc $codon;
    my (%genetic_code) = (
        'TCA' => 'S',    # Serine
        'TCC' => 'S',    # Serine
        'TCG' => 'S',    # Serine
        'TCT' => 'S',    # Serine
        'TTC' => 'F',    # Phenylalanine
        'TTT' => 'F',    # Phenylalanine
        'TTA' => 'L',    # Leucine
        'TTG' => 'L',    # Leucine
        'TAC' => 'Y',    # Tyrosine
        'TAT' => 'Y',    # Tyrosine
        'TAA' => '_',    # Stop
        'TAG' => '_',    # Stop
        'TGC' => 'C',    # Cysteine
        'TGT' => 'C',    # Cysteine
        'TGA' => '_',    # Stop
        'TGG' => 'W',    # Tryptophan
        'CTA' => 'L',    # Leucine
        'CTC' => 'L',    # Leucine
        'CTG' => 'L',    # Leucine
        'CTT' => 'L',    # Leucine
        'CCA' => 'P',    # Proline
        'CCC' => 'P',    # Proline
        'CCG' => 'P',    # Proline
        'CCT' => 'P',    # Proline
        'CAC' => 'H',    # Histidine
        'CAT' => 'H',    # Histidine
        'CAA' => 'Q',    # Glutamine
        'CAG' => 'Q',    # Glutamine
        'CGA' => 'R',    # Arginine
        'CGC' => 'R',    # Arginine
        'CGG' => 'R',    # Arginine
        'CGT' => 'R',    # Arginine
        'ATA' => 'I',    # Isoleucine
        'ATC' => 'I',    # Isoleucine
        'ATT' => 'I',    # Isoleucine
        'ATG' => 'M',    # Methionine
        'ACA' => 'T',    # Threonine
        'ACC' => 'T',    # Threonine
        'ACG' => 'T',    # Threonine
        'ACT' => 'T',    # Threonine
        'AAC' => 'N',    # Asparagine
        'AAT' => 'N',    # Asparagine
        'AAA' => 'K',    # Lysine
        'AAG' => 'K',    # Lysine
        'AGC' => 'S',    # Serine
        'AGT' => 'S',    # Serine
        'AGA' => 'R',    # Arginine
        'AGG' => 'R',    # Arginine
        'GTA' => 'V',    # Valine
        'GTC' => 'V',    # Valine
        'GTG' => 'V',    # Valine
        'GTT' => 'V',    # Valine
        'GCA' => 'A',    # Alanine
        'GCC' => 'A',    # Alanine
        'GCG' => 'A',    # Alanine
        'GCT' => 'A',    # Alanine
        'GAC' => 'D',    # Aspartic Acid
        'GAT' => 'D',    # Aspartic Acid
        'GAA' => 'E',    # Glutamic Acid
        'GAG' => 'E',    # Glutamic Acid
        'GGA' => 'G',    # Glycine
        'GGC' => 'G',    # Glycine
        'GGG' => 'G',    # Glycine
        'GGT' => 'G',    # Glycine
    );

    if ( exists $genetic_code{$codon} ) {
        return $genetic_code{$codon};
    }
    elsif ( $codon =~ /N/ ) {
        $genetic_code{$codon} = '*';  # If the codon contains N, * is returned
    }
}