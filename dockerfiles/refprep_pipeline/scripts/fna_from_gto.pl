#!/usr/bin/env perl

use strict;
use JSON qw( decode_json );

my $gto_input = $ARGV[0];
my $fna_output = $ARGV[1];

open GTO, $gto_input;
read(GTO, my $gto_input, (stat(GTO))[7]);
close GTO;
print "$file start processing...\n";
my $decoded = decode_json($gto);
my $seedid = $decoded->{'source_id'};
print "$seedid start processing...\n";

my %contigs = (); #collect contigs DNA

for( my $i = 0; $i < scalar(@{$decoded->{'contigs'}}); $i++ ){
    $contigs{$decoded->{'contigs'}[$i]{'id'}} = $decoded->{'contigs'}[$i]{'dna'};
}

open FNA, ">$fna_output";
    
foreach my $contig ( keys %contigs ){
    print FNA ">$contig\n$contigs{$contig}\n\n";
}

close FNA;
print "$seedid finished processing...\n";
