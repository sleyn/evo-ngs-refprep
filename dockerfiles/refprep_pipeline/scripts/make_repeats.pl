#!/usr/bin/env perl -w

use strict;

open IN, "$ARGV[0].aln";
open OUT, ">$ARGV[0].repeats";

my $s_count = 1;
my @coords = ();

while(<IN>){
	# skip headers
    if($s_count < 5){
        $s_count++;
        next;
    }
    
    chomp;
    my @hit = split /\t/, $_;
    
    # collect regions that have 85% or more identity
    if( $hit[6] > 85 ){
		if( $hit[11] eq $hit[12] ){
		    if( $hit[4] == $hit[7] ){
				next;
		    }
		    
		    if( $hit[0] == $hit[2] && $hit[1] == $hit[3] && $hit[6] == 100 ){
		    	next;
		    }
		}
		
		
		
        if( $hit[0] < $hit[1]){		#first region
            push @coords, [($hit[0],$hit[1],$hit[11])];
        }else{
            push @coords, [($hit[1],$hit[0],$hit[11])];
        }
        
        if( $hit[2] < $hit[3]){		#second region
            push @coords, [($hit[2],$hit[3],$hit[11])];
        }else{
            push @coords, [($hit[3],$hit[2],$hit[11])];
        }
    }
}

my %genome = ();
my %chromosomes = ();
my $min = 1;
my $max = 0;

for(my $i = 0; $i < (scalar @coords); $i++){
    if( $max < $coords[$i][1]){
        $max = $coords[$i][1];
    }
    
    if( $min > $coords[$i][0]){
        $min = $coords[$i][0];
    }
    
    # mark that there are some repeates in the contig
    $chromosomes{$coords[$i][2]} = 1;
    
    # write coverage by repeats to the contig
    for(my $j = $coords[$i][0]; $j <= $coords[$i][1]; $j++ ){
        $genome{$coords[$i][2]}{$j}++;
    }
}

my $left = 0;
my $right = 0;
my $check = 0;

open REPCOV, ">$ARGV[0].repcov";

foreach my $chrom ( keys %chromosomes ){
    for( my $i = $min; $i <= $max + 1; $i++ ){
        if( exists $genome{$chrom}{$i} ){
            print REPCOV "$chrom\t$i\t1\n";
            if( $check == 0 ){
                $left = $i;
                $check = 1;
            }
        }else{
            if( $check == 1){
                $right = $i - 1;
                print OUT "$chrom\t$left\t$right\n";
                $check = 0;
            }
        }
    }
}