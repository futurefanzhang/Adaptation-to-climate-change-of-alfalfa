#!/usr/bin/perl
use strict;
use warnings;

die "\nUSAGE:perl $0 raw.vcf  snp_name.txt  out.vcf\n" if(scalar @ARGV!=3);
my %SV;
open IN,$ARGV[0] or die $!;# raw vcf file
while(<IN>){
        if (/^#/) {
    next;
  }
        chomp;
        my ($chr, $pos, $id)=split /\t/,$_,3;
        $SV{$id}="$chr\t$pos\t$id";
}
close IN;

open IN,$ARGV[1] or die $!; #snp name file
open OUT,">$ARGV[2]" or die $!;
while(<IN>){
        chomp;
        my ($p_value,$indel,$fdr)=split '\t',$_,3;
        #print $SVid;
        if(exists $SV{$indel}){

                print OUT "$SV{$indel}\t$p_value\t$fdr\n";

        }
}
close IN;
close OUT;
