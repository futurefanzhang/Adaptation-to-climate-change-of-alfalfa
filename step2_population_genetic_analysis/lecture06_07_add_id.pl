#!/usr/bin/perl -w
if(@ARGV != 2){
	print("perl add_id.pl in.vcf out.vcf\n");
	die();
}
open(IN, $ARGV[0]);
open(OUT, ">$ARGV[1]");
while (<IN>) {
	chomp;
	next if(/^\s+$/);
	if(/^\#/){
		print OUT "$_\n";
		next;
	}
	my ($chr, $pos, $id, $others) = split(/\t/, $_, 4);
	$id = join("__", $chr, $pos);
	print OUT "$chr\t$pos\t$id\t$others\n";
}
close(IN);
close(OUT);
