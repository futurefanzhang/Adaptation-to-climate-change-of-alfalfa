#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
  if (/^#/) {
    print;
    next;
  }
  chomp;
  my @a = split (/\t/, $_);
  my $ref = $a[3];
  my $alt = $a[4];
  for (my $i=9; $i <= $#a; $i++) {
    if ($a[$i] eq "0/0") {
      $a[$i] = "0";
    } elsif ($a[$i] eq "1/1") {
      $a[$i] = "2";
    } elsif ($a[$i] eq "0/1") {
      $a[$i] = "1";
        } elsif ($a[$i] eq "1/0") {
      $a[$i] = "1";
    } else {
      $a[$i] = "NA";
    }
  }
  print join "\t", @a;
  print "\n";
}
