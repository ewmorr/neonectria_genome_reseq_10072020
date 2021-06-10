#!/usr/bin/perl
#Eric Morrison
#06072021

use strict;
use warnings;

#Main
{
    my $in = $ARGV[0];
    open(IN, "$in") || die "Can't open file.\n";
    
    my $counter = 0;
    my $average = 0;
    my @values;
    
    while(my $line = <IN>){
        $line =~ /DP=(\d+);/;
        $counter++;
        $average += $1;
        push(@values, $1);
    }
    #$average = $average/$counter;
    my $median;
    my $mid = int @values/2;
    my @sorted_values = sort {$a <=> $b} @values;
    if (@values % 2) {
        $median = $sorted_values[ $mid ];
    } else {
        $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
    }

    print "Average depth is ",$average/$counter, ". Median depth is ",$median, ". Total sites is ", $counter, "\n"
}
