#!/usr/bin/perl
#Eric Morrison
#11032021
#Usage: get_segments_by_pos_from_fasta.pl fasta seq_header pos1 pos2

use strict;
use warnings;

sub hash_fasta{
    my $fasta = $ARGV[0];
    open(FAS, "$fasta") || die "can't open fasta\n";
    chomp(my @fasta = <FAS>);
    my $fas = join("&&&&&&&&&&&&", @fasta);
    my @fas = split(">", $fas);
    shift(@fas);
    my %fas;
    foreach my $seq (@fas){
        my @seq = split("&&&&&&&&&&&&", $seq);
        my $header = shift(@seq);
        $fas{$header} = join("", @seq);
    }
    return(\%fas);
}

sub get_seqs{
    my($fasRef, $seqNam, $pos1, $pos2) = @_;
    my %fas = %$fasRef;
    my $seq = substr($fas{$seqNam}, $pos1-1, $pos2-$pos1+1);
    print ">".$seqNam."_".$pos1."-".$pos2, "\n", $seq, "\n";
}

#MAIN
{
    my $fasHashRef = hash_fasta();
    get_seqs($fasHashRef, @ARGV[1..3]);
    
}
