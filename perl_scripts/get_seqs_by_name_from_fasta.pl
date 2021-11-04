#!/usr/bin/perl
#Eric Morrison
#11032021
#Usage: get_seqs_by_name_from_fasta.pl fasta seq_header

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
    my($fasRef, $seqNam) = @_;
    my %fas = %$fasRef;
    print ">".$seqNam, "\n", $fas{$seqNam}, "\n";
}

#MAIN
{
    my $fasHashRef = hash_fasta();
    get_seqs($fasHashRef, $ARGV[1]);
}
