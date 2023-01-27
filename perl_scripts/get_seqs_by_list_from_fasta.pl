#!/usr/bin/perl
#Eric Morrison
#11032021
#Usage: get_seqs_by_name_from_fasta.pl fasta seq_header_file
#Note the script removes everything in the fasta header after the first space

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
        $header =~ /(.+?)\s/;
        $header = $1;
        #print $header, "\n";
        $fas{$header} = join("", @seq);
    }
    return(\%fas);
}

sub get_seqs{
    my($fasRef, $idsHash) = @_;
    my %fas = %$fasRef;
    my @ids = @$idsHash;
    
    foreach my $seqNam (@ids){
        print ">".$seqNam, "\n", $fas{$seqNam}, "\n";
    }
}

#MAIN
{
    my $fasHashRef = hash_fasta();
    my $idList = $ARGV[1];
    open(IDS, "$idList") || die "Can't open IDs\n";
    chomp(my @ids = <IDS>);
    get_seqs($fasHashRef, \@ids);
}
