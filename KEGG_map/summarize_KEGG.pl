#!/usr/bin/perl
#summarize_KEGG.pl
#Eric Morrison
#01292021
#Usage: perl summarize_KEGG.pl [KEGG_map] [query KO list]
#KEGG_map should point to KEGG database, query list should be a list of KOs in the format "KO:Kxxxxx" one per line, the output is a table of mappings with columns in order of level A, B, C, D, with one entry per KO mapping (some KOs may have multiple mappings)

use strict;
use warnings;

sub hash_KEGG_db{
    my $kegg_db = $_[0];
    open(KEGG, "$kegg_db") || die "Can't open KEGG database.\n";
    chomp(my @kegg_db = <KEGG>);
    
    my $A_level;
    my $B_level;
    my $C_level;
    my $D_level;
    my %kegg_db;
    #DB hash has D level as top level key
    #each entry then stores levels in reverse order so that key count of second hash level is equal to number of KO mappings
    foreach my $line (@kegg_db){
        #assign level to vars
        if($line =~ /^A\d+/){
            $A_level = $line;
        }
        if($line =~ /^B\s+/){
            $line =~ s/\s+/ /g;
            $B_level = $line;
        }
        if($line =~ /^C\s+/){
            $line =~ s/\s+/ /g;
            $C_level = $line;
        }
        if($line =~ /^D\s+/){
            $line =~ s/^D\s+(K\d+)\s//g;
            $D_level = $1;
            #if D level assign hash (all other levels should be assigned)
            my $AtoC = $A_level."\t".$B_level."\t".$C_level;
            $kegg_db{$D_level}{$AtoC} = 1;
        }

    }
    return(\%kegg_db);
}

sub filter_hash{
    my($kos_ref, $kegg_db_ref) = @_;
    my @kos = @$kos_ref;
    my %kegg_db = %$kegg_db_ref;
    
    foreach my $query (@kos){
        $query =~ s/^KO://;
        foreach my $mapping (sort{$a cmp $b} keys %{ $kegg_db{$query} } ){
            print $mapping, "\t", $query, "\n";
        }
    }
}

sub strip_win_line_feeds{
    my $kos_ref = $_[0];
    my @kos = @$kos_ref;
    
    if(scalar(@kos) == 1){
        $kos[0] = s/\r|\r\n|\n/\n/g;
        my @new_kos = split("\n", $kos[0]);
        return(\@new_kos);
    }else{
        return(\@kos);
    }
}

#MAIN
{
    my $kegg_db = $ARGV[0];
    my $ko_query = $ARGV[1];
    open(KOS, "$ko_query") || die "Can't open KO query list.\n";
    chomp(my @kos = <KOS>);
   
    my $kos_ref = strip_win_line_feeds(\@kos);
    my $kegg_db_ref = hash_KEGG_db($kegg_db);
    filter_hash($kos_ref, $kegg_db_ref)
}






