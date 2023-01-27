#!/usr/bin/perl
# Eric Morrison
# 01/25/2023
# Usage: B2GOslim2long.pl [input Blast2GO slim]
# This script takes a list of genes with slim terms computed by blast2GO and transforms it to "long" form, with each row representing a gene-slim term association instead of each row representing a gene with slim terms comma-separated in a single column.

use strict;
use warnings;

sub convert_line_feeds{
    my $in = $_[0];
    $in =~ s/\r|\r\n/\n/g;
    return($in);
}
                       
sub hash_GOs{
    my $inRef = $_[0];
    my @in = @$inRef;
    
    my %go;
    foreach my $gene (@in){
        
        my @gene = split("\t", $gene);
        #assign seqName and description
        my $seqName = $gene[2];
        my $descrip;
        if($gene[3] eq "" || $gene[3] eq "---NA---"){#some genes do not have description
            $descrip = "NA";
        }else{
            $descrip = $gene[3];
        }
        
        #For GO IDs and GO names need to check if defined
        if($gene[9] eq ""){
            %{ $go{$seqName} } = (
                "descrip" => $descrip,
                "aspect" => [ "NA" ],#assign as length 1 array for consistency
                "id" => [ "NA" ],
                "name" => [ "NA" ]
            );
            next;
        }
        
        #split GO IDs
        $gene[9] =~ s/F:|P:|C://g;
        my @goID = split("; ", $gene[9]);
        #split GO names
        my @goNames = split("; ", $gene[10]);
        
        #new arrays to separate aspect and names
        my @aspect;
        my @names;
        foreach my $goName (@goNames){
            $goName =~ /(F|C|P):(\w.+)/;
            push(@aspect, $1);
            push(@names, $2);
        }
        $go{$seqName}{"descrip"} = $descrip;
        $go{$seqName}{"aspect"} = [ @aspect ];
        $go{$seqName}{"id"} = [ @goID ];
        $go{$seqName}{"name"} = [ @names ];
    }
    return(\%go);
}

sub print_gos{
    my $goRef = $_[0];
    my %go = %$goRef;
    
    foreach my $gene (keys %go){
        for(my $i = 0; $i < scalar(@{ $go{$gene}{"aspect"} }); $i++){
            #print $i, "\n";
            print $gene, "\t", $go{$gene}{"descrip"}, "\t",
            ${ $go{$gene}{"aspect"} }[$i], "\t",
            ${ $go{$gene}{"id"} }[$i], "\t",
            ${ $go{$gene}{"name"} }[$i], "\n";
        }
    }
}

#MAIN
{
    my $in = $ARGV[0];
    open(IN, "$in") || die "Can't open input\n";
    chomp(my @in = <IN>);
    shift(@in); #remove header row
    print "SeqName\tDescription\tGOaspect\tGOID\tGOname\n";
    
    my $goRef = hash_GOs(\@in);
    print_gos($goRef);
}
