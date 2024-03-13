#!/usr/bin/env perl 
use strict;

my $a = scalar(@ARGV);
print "$a---";
die "perl $0 <ncbi.gff3>   <out.longestID>  <out.ncbi_longest.gff3>\n" unless (@ARGV == 3);

my $gff = shift;
my $otab = shift;
my $ogff = shift;

## read GFF and get ID info
my %mRNA2cds;
my %gene2mRNA;
my %cdslen;
open GFF, $gff or die  $!;
while(<GFF>){
	chomp;
	next if /^#/;
	my @t = split /\t/;
	next unless  ($t[8] =~ /ID=/ &&  $t[8] =~ /Parent=/);

	$t[8] =~ /ID=([^;]+)/;
	my $id1 = $1;
	$t[8] =~ /Parent=([^;]+)/;
	my $id2 = $1;
	
	if($t[2] eq "mRNA"){
		$gene2mRNA{$id2}{$id1} = 1;
	}
	if($t[2] eq "CDS"){
		$mRNA2cds{$id2} = $id1;
		## get cds length;
		$cdslen{$id1} += $t[4] - $t[3] + 1;
	}
}
close GFF;


## get longest cds and mRNA ID;
open O1, ">$otab" or die  $!; 
my %lmRNA;
foreach my $geneID (keys %gene2mRNA) {
	my @mRNA = keys %{$gene2mRNA{$geneID}};
	my $len = 0;
	my $lmRNA;
	foreach my $mRNA (@mRNA){
		my $cds = $mRNA2cds{$mRNA};
		next if (!exists $cdslen{$cds});
		my $lcds = $cdslen{$cds};
		if($lcds > $len){
			$len = $lcds;
			$lmRNA = $mRNA;
		}
	}
	next if $len == 0;
	print O1 "$geneID\t$lmRNA\t$mRNA2cds{$lmRNA}\n";
	$lmRNA{$lmRNA} = $mRNA2cds{$lmRNA};
}
close O1;

## filter gff to keep longest cds and mRNA;


open O2, ">$ogff" or die  $!; 
open GFF, $gff or die  $!;
while(my $line = <GFF>){
	chomp $line ;
	next if $line =~ /^#/;
	my @t = split /\t/, $line;
	next unless  ( @t == 9);

	my $id1 = "NA";
	my $id2 = "NA";
	if($t[8] =~ /ID=([^;]+)/){
		$id1 = $1;
	}

	if($t[8] =~ /Parent=([^;]+)/){
		$id2 = $1;
	}
	
	my $out = $line ;
	my $flag = 0;
	if($t[2] eq "gene"){
		next unless exists $gene2mRNA{$id1}; 
		$flag = 1;
	}

	if($t[2] eq "mRNA"){
		next unless exists $lmRNA{$id1};
		$flag = 1;
	}

	if($t[2] =~ /CDS/i || $t[2] =~ /exon/i || $t[2] =~ /UTR/i ){
		next unless exists $lmRNA{$id2};
		$flag = 1;
	}
	if($flag == 1){	
		print O2 "$out\n";
	}
}
close GFF;
close O2;


