#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# 定义变量
my ($maf_file, $output_file, $ref_species, $query_input, $tmp_file, $help);

# 解析参数
GetOptions(
    "i=s"     => \$maf_file,
    "o=s"     => \$output_file,
    "ref=s"   => \$ref_species,
    "query=s" => \$query_input,
    "tmp=s"   => \$tmp_file,
    "h|help"  => \$help,
) or die "Error in command line arguments\n";

# 如果什么参数都没给，或者显式要求帮助
if (!$maf_file && !$output_file && !$ref_species && !$query_input && !$tmp_file || $help) {
    print_help();
    exit;
}

# 检查必需参数
die "Missing -i input file\n" unless $maf_file;
die "Missing -o output file\n" unless $output_file;
die "Missing -ref reference species\n" unless $ref_species;
die "Missing -query query species\n" unless $query_input;

# 设定默认的tmp输出
$tmp_file ||= "no_ref.maf";

# 处理query输入（逗号分割 或 文件）
my @species;
if (-e $query_input) {
    open my $fh, "<", $query_input or die "Cannot open $query_input\n";
    while (<$fh>) {
        chomp;
        next if /^$/;
        push @species, $_;
    }
    close $fh;
} else {
    @species = split(/,/, $query_input);
}

# 打开输出文件
open my $O, ">", $output_file or die "Cannot open $output_file\n";
open my $E, ">", $tmp_file or die "Cannot open $tmp_file\n";

# 打印表头
my @head = ("chr", "pos", @species);
print $O join("\t", @head), "\n";

# 读取maf文件
open my $I, "-|", "cat $maf_file" or die "Cannot open $maf_file\n";
my $content = "";
while (<$I>) {
    next if(/^#/);
    if (/^a\s+score=/) {
        $content = $_;
        while (<$I>) {
            if (/^s\s+/) {
                $content .= $_;
            } else {
                last;
            }
        }
        analysis($content, $O, $E, $ref_species, \@species);
    }
}
close $I;
close $O;
close $E;

# 子程序：解析一段alignment block
sub analysis {
    my ($content, $O, $E, $ref, $species_ref) = @_;
    my @line = split(/\n/, $content);
    shift @line; # 删除a行
    my %pos;
    my $isRefFound = 0;
    my $ref_chr = "NA";

    foreach my $line (@line) {
        my @a = split(/\s+/, $line);
        my ($s, $chr, $start, $alignment_length, $strand, $sequence_length, $sequence) = @a;
        $chr =~ /^([^\.]+)\.(.*)/;
        my $species = $1;
        $chr = $2;
        if ($species eq $ref) {
            $ref_chr = $chr;
            $isRefFound = 1;
            my @base = split(//, $sequence);
            if ($strand eq "+") {
                my $pos = $start;
                for (my $i = 0; $i < @base; $i++) {
                    if ($base[$i] ne "-") {
                        $pos++;
                        $pos{$i} = $pos;
                    }
                }
            } elsif ($strand eq "-") {
                my $pos = $start;
                for (my $i = 0; $i < @base; $i++) {
                    if ($base[$i] ne "-") {
                        $pos++;
                        my $real_pos = $sequence_length - 1 - ($pos - 1) + 1;
                        $pos{$i} = $real_pos;
                    }
                }
            }
        }
    }

    if ($isRefFound == 0) {
        print $E $content;
        return;
    }

    my %data;
    foreach my $line (@line) {
        my @a = split(/\s+/, $line);
        my ($s, $chr, $start, $alignment_length, $strand, $sequence_length, $sequence) = @a;
        $chr =~ /^([^\.]+)\.(.*)/;
        my $species = $1;
        $chr = $2;
        my @base = split(//, $sequence);
        $sequence =~ tr/ATCGRYMK/TAGCYRKM/ if ($strand eq "-");
        for (my $i = 0; $i < @base; $i++) {
            next unless exists $pos{$i};
            my $pos = $pos{$i};
            $data{$pos}{$species} = $base[$i];
        }
    }

    foreach my $pos (sort { $a <=> $b } keys %data) {
        my @output_line = ($ref_chr, $pos);
        foreach my $species (@$species_ref) {
            my $base = "-";
            $base = $data{$pos}{$species} if exists $data{$pos}{$species};
            push @output_line, $base;
        }
        print $O join("\t", @output_line), "\n";
    }
}

# 打印帮助信息
sub print_help {
    print <<EOF;

Usage:
    perl 01.convertMaf2List.pl -i <input.maf> -o <output.lst> -ref <reference_species> -query <species1,species2,...|species.list> [-tmp <no_ref_output.maf>]

Parameters:
    -i       Input MAF file (required)
    -o       Output LST file (required)
    -ref     Reference species name (required)
    -query   Query species names (comma-separated) or a file with species list (required)
    -tmp     Output file for MAF blocks without reference species (optional, default: no_ref.maf)
    -h, --help  Show this help message

Examples:
    perl 01.convertMaf2List.pl -i 5zokors.maf -o 5zokors.maf.lst -ref Efontanierii -query Ebaileyi,Esmithi,Erufescens,Erothschildi
    perl 01.convertMaf2List.pl -i 5zokors.maf -o 5zokors.maf.lst -ref Efontanierii -query species.list

EOF
}

