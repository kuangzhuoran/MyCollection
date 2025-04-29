#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# 定义变量
my ($lst_file, $gff_file, $outdir, $help);

# 参数解析
GetOptions(
    "lst=s"   => \$lst_file,
    "gff=s"   => \$gff_file,
    "outdir=s" => \$outdir,
    "h|help"  => \$help,
) or die "Error in command line arguments\n";

# 如果没给参数或者要求帮助
if (!$lst_file && !$gff_file || $help) {
    print_help();
    exit;
}

# 默认输出目录
$outdir ||= "genes";

# 创建输出目录
`mkdir -p $outdir` unless -e $outdir;

# 读取GFF
my %gene;
my %pos;
open my $GFF, "<", $gff_file or die "Cannot open $gff_file\n";
while (<$GFF>) {
    chomp;
    my @a = split(/\s+/, $_);
    next unless ($a[2] eq "CDS");
    my ($chr, $start, $end, $strand) = ($a[0], $a[3], $a[4], $a[6]);
    $a[8] =~ /Parent=(\w.*)/;
    my $id = $1;
    for (my $i = $start; $i <= $end; $i++) {
        $pos{$chr}{$i} = $id;
        $gene{$id}{pos}{$i} = 0;
    }
    $gene{$id}{strand} = $strand;
}
close $GFF;
print STDERR "$gff_file loaded...\n";

# 读取.lst文件
my %seq;
open my $LST, "<", $lst_file or die "Cannot open $lst_file\n";
my $head = <$LST>;
my @head = split(/\s+/, $head);
my $control = 0;
while (<$LST>) {
    my @a = split(/\s+/, $_);
    my ($chr, $pos) = ($a[0], $a[1]);
    next unless exists $pos{$chr}{$pos};
    print STDERR "[ $control ] sites loaded...\r" if ($control % 1000 == 0);
    my $id = $pos{$chr}{$pos};
    for (my $i = 2; $i < @a; $i++) {
        my $species = $head[$i];
        $a[$i] = uc($a[$i]);
        $seq{$id}{$species}{$pos} = $a[$i];
    }
    $control++;
}
close $LST;
print STDERR "$lst_file loaded...\n";

# 输出每个基因
foreach my $id (sort keys %seq) {
    print STDERR "$id\r";
    open my $OUT, ">", "$outdir/$id.fa" or die "Cannot write to $outdir/$id.fa\n";
    foreach my $species (sort keys %{$seq{$id}}) {
        my $nucl;
        foreach my $pos (sort { $a <=> $b } keys %{$gene{$id}{pos}}) {
            my $base = "-";
            $base = $seq{$id}{$species}{$pos} if exists $seq{$id}{$species}{$pos};
            $nucl .= $base;
        }
        my $strand = $gene{$id}{strand};
        if ($strand eq "-") {
            $nucl =~ tr/ATCGRYMK/TAGCYRKM/;
            $nucl = reverse($nucl);
        }
        print $OUT ">$species\n$nucl\n";
    }
    close $OUT;
}
print STDERR "\nDone\n";

# 帮助信息
sub print_help {
    print <<EOF;

Usage:
    perl 02.lst2gene.pl -lst <input.maf.lst> -gff <reference.gff> [-outdir <output_directory>]

Parameters:
    -lst       Input MAF-LST file (required)
    -gff       GFF annotation file (required)
    -outdir    Output directory for gene fasta files (optional, default: genes/)
    -h, --help Show this help message

Examples:
    perl 02.lst2gene.pl -lst 5zokors.maf.lst -gff Efontanierii.sorted.rename.gff
    perl 02.lst2gene.pl -lst 5zokors.maf.lst -gff Efontanierii.sorted.rename.gff -outdir genes_output

EOF
}

