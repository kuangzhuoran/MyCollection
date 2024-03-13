#perl calculate_base_content.pl 染色体.fasta
#
#!/usr/bin/perl

use strict;
use warnings;

# 获取FASTA文件路径
my $file = shift @ARGV;
die "未指定FASTA文件。请提供文件路径作为参数。\n" unless defined $file;
die "无法打开文件 '$file'：$!\n" unless -e $file;

# 读取FASTA文件
open(my $fh, '<', $file) or die "无法打开文件 '$file'：$!";
my $current_sequence = '';
while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^>/) {
        # 忽略注释行
        next;
    } else {
        $current_sequence .= $line;
    }
}
close $fh;

# 计算碱基含量
my $total_length = length($current_sequence);
my $a_count = ($current_sequence =~ tr/Aa//);
my $t_count = ($current_sequence =~ tr/Tt//);
my $c_count = ($current_sequence =~ tr/Cc//);
my $g_count = ($current_sequence =~ tr/Gg//);

# 输出结果
print "A含量: " . ($a_count / $total_length) * 100 . "%\n";
print "T含量: " . ($t_count / $total_length) * 100 . "%\n";
print "C含量: " . ($c_count / $total_length) * 100 . "%\n";
print "G含量: " . ($g_count / $total_length) * 100 . "%\n";
