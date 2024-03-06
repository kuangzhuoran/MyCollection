#usage perl <fasta> <filter_file>
#过滤多序列(CDS)对齐中的终止密码子
use warnings;
open($I,"<$ARGV[0]")  or die "No input fasta file!!\n"."Usage :perl filter.pl <fasta> <filter_file>; ";
open($O,">$ARGV[1]")  or die "Please assign output_file!!\n"."Usage :perl filter.pl <fasta> <filter_file>; ";
@array1 = <$I>;
@count = undef;
##捕获位点
for($i=0;$i<=$#array1;$i++)
	{
	if($array1[$i] =~ /\>/)
		{
		}	
	else	
	{	
		my @a = $array1[$i] =~ /(.{3})/g;
		for($n=0;$n<=$#a;$n++)
		{	
			my $c = $a[$n];
			chomp $c;
			if ($c eq "TGA" || $c eq "TAA" || $c eq "TAG" )
			{
			push @count,"$n";
			}
		}
		
	}
	}
##过滤密码子
for($i=0;$i<=$#array1;$i++)
	{
	if($array1[$i] =~ /\>/)
		{
		print {$O} "\n$array1[$i]";
		}	
	else	
	{	
		my @a = $array1[$i] =~ /(.{3})/g;
		for ($e=0;$e<=$#a;$e++)
		{	
			for ($g=1;$g<=$#count;$g++)
			{
			$a[$count[$g]] = "";
			}
			$c = $a[$e];
                        print {$O} $c;
		}
		}
		
	}
system("sed -i '1d' $ARGV[1]");
close $I;
close $O;
