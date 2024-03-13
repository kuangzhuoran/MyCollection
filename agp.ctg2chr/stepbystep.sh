#参考https://hub.docker.com/r/informationsea/transanno
#参考https://github.com/informationsea/transanno

#最基本的使用是下面三行命令
minimap2 -cx asm5 --cs query.fa ref.fa > out.paf  #query.fa就是我们的目标，ref.fa是转换之前的vcf，我们希望把vcf的基因组位置信息转成query.fa
transanno minimap2-to-chain out.paf --output out.chain #当然不光可以转vcf，总的来说是基因组坐标的转换
transanno liftvcf --chain out.chain --output final.vcf --fail fail.vcf --new-assembly query.fa --original-assembly ref.fa --vcf old.vcf

#这个脚本处理3DDNA输出的结果，可以得到 *.agp，也就是contig级别和染色体级别基因组的位置对应关系
python juicer_assembly2agp_fa.py review.assembly ctg.fa Chr 32 0
#.agp文件转换vcf
grep -w -v 'N' Galili.Chr.agp | grep -v '#' > agp.1
python 1.py galili.ctg.length agp.1 > agp.2  #contig级别基因组的长度信息, 每一行是这样的:   ptg00001    101241
python 2.py Galili.Chr.length agp.2 > agp.3  #染色体级别基因组的长度信息, 这样制取bioawk -c fastx '{print $name, length($seq)}' galili.ctg.fa > galili.ctg.length
awk 'BEGIN{ FS="\t";OFS="\t" }{print $6,$11,$7-1,$8,$9,$1,$13,$2-1,$3,$8-$7+1,$8-$7+1,'60',"cg:Z:"$8-$7+1"M","cs:Z::"$8-$7+1}' agp.3 | sort > out.paf
rm agp.1 agp.2 agp.3
/data/01/user164/software/transanno-x86_64-unknown-linux-musl-v0.3.0/transanno minimap2chain out.paf --output out.chain
transanno liftvcf --chain out.chain --output final.vcf --fail fail.vcf --new-assembly query.fa --original-assembly ref.fa --vcf old.vcf

CrossMap.py gff carmeli.chain carmeli.ctg.gff | awk -F "->" '{print $2}' | sed 's/^[\t ]\+//' > carmeli.chr.tmp
bedtools sort -i  carmeli.chr.tmp > carmeli.chr.gff


#下面是实际运行的命令行
minimap2 -cx asm5 --cs /data/01/user164/workspace_for_whole_object/spalax_chrSV/01.Assembly/4.Mount/galili/Galili.Chr.fasta /data/01/user164/workspace_for_whole_object/spalax_chrSV/01.Assembly/4.Mount/galili/galili.ctg.fa > galili.paf
/data/01/user164/software/transanno-x86_64-unknown-linux-musl-v0.3.0/transanno minimap2chain galili.paf --output galili.chain
/data/01/user164/software/transanno-x86_64-unknown-linux-musl-v0.3.0/transanno liftvcf --chain galili.chain --output RefGalili.NoOut.onlygalili.chr.vcf --fail galili.fail.vcf --new-assembly /data/01/user164/workspace_for_whole_object/spalax_chrSV/01.Assembly/4.Mount/galili/Galili.Chr.fasta --original-assembly /data/01/user164/workspace_for_whole_object/spalax_chrSV/01.Assembly/4.Mount/galili/galili.ctg.fa --vcf /data/01/user164/workspace_for_whole_object/spalax_chrSV/00.tidy/SNP/RefGalili.NoOut.onlygalili.ctg.vcf
