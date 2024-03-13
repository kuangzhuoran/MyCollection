nucmer --mum --mincluster 500 -t 30 ref.fasta  query.fasta -p ref_query
delta-filter -m -i 90 -l 100 ref_query.delta > ref_query.filter.delta
show-coords -c -r ref_query.filter.delta > ref_query.filter.coords

sed '1,5d' galiliv3_golaniv4.filter.Chr.coords | awk 'BEGIN{ OFS="\t" }{if($4<$5) print $16,"q.length",$4-1,$5,"+",$15,"r.length",$1-1,$2,$5-$4+1,$5-$4+1,'60',"cg:Z:"$5-$4+1"M","cs:Z::"$5-$4+1; else print $16,"q.length",$4-1,$5,"-",$15,"r.length",$1-1,$2,($5-$4+1)*-1,($5-$4+1)*-1,'60',"cg:Z:"($5-$4+1)*-1"M","cs:Z::"($5-$4+1)*-1}' > tmp1

#tmp1总共14列, 此时还有2列信息(第2列、、第7列)需要补充
bioawk -c fastx '{print $name, length($seq)}' ref.fasta > ref.length
bioawk -c fastx '{print $name, length($seq)}' query.fasta > query.length

#填充了第2列
python 1.py query.length tmp1  | awk '{$3="";$4="";print $0}' | sed 's/\s\+/ /g' | tr " " "\t" > tmp2
#填充了第7列
python 2.py ref.length tmp2 | awk 'BEGIN{ FS="\t";OFS="\t" }{print $3,$4,$5,$6,$7,$1,$2,$10,$11,$12,$13,$14,$15,$16}' > nucmer.paf

transanno minimap2chain nucmer.paf --output query2ref.nucmer.chain
