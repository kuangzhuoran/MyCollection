#!/bin/bash
set -e

# 参数
lst_file="$1"     # maf.lst文件
gff_file="$2"     # gff注释文件
outdir="$3"       # 输出目录
cpu_limit="${4:-8}"   # 并发CPU数，默认8

if [[ -z "$lst_file" || -z "$gff_file" || -z "$outdir" ]]; then
    echo "Usage: bash run_lst2gene_parallel.sh input.maf.lst input.gff output_dir [CPU]"
    echo "Example: bash run_lst2gene_parallel.sh 5zokors.maf.lst Efontanierii.sorted.rename.gff genes_out 8"
    exit 1
fi

# 创建输出目录
mkdir -p "$outdir"

# 预处理：把gff里面的CDS位置全部提取出来
echo "Step1: Preprocessing GFF file..."
awk '$3=="CDS" {split($9,a,";"); for(i in a) if(a[i]~"Parent=") {sub("Parent=","",a[i]); print $1,$4,$5,$6,a[i]}}' "$gff_file" > cds_info.tmp

# 预处理：把maf.lst读出来，建立位点与物种对应关系
echo "Step2: Preprocessing LST file..."
head -n 1 "$lst_file" > header.tmp
tail -n +2 "$lst_file" > lst_body.tmp

# 解析maf.lst文件，找到每个位点属于哪个基因
echo "Step3: Mapping sites to genes..."

# awk超快速处理
awk 'NR==FNR {pos[$1":"$2]=$5; next} 
     {key=$1":"$2; if(key in pos) print pos[key],$0}' cds_info.tmp lst_body.tmp > mapped_sites.tmp

# 拆分每个基因的记录
echo "Step4: Splitting per-gene data..."

mkdir -p gene_tmp
awk '{print > "gene_tmp/"$1".lst"}' mapped_sites.tmp

echo "Found $(ls gene_tmp | wc -l) genes with data."

# 并行处理每个基因
echo "Step5: Running parallel generation..."

current_jobs=0

for gene_lst in gene_tmp/*.lst; do
    gene_id=$(basename "$gene_lst" .lst)
    (
        # 构建fasta文件
        species=($(head -n 1 header.tmp | awk '{for(i=3;i<=NF;i++) print $i}'))
        {
            for sp in "${species[@]}"; do
                echo ">$sp"
                awk -v sp="$sp" '
                    BEGIN{FS=OFS="\t"}
                    {
                        for(i=3;i<=NF;i++) {
                            if(FNR==1) head[i]=$i
                            else if(head[i]==sp) {
                                if($i=="") $i="-";
                                printf "%s",toupper($i)
                            }
                        }
                        printf "\n"
                    }
                ' header.tmp "$gene_lst" | tail -n +2 | tr -d '\n'
                echo ""
            done
        } > "$outdir/$gene_id.fa"
    ) &

    ((current_jobs++))
    if [[ "$current_jobs" -ge "$cpu_limit" ]]; then
        wait -n
        ((current_jobs--))
    fi
done

wait

echo "All genes finished!"

# 清理临时文件
rm -rf gene_tmp cds_info.tmp lst_body.tmp mapped_sites.tmp header.tmp

echo "Done!"

