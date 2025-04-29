#!/bin/bash
set -e

# 参数
big_maf="$1"        # 大的maf文件，比如 5zokors.maf
ref_species="$2"    # 参考物种，比如 Efontanierii
query_list="$3"     # 查询物种，比如 "Ebaileyi,Esmithi,Erufescens,Erothschildi" 或 species.list文件
final_output="$4"   # 最终合并的输出，比如 5zokors.final.lst
block_per_file="${5:-5000}"  # 每多少个block切一份，默认5000
cpu_limit="${6:-8}"  # 默认最多并行8个线程，如果第6个参数指定了，就用指定的

# 检查
if [[ -z "$big_maf" || -z "$ref_species" || -z "$query_list" || -z "$final_output" ]]; then
    echo "Usage: bash run_convertMaf_parallel.sh big.maf ref_species query_list final_output [block_per_file] [CPU]"
    echo "Example: bash run_convertMaf_parallel.sh 5zokors.maf Efontanierii Ebaileyi,Esmithi,Erufescens,Erothschildi 5zokors.merged.lst 5000 8"
    exit 1
fi

# 创建临时工作目录
mkdir -p split_maf
mkdir -p split_output
mkdir -p split_tmp

# 切分maf
echo "Step1: Splitting $big_maf into blocks..."
block_count=0
file_count=1
out="split_maf/part_$(printf "%04d" $file_count).maf"
> "$out"

while IFS= read -r line; do
    if [[ "$line" =~ ^a\ score= ]]; then
        if (( block_count >= block_per_file )); then
            ((file_count++))
            out="split_maf/part_$(printf "%04d" $file_count).maf"
            > "$out"
            block_count=0
        fi
        ((block_count++))
    fi
    echo "$line" >> "$out"
done < "$big_maf"

echo "Total split into $file_count parts."

# 跑子脚本
echo "Step2: Running perl 01.convertMaf2List.pl in parallel with CPU limit $cpu_limit..."

# 当前并发数
current_jobs=0

for maf_part in split_maf/*.maf; do
    base=$(basename "$maf_part" .maf)
    perl 01.convertMaf2List.pl -i "$maf_part" \
                               -o "split_output/${base}.lst" \
                               -ref "$ref_species" \
                               -query "$query_list" \
                               -tmp "split_tmp/${base}.no_ref.maf" &
    ((current_jobs++))
    if [[ "$current_jobs" -ge "$cpu_limit" ]]; then
        wait -n   # 等任意一个子任务完成
        ((current_jobs--))
    fi
done

# 等待剩下的所有子进程
wait
echo "All parts finished."

# 合并输出
echo "Step3: Merging outputs..."

first=1
for file in split_output/*.lst; do
    if [[ $first -eq 1 ]]; then
        cat "$file" > "$final_output"
        first=0
    else
        tail -n +2 "$file" >> "$final_output"
    fi
done

echo "Final merged output: $final_output"

# 结束
echo "Done!"

