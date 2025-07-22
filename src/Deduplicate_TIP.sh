#!/usr/bin/env bash

# 显示帮助信息
show_help() {
    echo "Usage: $0 -i <input_bam> -o <output_bam>"
    echo
    echo "Arguments:"
    echo "  -i    Input BAM file (unsorted)"
    echo "  -o    Output deduplicated BAM file"
}

# 解析命令行参数
while getopts "i:o:" opt; do
    case $opt in
        i) input_bam=$OPTARG ;;
        o) output_bam=$OPTARG ;;
        *) show_help; exit 1 ;;
    esac
done

# 参数检查
if [[ -z "$input_bam" || -z "$output_bam" ]]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# 提取前缀
basename=$(basename "$input_bam" .bam)

# 创建临时目录
tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

# 排序并生成 index
samtools sort "$input_bam" -o "$tmpdir/${basename}.sort.bam"
samtools index "$tmpdir/${basename}.sort.bam"

# 提取 SAM 并分离 header / body
samtools view -h "$tmpdir/${basename}.sort.bam" > "$tmpdir/sort.sam"
grep '^@' "$tmpdir/sort.sam" > "$tmpdir/header.sam"
grep -v '^@' "$tmpdir/sort.sam" > "$tmpdir/body.sam"

# 去重 based on 5' end
awk '{
    flag = $2
    chr = $3
    start = $4
    seq = $10
    len = length(seq)

    if (and(flag, 16)) {
        strand = "-"
        pos5 = start + len - 1
    } else {
        strand = "+"
        pos5 = start
    }

    key = chr"\t"strand"\t"pos5
    print key"\t"len"\t"$0
}' "$tmpdir/body.sam" \
| sort -k1,1 -k2,2 -k3,3n -k4,4nr \
| awk '!seen[$1"\t"$2"\t"$3]++ { for(i=5; i<=NF; i++) printf "%s%s", $i, (i==NF ? ORS : OFS) }' OFS="\t" \
> "$tmpdir/unique.sam"

# 合并 header 和 unique body
cat "$tmpdir/header.sam" "$tmpdir/unique.sam" > "$tmpdir/final.sam"

# 生成最终 clean BAM
samtools view -bS "$tmpdir/final.sam" > "$tmpdir/clean_TIP.bam"
samtools sort "$tmpdir/clean_TIP.bam" -o "$output_bam"
samtools index "$output_bam"

# 生成去重报告
all_reads=$(samtools view -c "$tmpdir/${basename}.sort.bam")
nondup_reads=$(samtools view -c "$output_bam")
duplicate_ratio=$(echo "1-($nondup_reads/$all_reads)" | bc -l)

report_file="${output_bam%.bam}.TIP_dedup_report.txt"
{
  echo -e "Sample:\t$basename"
  echo -e "All_reads:\t$all_reads"
  echo -e "Nondup_reads:\t$nondup_reads"
  echo -e "Duplicate_ratio:\t$duplicate_ratio"
} > "$report_file"

echo "Deduplication complete. Output: $output_bam"
echo "Report written to: $report_file"
