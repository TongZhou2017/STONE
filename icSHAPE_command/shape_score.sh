#!/bin/bash

# 指定输入文件夹
input_directory="/home/tmp/data/a375/02_count_stop"

# 指定输出文件夹
output_directory="/home/tmp/data/a375/03_shape_score"

# 指定参考序列长度文件路径
ref_length_file="/home/tmp/data/a375/human_rRNA_tRNA_mtRNA.len"

# 遍历所有 .tab 文件并执行 genSHAPEToTransSHAPE 命令
for input_file in "$input_directory"/*.gTab; do
    # 获取文件名（不包含路径和扩展名）
    filename=$(basename -- "$input_file")
    filename_no_extension="${filename%.*}"
    
    # 执行 genSHAPEToTransSHAPE 命令
    echo "Running: icSHAPE-pipe genSHAPEToTransSHAPE -i $input_file -o $output_directory/${filename_no_extension}.out -s $ref_length_file -p 80"
    icSHAPE-pipe genSHAPEToTransSHAPE -i "$input_file" -o "$output_directory/${filename_no_extension}.out" -s "$ref_length_file" -p 80
done
