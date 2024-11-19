#!/bin/bash

# 指定SAM文件所在目录
sam_directory="/home/tmp/data/a375/18s_sam"

# 定义输出文件夹
output_directory="/home/tmp/data/a375/02_count_stop"

# 处理SAM文件并计算RT值
for sam_file in "$sam_directory"/*.sam; do
    # 获取SAM文件名（不包含路径）
    sam_filename=$(basename "$sam_file")
    
    # 构建对应的TAB文件路径
    tab_file="$output_directory/${sam_filename%.sam}.tab"
    
    # 执行icSHAPE-pipe sam2tab命令
    echo "Running: icSHAPE-pipe sam2tab -in $sam_file -out $tab_file"
    icSHAPE-pipe sam2tab -in "$sam_file" -out "$tab_file"
    
    # 构建对应的输出文件路径
    output_file="$output_directory/${sam_filename%.sam}_countRT.c"
    
    # 执行icshape countRT命令
    echo "Running: icSHAPE-pipe countRT -in $tab_file -size /home/tmp/data/a375/human_rRNA_tRNA_mtRNA.len -out $output_file"
    icSHAPE-pipe countRT -in "$tab_file" -size "/home/tmp/data/a375/human_rRNA_tRNA_mtRNA.len" -out "$output_file"
    
    # 构建对应的SHAPE文件路径
    shape_file="$output_directory/${sam_filename%.sam}_shape.gTab"

    # 执行icSHAPE-pipe calcSHAPENoCont命令
    echo "Running: icSHAPE-pipe calcSHAPENoCont -N $tab_file -size /home/tmp/data/a375/human_rRNA_tRNA_mtRNA.fa -genome /home/tmp/data/a375/srp_7sk/7sk_srp.fasta -bases A,T,C,G -out $shape_file -non-sliding"
    icSHAPE-pipe calcSHAPENoCont -N "$tab_file" -size /home/tmp/data/a375/human_rRNA_tRNA_mtRNA.len -genome /home/tmp/data/a375/human_rRNA_tRNA_mtRNA.fa -bases A,T,C,G -out "$shape_file" -non-sliding
    
    # 将每个输出文件转换为CSV格式
    # 获取文件名（不包含路径）
    result_filename=$(basename "$output_file")
    
    # 构建对应的CSV文件路径
    csv_file="${output_file%.c}.csv"
    
    # 执行sed命令将TAB文件转换为CSV文件
    echo "Converting to CSV: sed 's/\t/,/g' $output_file > $csv_file"
    sed 's/\t/,/g' "$output_file" > "$csv_file"
done