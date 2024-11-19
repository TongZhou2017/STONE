input_dir="/data1/bioinfo/zhoutong/zhoutong/data/RASP/PRJNA301635"
output_dir="/data1/bioinfo/zhoutong/zhoutong/data/RASP/PRJNA301635/fastp_data"

for input_file_R1 in "$input_dir"/*1.fastq.gz; do
    file_name=$(basename "${input_file_R1}" | sed 's/1.fastq.gz//')
    input_file_R2="${input_dir}/${file_name}2.fastq.gz"

    fastp -i "$input_file_R1" -I "$input_file_R2" -o "$output_dir/${file_name}_R1_fastp.fq.gz" -O "$output_dir/${file_name}_R2_fastp.fq.gz" --dedup -q 20 -u 20 -n 5 -U --umi_loc=per_read --umi_len=5 --adapter_fasta adapter.fasta -h "$output_dir/${file_name}_fastp.html"
done
