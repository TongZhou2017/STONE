input_dir="/data1/bioinfo/zhoutong/zhoutong/data/RASP/PRJNA301635/fastp_data/"
output_dir="/data1/bioinfo/zhoutong/zhoutong/data/RASP/PRJNA301635/mapping_kown/"
index="/data1/bioinfo/zhoutong/zhoutong/data/RASP/PRJNA301635/bw_index/5s_SRP_U1"
threads=100

for fastp_file in "${input_dir}"*R1_fastp.fq.gz; do
    basename=$(basename "$fastp_file" _R1.fq)
    fastp_file_R1="${input_dir}${basename}_R1_fastp.fq.gz"
    fastp_file_R2="${input_dir}${basename}_R2_fastp.fq.gz"

    bowtie2 -x "${index}" -1 "${fastp_file_R1}" -2 "${fastp_file_R2}" -S "${output_dir}${basename}.sam" \
    -p "${threads}" --local 2>"${output_dir}${basename}_Align.summary"
done
