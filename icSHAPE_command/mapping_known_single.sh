input_dir="/data1/bioinfo/zhoutong/data/RASP/PRJNA608297_icSHAPE2021b/fastp_data/"
output_dir="/data1/bioinfo/zhoutong/data/RASP/PRJNA608297_icSHAPE2021b/mapping_known/human_rbp/"
index="/data1/bioinfo/zhoutong/data/RASP/human_annotation_files/index_human/human_RBP/human_RBP"
threads=100

for fastp_file in "${input_dir}"*.fq; do
    basename=$(basename "$fastp_file" .fq)

    bowtie2 -x "${index}" -U "${fastp_file}" -S "${output_dir}${basename}.sam" \
    -p "${threads}" --local 2>"${output_dir}${basename}_Align.summary"
done
