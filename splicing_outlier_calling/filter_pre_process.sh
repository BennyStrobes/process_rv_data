
pre_process_output_dir="$1"
snaptron_gtex_samples_file="$2"
filter_pre_process_output_dir="$3"

python filter_pre_process.py $pre_process_output_dir $snaptron_gtex_samples_file $filter_pre_process_output_dir