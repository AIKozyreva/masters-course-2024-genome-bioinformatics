#!/bin/bash

input_storage="/mnt/projects/users/aalayeva/genomics/raw"
output_dir="$(pwd)"
working_dir="$(pwd)"

# Check if output_dir exists, if not, create it
#if [ ! -d "$working_dir/$output_dir" ]; then
#    mkdir -p "$working_dir/$output_dir"
#fi

for file in $input_storage/*.{fa.gz,fna.gz,fastq.gz,fasta.gz}; do
    gzip -d "$file"
    echo "Decompressed: $file"
    echo "\n"
done

# Process each file in input_storage folder
for file_path in $input_storage/*.{fa,fna,fastq,fasta}; do
    if [ -e "$file_path" ]; then
        file_name=$(basename "$file_path")
        file_size=$(du -h "$file_path" | cut -f1)
        file_line_count=$(wc -l < "$file_path")
        file_header_count=$(grep -c '^>' "$file_path")
        file_header_info=$(grep '^>' "$file_path" | head -n 1 | sed 's/^>//')

        # Output to console and save to out_gen_info.txt
        echo "Processing $file_name"
        echo "path to file is $file_path"
        echo "File size is $file_size"
        echo "$file_name amount of lines: $file_line_count"
        echo "$file_name amount of headers: $file_header_count"
        echo "$file_name header information: $file_header_info"
        echo "\n"

        # Save to out_gen_info.txt
        echo "Processing $file_name" >> "$working_dir/out_gen_info.txt"
        echo "path to file is $file_path" >> "$working_dir/out_gen_info.txt"
        echo "File size is $file_size" >> "$working_dir/out_gen_info.txt"
        echo "$file_name amount of lines: $file_line_count" >> "$working_dir/out_gen_info.txt"
        echo "$file_name amount of headers: $file_header_count" >> "$working_dir/out_gen_info.txt"
        echo "$file_name header information: $file_header_info" >> "$working_dir/out_gen_info.txt"
        echo "" >> "$working_dir/out_gen_info.txt"
    fi
done
