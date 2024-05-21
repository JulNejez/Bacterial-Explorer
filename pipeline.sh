#!/bin/bash

# Variables from widget
input_file=$4
reference_path=$3
py_directory='/julie/static'
methods=$1
threshold=$2
analysis=$5

echo "$methods"


#Definding working directory
working_directory='/julie/static/process_folder'
echo $working_directory
cd "$working_directory"

# Make dir to output
mkdir "$working_directory/Bacterial_Explorer"
mkdir "$working_directory"/final_output

# Copy README.txt to Bacterial_Explorer
cp "$py_directory"/README.txt Bacterial_Explorer/README.txt


# Make multifasta
cat "$reference_path"/*.fasta > "$working_directory"/multifasta.fasta

#Do text file with path to the genomes (genomes_list)
find "$reference_path" -type f -name "*.fasta" -exec readlink -f {} \; > "$working_directory"/genomes_list.txt

file_name="genomes_list.txt"
ref="$working_directory/$file_name"



if [[ "$analysis" == "16SrRNA" ]]; then

  #Barrnap
  barrnap -i "$input_file" -kingdom bac -evalue 1e-6 -lencutoff 0.8 -reject 0.25 -outseq "$working_directory"/rRNA_all.fasta -threads 6
  output_rRNA="$working_directory/rRNA_references.fasta"
  > "$output_rRNA"

  if [ "$reference_path" == "julie/static/own_database" ]; then
    for file in "$reference_path"/*; do
      if [[ -f "$file" ]]; then       
        barrnap -kingdom bac -evalue 1e-6 -lencutoff 0.8 -reject 0.25 -outseq "${file}.temp.fasta" -threads 6 "$file"
          cat "${file}.temp.fasta" >> "$output_rRNA"
          rm "${file}.temp.fasta"
      fi
    done
    grep -A 1 '^>16S_rRNA' rRNA_references.fasta | grep -v '^--' > "$working_directory"/16SrRNA_references.fasta
  fi



  # Select only 16S rRNA
  grep -A 1 '^>16S_rRNA' rRNA_all.fasta | grep -v '^--' > "$working_directory"/rRNA.fasta
  #grep -A 1 '^>16S_rRNA' rRNA_references.fasta | grep -v '^--' > "$working_directory"/16SrRNA_references.fasta

  # Remove .fai files
  rm "$reference_path"/*.fai

  # Is rRNA  empty?
  if [ -s "$working_directory"/rRNA.fasta ]; then
    echo "Known rRNAs have been identified!"
    rRNA_file="$working_directory"/rRNA.fasta

    if [ "$reference_path" == "julie/static/own_database" ]; then
      rRNA_references="$working_directory"/16SrRNA_references.fasta
    else
      rRNA_references="$py_directory"/16SrRNA_references.fasta
    fi

  else
    echo "No known rRNA has been identified!"
    exit 0
  fi

  ### Cd-hit
  if [ "$methods" = "cd_hit" ] || [ "$methods" = "cd_fast" ] || [ "$methods" = "cd_blat" ] || [ "$methods" = "cd_blast" ] || [ "$methods" = "all" ] || [ "$methods" = "cd_fast_blat" ] || [ "$methods" = "cd_blast_blat" ]; then

    cd-hit-est-2d -i "$rRNA_file" -i2 "$rRNA_references" -o "$working_directory"/cd-hit_output -T 6 -M 0 -c $threshold

    # Analysis of cd-hit results
    python3 "$py_directory"/cd_hit_analysis_new.py "$ref"
    cp cd-hit_output.clstr Bacterial_Explorer/cd-hit_output.clstr
    cp cd-hit.txt "$working_directory"/final_output/cd-hit.txt
  fi

  # If method cd_hit, save cd-hit.txt to output.txt
  if [ "$methods" = "cd_hit" ]; then
    cp cd-hit.txt output_terminal.txt
    cp output_terminal.txt Bacterial_Explorer/output.txt
  fi
    
  ### BLAST
  if [ "$methods" = "blast" ] || [ "$methods" = "cd_blast" ] || [ "$methods" = "fast_blast" ] || [ "$methods" = "all" ] || [ "$methods" = "fast_blast_blat" ] || [ "$methods" = "cd_blast_blat" ] || [ "$methods" = "blast_blat" ]; then
    makeblastdb -in "$rRNA_references" -dbtype nucl -out blast_database
    blastn -query "$rRNA_file" -db blast_database -out "$working_directory"/blast_output.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 6 -perc_identity $threshold
      
    # Analysis of BLAST results
    python3 "$py_directory"/blast_analysis_16SrRNA.py "$ref" "$threshold"
    cp blast_output.txt Bacterial_Explorer/blast_output.txt
    cp blast.txt "$working_directory"/final_output/blast.txt
  fi

  if [ "$methods" = "blast" ]; then
    cp blast.txt output_terminal.txt
    cp output­_terminal.txt Bacterial_Explorer/output.txt
  fi

  ### BLAT
  cd "$working_directory"
  if [ "$methods" = "blat" ] || [ "$methods" = "all" ] || [ "$methods" = "cd_blat" ] || [ "$methods" = "fast_blat" ] || [ "$methods" = "cd_fast_blat" ] || [ "$methods" = "fast_blast_blat" ] || [ "$methods" = "cd_blast_blat" ] || [ "$methods" = "blast_blat" ]; then
    makeblastdb -in "$rRNA_references" -dbtype nucl -out reference_db
    blastdbcmd -db reference_db -entry all -outfmt "%f" -out converted_database.fa
    blat converted_database.fa "$rRNA_file" blat_output.blast -out=blast8 -minMatch=50
    python3 "$py_directory"/blat_analysis_new.py "$threshold" "$ref"
    cp blat_output.blast Bacterial_Explorer/blat_output.blast
    cp blat.txt "$working_directory"/final_output/blat.txt
  fi

  if [ "$methods" = "blat" ]; then
    cp blat.txt output_terminal.txt
    cp output_terminal.txt Bacterial_Explorer/output.txt
  fi

elif [[ "$analysis" == "whole" ]]; then
    
  ### Cd-hit
  cd_hit_input_file="$input_file"
  cd_hit_references="$working_directory"/multifasta.fasta

  if [ "$methods" = "cd_hit" ] || [ "$methods" = "cd_fast" ] || [ "$methods" = "cd_blat" ] || [ "$methods" = "cd_blast" ] || [ "$methods" = "all" ] || [ "$methods" = "cd_fast_blat" ] || [ "$methods" = "cd_blast_blat" ]; then

    cd-hit-est-2d -i "$cd_hit_input_file" -i2 "$cd_hit_references" -o "$working_directory"/cd-hit_output -T 6 -M 0 -c $threshold

    # Analysis of cd-hit results
    python3 "$py_directory"/cd_hit_analysis_new.py "$ref"
    cp cd-hit_output.clstr Bacterial_Explorer/cd-hit_output.clstr
    cp cd-hit.txt "$working_directory"/final_output/cd-hit.txt
  fi

  # If method cd_hit, save cd-hit.txt to output.txt
  if [ "$methods" = "cd_hit" ]; then
    cp cd-hit.txt output_terminal.txt
    cp output_terminal.txt Bacterial_Explorer/output.txt
  fi


  ### FastANI
  if [ "$methods" = "fastAni" ] || [ "$methods" = "cd_fast" ] || [ "$methods" = "fast_blast" ] || [ "$methods" = "all" ] || [ "$methods" = "cd_fast_blat" ] || [ "$methods" = "fast_blast_blat" ] || [ "$methods" = "fast_blat" ]; then
    cd "$reference_path"
    fastANI -q "$input_file" --rl "$working_directory"/genomes_list.txt -o "$working_directory"/fast_ANI_output.txt -t 6
    cd "$working_directory"

    # Analysis of fastANI
    python3 "$py_directory"/fastANI_analysis_new.py "$threshold"
    cp fast_ANI_output.txt Bacterial_Explorer/fast_ANI_output.txt
    cp fastANI.txt "$working_directory"/final_output/fastANI.txt
  fi

  # If method is fastANI
  if [ "$methods" = "fastAni" ]; then
    cp fastANI.txt output_terminal.txt
    cp output_terminal.txt Bacterial_Explorer/output.txt
  fi

  ### BLAST
  cd "$working_directory"
  blast_input="$input_file"
  echo "$threshold" > threshold_value.txt
  blast_reference="$working_directory"/multifasta.fasta
  if [ "$methods" = "blast" ] || [ "$methods" = "cd_blast" ] || [ "$methods" = "fast_blast" ] || [ "$methods" = "all" ] || [ "$methods" = "fast_blast_blat" ] || [ "$methods" = "cd_blast_blat" ] || [ "$methods" = "blast_blat" ]; then
    makeblastdb -in "$blast_reference" -dbtype nucl -out blast_database
    blastn -query "$blast_input" -db blast_database -out "$working_directory"/blast_output_orig.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 6 -perc_identity $threshold
    cp blast_output_orig.txt blast_out.txt
    awk '$4 >= 3000' blast_out.txt > blast_output.txt

    # Analysis of BLAST results
    python3 "$py_directory"/blast_analysis_whole.py "$ref" "$threshold"
    cp blast_output.txt Bacterial_Explorer/blast_output.txt
    cp blast.txt "$working_directory"/final_output/blast.txt
  fi

  if [ "$methods" = "blast" ]; then
    cp blast.txt output_terminal.txt
    cp output­_terminal.txt Bacterial_Explorer/output.txt
  fi

  ### BLAT
  cd "$working_directory"
  if [ "$methods" = "blat" ] || [ "$methods" = "all" ] || [ "$methods" = "cd_blat" ] || [ "$methods" = "fast_blat" ] || [ "$methods" = "cd_fast_blat" ] || [ "$methods" = "fast_blast_blat" ] || [ "$methods" = "cd_blast_blat" ] || [ "$methods" = "blast_blat" ]; then
    makeblastdb -in "$working_directory"/multifasta.fasta -dbtype nucl -out reference_db
    blastdbcmd -db reference_db -entry all -outfmt "%f" -out converted_database.fa
    blat converted_database.fa "$input_file" blat_output.blast -out=blast8 -minMatch=300
    python3 "$py_directory"/blat_analysis_new.py "$threshold" "$ref"
    cp blat_output.blast Bacterial_Explorer/blat_output.blast
    cp blat.txt "$working_directory"/final_output/blat.txt
  fi

  if [ "$methods" = "blat" ]; then
    cp blat.txt output_terminal.txt
    cp output_terminal.txt Bacterial_Explorer/output.txt
  fi

else
  echo "There is some problem!"
fi

### Final analysis
# If output_terminal.txt -> stop
cd "$working_directory"
if [ -f "output_terminal.txt" ]; then
  exit 0
else
  output_file="output_terminal.txt"
  all_bacteria=()
  
  for file in "$working_directory"/final_output/*; do
    if [ -f "$file" ]; then
      all_bacteria+=($(sort "$file"))
    fi
  done

  # Save unique bacteria
  unique_bacteria=($(printf "%s\n" "${all_bacteria[@]}" | sort -u))
  common_bacteria=()

  for bacteria in "${unique_bacteria[@]}"; do
    count=$(grep -l "$bacteria" "$working_directory"/final_output/* | wc -l)
    
    if [ "$count" -eq "$(ls -1 "$working_directory"/final_output | wc -l)" ]; then
      common_bacteria+=("$bacteria")
    fi
  done

  printf "%s\n" "${common_bacteria[@]}" > "$output_file"
  
  output_file="output.txt"
  for file in "$working_directory"/final_output/*; do
    if [ -f "$file" ]; then
      filename=$(basename "$file" .txt)  
      echo "$filename" >> "$output_file"
      echo "" >> "$output_file"  
      cat "$file" >> "$output_file"  
      echo "" >> "$output_file"  
      echo "------------------------------------------" >> "$output_file" 
      echo "" >> "$output_file"  
    fi
  done
  cp output.txt Bacterial_Explorer/output.txt
fi

exit 0



