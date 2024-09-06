# path to reference virus
virus_path=$1

# name of virus
Virus=$2

# path to a file with a list of accessions, one accession per line
acc_list=$3

# Creating all of my folers and directories
mkdir work

mkdir work/virus_ref
mkdir work/plotting
mkdir work/plotting/original_depth
mkdir work/plotting/mismatches
mkdir work/bam_read_seperation/

mkdir work/read_start_end

mkdir work/infer_exp

mkdir work/plotting/single

mkdir work/plotting/single/positive_depth
mkdir work/plotting/single/negative_depth
mkdir work/plotting/single/negative_start
mkdir work/plotting/single/positive_start
mkdir work/plotting/single/negative_end
mkdir work/plotting/single/positive_end

mkdir work/plotting/paired

mkdir work/plotting/paired/positive_depth
mkdir work/plotting/paired/negative_depth
mkdir work/plotting/paired/negative_start
mkdir work/plotting/paired/positive_start
mkdir work/plotting/paired/negative_end
mkdir work/plotting/paired/positive_end

mkdir $Virus"_results"

mkdir work/mismatch
# downloading the script that is required to parse the output of bam-readcount for the mismatch plot
curl -o work/parse_brc.py "https://raw.githubusercontent.com/genome/bam-readcount/master/tutorial/scripts/parse_brc.py"
chmod +x work/parse_brc.py

mkdir work/fastq
mkdir work/trimmed



# while read p; do
# echo "$p"

# six_p="${p:0:6}"
# type="${p:0:3}"

# if [ "$type" = "ERR" ]; then
#     last="${p: -1}"
#     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$six_p/00$last/$p/*.gz
#     mv $p* work/fastq/
# elif [ "$type" = "SRR" ]; then
#     last="${p: -2}"
#     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$six_p/0$last/$p/*.gz
#     mv $p* work/fastq/
# elif [ "$type" = "DRR" ]; then
#     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$six_p/$p/*.gz
#     mv $p* work/fastq/
# fi

# done < $acc_list

# while read p; do
# echo "$p"

# pigz -d work/fastq/$p*

# done < $acc_list



# while read p; do
# echo "$p" 

# prefetch $p --max-size u 
# fasterq-dump --progress -O work/fastq/ $p 

# rm -r "$p"

# done < $acc_list



# -----------------------------------------------------------------------------------------------------------------------------------------

# First step is to download the fastq files

while read p; do
echo "$p" 

# Checking some conditions so that I don't download files unnecessarily
if [ -f work/fastq/$p"_1.fastq" ] || [ -f work/fastq/$p"_2.fastq" ] || [ -f work/fastq/$p".fastq" ]; then
    echo $p already downloaded

elif ([ -f work/trimmed/$p"_1_val_1.fq" ] && [ -f work/trimmed/$p"_2_val_2.fq" ]) || [ -f work/trimmed/$p"_trimmed.fq" ]; then
    echo $p already trimmed
    
else
    # enaBrowserTools is the tool that I am using to download the fastq files. It can run with aspera or with ftp
    # I have it running with ftp currently as that is faster

    # this is the option to use aspera
    # enaBrowserTools-1.7.1/python3/enaDataGet -f fastq -as aspera_settings.ini $p -d work/fastq/
    enaBrowserTools-1.7.1/python3/enaDataGet -f fastq $p -d work/fastq/ 
    mv work/fastq/$p/*.gz work/fastq/

    pigz -d work/fastq/$p*

    rm -r work/fastq/$p/

fi


done < $acc_list



# Trimming the fastq files and moving them to work/trimmed

while read p; do
echo "$p" 

if ([ -f work/trimmed/$p"_1_val_1.fq" ] && [ -f work/trimmed/$p"_2_val_2.fq" ]) || [ -f work/trimmed/$p"_trimmed.fq" ]; then
    echo $p already trimmed

elif [ -f work/fastq/$p"_1.fastq" ] && [ -f work/fastq/$p"_2.fastq" ]; then
    trim_galore --paired -o work/trimmed/ work/fastq/$p"_1.fastq" work/fastq/$p"_2.fastq"
    rm work/fastq/$p*
else
    trim_galore -o work/trimmed/ work/fastq/$p".fastq"
    rm work/fastq/$p*
fi

done < $acc_list


bowtie2-build --quiet -f $virus_path work/virus_ref/$Virus


# -----------------------------------------------------------------------------------------------------------------------------------------


filtering() {

    # -F 260 removes all unmapped and secondary mapping reads. Greatly reduces bam file size
    samtools view -S -b -@ 8 -F 260 $Virus"_"$p > $Virus"_"$p".bam"

    rm $Virus"_"$p

    samtools sort -@ 8 $Virus"_"$p".bam" > $Virus"_"$p"_sorted.bam"

    rm $Virus"_"$p".bam"

    samtools index $Virus"_"$p"_sorted.bam"
    
    # creating num_reads to count number of reads in the bam file.
    num_reads=$(samtools view -c $Virus"_"$p"_sorted.bam")
}

step_one() {

    # setting up my variables for subsequent use
    ref_name=$(samtools view $Virus"_"$p"_sorted.bam" | cut -f 3 | head -1)
    ref_len=$(seqkit fx2tab --length --name $virus_path | cut -f 2)

    # Creating my bed file for infer_experiment
    echo -e $ref_name"\t0\t"$ref_len"\tp\tp\t+\tp\tp\tp\tp\tp\tp" > work/temp.bed

    infer_experiment.py -r work/temp.bed -i $Virus"_"$p"_sorted.bam" > work/infer_exp/$p".infer_experiment.txt"

    rm work/temp.bed

    # bedtools genomecov -d -ibam is giving the coverage per position for the bam file.
    # Check if the variable is set to "paired" or if it is not and then run the appropriate bedtools genomecov parameters
    if [ "$reads_type" = "paired" ]; then
        bedtools genomecov -d -pc -ibam $Virus"_"$p"_sorted.bam" > work/plotting/original_depth/$p".tsv"
    else
        bedtools genomecov -d -ibam $Virus"_"$p"_sorted.bam" > work/plotting/original_depth/$p".tsv"
    fi


    #Getting the Mismatches Plot ready
    bam-readcount -w1 -f $virus_path $Virus"_"$p"_sorted.bam" > work/mismatch/$p".tsv"
    python work/parse_brc.py work/mismatch/$p".tsv" > work/mismatch/$p"_filtered.tsv"
    
    #quality filtering
    samtools view -@ 8 -b -q 10 $Virus"_"$p"_sorted.bam" > work/bam_read_seperation/$p"_filtered.bam"
}


# -----------------------------------------------------------------------------------------------------------------------------------------


single() {

# setting up my variables for subsequent use
tolerance=$(echo "$strand_one * 0.1" | bc -l)
lower_bound=$(echo "$strand_one - $tolerance" | bc -l)
upper_bound=$(echo "$strand_one + $tolerance" | bc -l)

# The first check is if the fastq file is unstranded. If it is then usually strand_one and strand_two are within 10% of each other
# If they are within 10% of each other, I want to default to ++,-- parameters.
if (( $(echo "$strand_two >= $lower_bound" | bc -l) && $(echo "$strand_two <= $upper_bound" | bc -l) )); then
    echo "strand_one: $strand_one"
    echo "strand_two: $strand_two"
    echo "++,--"

    samtools view -@ 8 -b -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense.bam"
    samtools view -@ 8 -b -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense.bam"

# Infer_experiment tells us that this is positive stranded, so parameters appropriate for that. 
elif (( $(echo "$strand_one > $strand_two" | bc -l) )); then
    echo "strand_one: $strand_one"
    echo "strand_two: $strand_two"
    echo "++,--"

    samtools view -@ 8 -b -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense.bam"
    samtools view -@ 8 -b -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense.bam"

# negative stranded
else
    echo "strand_one: $strand_one"
    echo "strand_two: $strand_two"
    echo "+-,-+"

    samtools view -@ 8 -b -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense.bam"
    samtools view -@ 8 -b -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense.bam"

fi 

samtools sort -@ 8 work/bam_read_seperation/$p"_sense.bam" > work/bam_read_seperation/$p"_sense_sorted.bam"
samtools sort -@ 8 work/bam_read_seperation/$p"_antisense.bam" > work/bam_read_seperation/$p"_antisense_sorted.bam"

rm work/bam_read_seperation/$p"_filtered.bam"
rm work/bam_read_seperation/$p"_sense.bam"
rm work/bam_read_seperation/$p"_antisense.bam"

}


paired() {

tolerance=$(echo "$strand_one * 0.1" | bc -l)
lower_bound=$(echo "$strand_one - $tolerance" | bc -l)
upper_bound=$(echo "$strand_one + $tolerance" | bc -l)

# If the file is unstranded, then I want to default to positive stranded parameters.
if (( $(echo "$strand_two >= $lower_bound" | bc -l) && $(echo "$strand_two <= $upper_bound" | bc -l) )); then
    echo "strand_one: $strand_one"
    echo "strand_two: $strand_two"
    echo "1++,1--,2+-,2-+"
    samtools view -@ 8 -b -f 0x40 -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense_read1.bam"

    samtools view -@ 8 -b -f 0x40 -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense_read1.bam"

    samtools view -@ 8 -b -f 0x80 -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense_read2.bam"

    samtools view -@ 8 -b -f 0x80 -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense_read2.bam"


# Positive strandedness
elif (( $(echo "$strand_one > $strand_two" | bc -l) )); then
    echo "strand_one: $strand_one"
    echo "strand_two: $strand_two"
    echo "1++,1--,2+-,2-+"
    samtools view -@ 8 -b -f 0x40 -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense_read1.bam"

    samtools view -@ 8 -b -f 0x40 -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense_read1.bam"

    samtools view -@ 8 -b -f 0x80 -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense_read2.bam"

    samtools view -@ 8 -b -f 0x80 -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense_read2.bam"


# Negative strandedness
else
    echo "strand_one: $strand_one"
    echo "strand_two: $strand_two"
    echo "1+-,1-+,2++,2--"

    samtools view -@ 8 -b -f 0x40 -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense_read1.bam"

    samtools view -@ 8 -b -f 0x40 -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense_read1.bam"

    samtools view -@ 8 -b -f 0x80 -F 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_sense_read2.bam"

    samtools view -@ 8 -b -f 0x80 -f 0x10 work/bam_read_seperation/$p"_filtered.bam" > work/bam_read_seperation/$p"_antisense_read2.bam"


fi

# combining read1 and read2 files into one file
samtools cat -o work/bam_read_seperation/$p"_sense.bam" work/bam_read_seperation/$p"_sense_read1.bam" work/bam_read_seperation/$p"_sense_read2.bam"
samtools cat -o work/bam_read_seperation/$p"_antisense.bam" work/bam_read_seperation/$p"_antisense_read1.bam" work/bam_read_seperation/$p"_antisense_read2.bam"

samtools sort -@ 8 work/bam_read_seperation/$p"_sense.bam" > work/bam_read_seperation/$p"_sense_sorted.bam"
samtools sort -@ 8 work/bam_read_seperation/$p"_antisense.bam" > work/bam_read_seperation/$p"_antisense_sorted.bam"


samtools sort -@ 8 work/bam_read_seperation/$p"_sense_read1.bam" > work/bam_read_seperation/$p"_sense_read1_sorted.bam"
samtools sort -@ 8 work/bam_read_seperation/$p"_antisense_read1.bam" > work/bam_read_seperation/$p"_antisense_read1_sorted.bam"

samtools sort -@ 8 work/bam_read_seperation/$p"_sense_read2.bam" > work/bam_read_seperation/$p"_sense_read2_sorted.bam"
samtools sort -@ 8 work/bam_read_seperation/$p"_antisense_read2.bam" > work/bam_read_seperation/$p"_antisense_read2_sorted.bam"


rm work/bam_read_seperation/$p"_sense.bam"
rm work/bam_read_seperation/$p"_antisense.bam"
rm work/bam_read_seperation/$p"_filtered.bam"
rm work/bam_read_seperation/$p"_sense_read1.bam"
rm work/bam_read_seperation/$p"_sense_read2.bam"
rm work/bam_read_seperation/$p"_antisense_read1.bam"
rm work/bam_read_seperation/$p"_antisense_read2.bam"


}


# -----------------------------------------------------------------------------------------------------------------------------------------


step_two_single() {
    # Creating the sense and antisense reads coverage plot
    bedtools genomecov -d -ibam work/bam_read_seperation/$p"_sense_sorted.bam" > work/plotting/single/positive_depth/$p".tsv"
    bedtools genomecov -d -ibam work/bam_read_seperation/$p"_antisense_sorted.bam" > work/plotting/single/negative_depth/$p".tsv"

    # Creating a temporary file for getting the start and cigar string of each read. 
    samtools view work/bam_read_seperation/$p"_sense_sorted.bam" | cut -f 4,6 > work/read_start_end/$p"_sense.tsv"
    samtools view work/bam_read_seperation/$p"_antisense_sorted.bam" | cut -f 4,6 > work/read_start_end/$p"_antisense.tsv"
}

step_two_paired() {
    # Creating the sense and antisense reads coverage plot
    bedtools genomecov -d -pc -ibam work/bam_read_seperation/$p"_sense_sorted.bam" > work/plotting/paired/positive_depth/$p".tsv"
    bedtools genomecov -d -pc -ibam work/bam_read_seperation/$p"_antisense_sorted.bam" > work/plotting/paired/negative_depth/$p".tsv"

    # Creating a temporary file for getting the start and cigar string of each read. 
    samtools view work/bam_read_seperation/$p"_sense_read1_sorted.bam" | cut -f 4,6 > work/read_start_end/$p"_sense_read1.tsv"
    samtools view work/bam_read_seperation/$p"_antisense_read1_sorted.bam" | cut -f 4,6 > work/read_start_end/$p"_antisense_read1.tsv"

    samtools view work/bam_read_seperation/$p"_sense_read2_sorted.bam" | cut -f 4,6 > work/read_start_end/$p"_sense_read2.tsv"
    samtools view work/bam_read_seperation/$p"_antisense_read2_sorted.bam" | cut -f 4,6 > work/read_start_end/$p"_antisense_read2.tsv"
}


# -----------------------------------------------------------------------------------------------------------------------------------------



while read p; do
echo $p

# I need to first distinguish between single and paired libraries, which is why I have run these if-else statements each time.
if [ -f work/trimmed/$p"_1_val_1.fq" ] && [ -f work/trimmed/$p"_2_val_2.fq" ]; then
    # read mapping
    time bowtie2 --local --very-sensitive-local --threads 8 -k 10 -x work/virus_ref/$Virus -1 work/trimmed/$p"_1_val_1.fq" -2 work/trimmed/$p"_2_val_2.fq" -S  $Virus"_"$p
    
    filtering

    if [ "$num_reads" -eq 0 ]; then
        sed -i '/'$p'/d' $acc_list

        echo $Virus is not present in Library $p as 0 reads mapped to your reference virus. Removing accession from accession file.
    else

        reads_type="paired"

        step_one

        # reading the infer_experiment file and using regex, getting the percentages for each strand.
        strand_one=$(grep '1++,1--,2+-,2-+' work/infer_exp/$p".infer_experiment.txt" | sed 's/.*1++,1--,2+-,2-+": \([0-9.]*\).*/\1/')
        strand_two=$(grep '1+-,1-+,2++,2--' work/infer_exp/$p".infer_experiment.txt" | sed 's/.*1+-,1-+,2++,2--": \([0-9.]*\).*/\1/')

        paired

        step_two_paired
    fi

    
    
else
    # read mapping
    time bowtie2 --local --very-sensitive-local --threads 8 -k 10 -x work/virus_ref/$Virus -U work/trimmed/$p"_trimmed.fq" -S  $Virus"_"$p

    filtering

    if [ "$num_reads" -eq 0 ]; then
        sed -i '/'$p'/d' $acc_list

        echo $Virus is not present in Library $p as 0 reads mapped to your reference virus. Removing accession from accession file.
    else

        reads_type="single"

        step_one

        # reading the infer_experiment file and using regex, getting the percentages for each strand.
        strand_one=$(grep '++,--' work/infer_exp/$p".infer_experiment.txt" | sed 's/.*++,--": \([0-9.]*\).*/\1/')
        strand_two=$(grep '+-,-+' work/infer_exp/$p".infer_experiment.txt" | sed 's/.*+-,-+": \([0-9.]*\).*/\1/')

        single

        step_two_single
    fi
fi


done < $acc_list


# -----------------------------------------------------------------------------------------------------------------------------------------


# run python file with the necessary variables
python final_script2_new.py $acc_list $Virus