# path to reference virus
virus_path='/drive/13xrrna/Hepatitis_A_virus.fasta'

# name of virus
Virus='Hepatitis_A'
Virus=$(tr -s ' ' '_' <<< $Virus)

# path to a file with a list of accessions, one accession per line
acc_list='/drive/13xrrna/Hepatitis_A_Virus_acc.txt'

./final_script1_new.sh $virus_path $Virus $acc_list