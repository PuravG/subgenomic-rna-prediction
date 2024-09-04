# path to reference virus
virus_path='/drive/13xrrna/Rodent_hepacivirus.fasta'

# name of virus
Virus='Rodent_hepacivirus'
Virus=$(tr -s ' ' '_' <<< $Virus)

# path to a file with a list of accessions, one accession per line
acc_list='/drive/13xrrna/hep_rodent_acc.txt'

./final_script1_new.sh $virus_path $Virus $acc_list