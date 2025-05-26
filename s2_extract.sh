#!/bin/bash

# --------------------------------------------------------
# questo script è diviso in due 6 parti
#     1) mi dici quale cartella ha il file ha file Genomici 
#     2) mi estrae alcune informazioni sui vari batch 
#        (densita come si chiama il file con i genotipi ecc...)
#     3) da le info del file punto due mi trasforma i file in formati plink 
#     4) creo panelli hd e md 
#     5) mi cero un altro file per ciascun animale attacco alcune info e sopratutto mi estraggo il correspettivo ID e MATRICOLA e RAZZA
#-----------------------------------------------------------

############# INPUT ############################
#     ANA0622.DBF    /check_parentage/SCRIPT  pedda_row.py  
#     /Scarico_GENOTI
#     /mappe_reference    /script
#----------------------------------------------------------#
#############################################
#         OUTPUT 
#     (file info)
#  info_batch_all.csv
#  info_batch.csv
#  all_ani_info.csv
#  n_hd;m_md
#----- (file da usare)---------------
#   panelli mergiati hd e md           
########################################

pedda_row_dir="/mnt/user_data/enrico/code_and_sftw/pedda_row.py"

#ls *.zip|awk -F'.zip' '{print "unzip "$0" -d "$1}'|sh

printf "#########################################\n
                      parte uno \n batch info se c'è il file zip
##########################################"



echo "batch;ZIP;geno" > info_batch_all.csv

cd Scarico_*

v="$(ls -d * )"
for i in $v 
do
  cd  $i
     MAP="$(find . -iname "*_SNP_Map.csv" -print)"
     ZIP="$(find . -iname "*.zip" -print)"
	  if test -f "$MAP"; then
        gen="si"
        else
          gen="no"
    fi
    if test -f "$ZIP"; then
          zip="si"
      else
          zip="no"
    fi
    echo $i";"$zip";"$gen
    cd ..
done  >> ../info_batch_all.csv

cd ..

printf "###########################################################
 parte due  estrai dentro ciascuna cartella i propi info *.zip
##########################################################"


awk -F ";" '{if($2=="no" && $3=="no") print $1}' info_batch_all.csv > inutili
awk -F ";" '{if($2=="si" ) print $1}' info_batch_all.csv > batch_to_unzip
#gawk -F ";" '{if($3=="si" ) print $1}' info_batch_all.csv > batch_to_unzip

n_batch="$(awk 'END{print NR}' batch_to_unzip )"
sequenza="$(seq -s " " 1 $n_batch)"

# NB..............
# future implementazione estri e aggrega solo quello che mi serve

cd Scarico_*

for it in $sequenza
  do
    i=$(awk -v nrow=$it 'NR==nrow' ../batch_to_unzip) # la matricola del file in variabile
    cd $i
    echo "-----------------------"$i"----------------------"

    zip="$(ls | grep zip)"
    unzip -o $zip
    cd ..
  done
cd ..


printf "###########################################################
 parte tre  converti con peddarow
##########################################################"

echo "batch;zip_name;map_file;file_snp;nid;nsnp;plink_fmt" > info_batch.csv

n_batch="$(awk 'END{print NR}' batch_to_unzip )"
sequenza="$(seq -s " " 1 $n_batch)"

# controlla a mano !!!!!
# qua 
# implementate se c'è zip srgna che c'+ se mpm so rstrae èrpbrlmi 

cd Scarico_*

for it in $sequenza
    do
    i=$(awk -v nrow=$it 'NR==nrow' ../batch_to_unzip) # la matricola del file in variabile
    cd $i
    echo "-----------------------"$i"----------------------"
    MAP=$(ls *SNP_Map.csv)
    full="$(find . -iname "*fullcompact.csv" -print) "
    final_rep="$(find . -iname "*FinalReport.csv" -print)"

    if [[ $full =~ "mpact" ]]; then
           FILE_SNP="$(echo $full | sed 's|./||g' )"
          echo "compact"
    elif  [ -e $final_rep  ]; then            #lo cambio
             FILE_SNP="$(ls *_FinalReport.csv)"
             echo   "final report"
    else
        echo "not found..."
    fi

  nid=$(awk 'END{print NR}' *_RelMatr.csv  )
  nsnp=$(awk 'END{print NR}' $MAP )

  plinkped=$(find . -iname "*PLINK*" -print)
  if [[ -e $plinkped ]]; then
       plink_fmt="si"
  elif [[ ! -e $plinkped ]]; then
       plink_fmt="no"
  fi
  echo "di densità: "$nsnp
  echo $i";"$zip";"$MAP";"$FILE_SNP";"$nid";"$nsnp";"$plink_fmt >>  ../../info_batch.csv
  cd ..
done

cd ..

# piccole statistiche descrittive
sed -i '1d' info_batch.csv
awk -F ";" 'NR > 1 {sum+=$5;} END{print "number of animals: " sum;}' info_batch.csv
echo "gentotipizazioni Roberto: "
grep  "MANTOVANI" info_batch.csv | awk 'END {print "number batch: " NR}'
grep "MANTOVANI" info_batch.csv  | awk -F ';' '{sum+=$5;} END{print "number of animals: " sum;}'



echo "#final report
finrep=NOMESNP

### SNP map (orinal from Illumina)
snpmap=NOMEMAP

### Often there are multiple allele codings in the row format files, chose the one you wish on your PED file
### Options allowed: 'top', 'forward','ab'.
allele='top'

### Position of the SNP ID in the file (usually is the first field, but may change)
SNPid_pos='2'

### Position of the INDIVIDUAL ID in the file (usually is the secodn field)
INDid_pos='1'

### Name of output PED and MAP files
outname=OUTPUT

### This will be used on the Fid column (first column in the PED)
brdcode='TEST'

# Options: ',' (for CSV) / ' ' (for TXT) / '\t' (for TSV)
sep=',' " > peddar.param



n_batch="$(awk 'END{print NR}' info_batch.csv )"
sequenza="$(seq -s " " 1 $n_batch)"
echo "" > merge_plink.txt

cd Scarico*

for it in $sequenza
  do
    i=$(awk -v nrow=$it 'NR==nrow' ../info_batch.csv) # la matricola del file in variabile
    fold=$(echo $i | awk -F ";" '{print $1}' )
    map=$(echo $i | awk -F ";" '{print $3}')
    ped=$(echo $i | awk -F ";" '{print $4}')
    echo "-----------------------"$fold"----------------------"
    cd $fold
    cp ../../peddar.param .
    sed -i "s/NOMESNP/$ped/g" peddar.param
    sed -i "s/NOMEMAP/$map/g" peddar.param
    sed -i "s/OUTPUT/$fold/g" peddar.param

    python3 $pedda_row_dir > ped.log
     # create map
    awk -F ',' '{print $3, $2, $1, $4, "0" }' $map  | sed '1d' > $fold'.map'
    echo  $(pwd)"/"$fold.ped  $(pwd)"/"$fold.map >> ../../merge_plink.txt
    cd ..
  done
cd ..

sed -i "1d" merge_plink.txt
wc -l merge_plink.txt
wc -l info_batch.csv # uni in piu perche ha l header

printf "###################################################################
 4.) FAI REPORT PER CONFTONTRLO CON I DATI ANARE e id matricola
##################################################################"

n_batch="$(awk 'END{print NR}' batch_to_unzip)"
sequenza="$(seq -s " " 1 $n_batch)"

echo "ID;MATR;BATCH;DENSITY;PANNEL" > all_ani_info.csv
cd Scarico_*

for it in $sequenza
do
  i=$(awk -v nrow=$it 'NR==nrow' ../batch_to_unzip | awk -F ";" '{print $1}')
  fold=$(echo "$i" | awk -F ";" '{print $1}')
  cd "$fold" || exit 1
  nsnp=$(awk 'END {print NR}' "$fold".map)
  genotype=$(grep ProjectName *_DNAReport.csv | awk -F ',' '{print $4}' | awk -F '=' '{print $2}' | awk -F '_' '{print $2$3}')
  gawk -F "," -v genotype="$genotype" -v batch="$fold" -v density="$nsnp" '{print $1";"$2";"batch";"density";"genotype}' *_RelMatr.csv
  cd ..
done >> ../all_ani_info.csv

cd ..
# REPLACE BIANCO CON TRATTINO BASSO

awk -F ";" '{gsub(" ","_",$2); print $1";"$2";"$3";"$4";"$5}' all_ani_info.csv > all_ani_info_trattbasso.csv
awk -F ';' '{print $5}' all_ani_info_trattbasso.csv | sort -k  1 | uniq -c

# fael