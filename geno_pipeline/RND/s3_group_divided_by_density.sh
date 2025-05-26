##############################################################
# ------------- sostiruirlo con quello anaborare ------------#
#############################################################
#07/05/2025
####### script da sistemare #######

#--------------------------------------------------------#
#  PARTE 5 CREA UN UNICO PANELLO HD E MD 
#  quando arrivano i 100k trovare una strategia
#  fai per numero di chips perche nomi è un casino 
#-------------------------------------------------------#

awk -F ";"  '$6="type" {}1'  all_ani_info_trattbasso.csv > tmp
mv tmp all_ani_info_trattbasso.csv

awk -F ";" '{
    if ($5 ~ /LD/)   { $6 = "33K" }
    if ($5 ~ /HD/)   { $6 = "150K" }
    if ($5 ~ /LGSB/) { $6 = "60K" }
    print $0
}' all_ani_info_trattbasso.csv > tmp

mv  tmp all_ani_info_trattbasso.csv


# PARTE UNO 
# dividi per pannel !
####

mkdir pannel_2025
cd pannel_2025 

printf "CONTROLLA SEMPRE A MANO PERCHE QUELLI DI AGROTIS
FANNO CASINO ALCUNI 150K LI CHIAMANO LD\n" 

# step 1 crea file di rifermeinto

awk -F " " 'NR>1{print $6}' ../all_ani_info_trattbasso.csv | sort -k  1 | uniq | sed -e '1d' > chip.txt

# io ho modificato questo intanto
awk '{if($3=="19CA01173") {$6="150K"} print $0}' ../all_ani_info_trattbasso.csv > tmp
rm ../all_ani_info_trattbasso.csv
mv tmp ../all_ani_info_trattbasso.csv

while IFS= read -r i; do
    echo "#########################"$i"##################################"
    awk -v chips="$i" '$6 == chips {print $3, $4}' ../all_ani_info_trattbasso.csv | sort -u |
    while IFS= read -r ix; do
        echo "$ix"
    done
done < chip.txt

#poi copiali dal file e metteri nella cartele per plink
where=$(ls -d ../Scarico_*)
while IFS= read -r i; do
    echo "$i"
    awk -v dir="$where" -v chips="$i" '{if($6==chips) print dir"/"$3"/"$3".ped", dir"/"$3"/"$3".map"}' ../all_ani_info_trattbasso.csv | sort -k 1 | uniq > "${i}_merged.txt"
done < chip.txt


# step 2 crea file per densità

while IFS= read -r i; do
     echo "@@@@@@@@@@@@@@@@@@@@@@@ $i @@@@@@@@@@@@@@@@@@@@@@@@@"
     plink --cow --allow-extra-chr  --merge-list "${i}_merged.txt" --recode --out "${i}_merged"
     plink --cow --allow-extra-chr  --file "${i}_merged" --missing --out "${i}_merged_call_rate"
done < chip.txt

wc -l *.ped

echo "" > merge_all.txt
while IFS= read -r i; do
    echo "$i"
    awk '{print $2}' $i"_merged.ped" >> merge_all.txt
done < chip.txt

sed -i '1d'  merge_all.txt

printf "controlla queli che ho perso\n"

awk '{print $1}' ../all_ani_info_trattbasso.csv > tmp
grep -vxFf tmp merge_all.txt 
 
echo "NB hanno due numeri diversi perche ci sono duplicati.-.."

awk '{count[$1]++; lines[$1] = lines[$1] $0 RS} END {for (key in count) if (count[key] > 1) printf "%s", lines[key]}' \
   ../all_ani_info_trattbasso.csv > ../id_duplicati.txt

cd ..



#
# SECONDA PARTE DOPO AVERE DIVISO PER DENSITA RECUPERO LA RAZZA 
#


cd Scarico_*/
rm ../matr_id.txt

awk -F ";" '{print $1}' ../info_batch_all.csv | while IFS= read -r i; do
    if [[ -n $(ls "$i/"*.map 2>/dev/null) ]]; then
        echo "$i"
        if [[ -n $(ls "$i/"*_BGC.txt 2>/dev/null) ]]; then
            awk -v c="$i" -F " " '{for (i=1; i<=NF; i++) if ($i ~ /\.ZIP/) print $i, $3, c, "BGC"}' "$i/"*_BGC.txt >> ../matr_id.txt
        elif [[ -n $(ls "$i/"*_BL3.txt 2>/dev/null) ]]; then
            awk -v c="$i" -F " " '{for (i=1; i<=NF; i++) if ($i ~ /\.ZIP/) print $i, $3, c, "BL3"}' "$i/"*_BL3.txt >> ../matr_id.txt
        elif [[ -n $(ls "$i/"*_B55.txt 2>/dev/null) ]]; then
            awk -v c="$i" -F " " '{for (i=1; i<=NF; i++) if ($i ~ /\.ZIP/) print $i, $3, c, "B55"}' "$i/"*_B55.txt >> ../matr_id.txt
        else
            printf "No matching *_BGC.txt, *_BL3.txt, or *_B55.txt file found in $i/ assignin rendena\n"
            awk -v c="$i" -F "," '{print $1,"10",c,"None"}' "$i/"*RelMatr* >> ../matr_id.txt
            
        fi
    fi
done

cd ..

# ok asesso perche sia definiitpo estraggo l'id id laboratori che e la riga 1
awk -F " " '{new_col1 = substr($1, 1, 12); $1 = new_col1; print}'   matr_id.txt > id_breed.txt

# perche cìè questa disparita ?? capire questo e capire
# dove sono gli animali che mancano
wc -l  matr_id.txt
wc -l  all_ani_info_trattbasso.csv

###### FINE RECUPERO RAZZE

awk '$2==10 {print $3}' id_breed.txt | sort | uniq -c
awk '$2==10 {print $3}' id_breed.txt | sort | uniq -c | wc -l


#
#     PARTE 3 DIVIDO PER DENSITA E RAZZA
#


cd pannel_2025

awk '{print $2}' ../id_breed.txt | sort | uniq > breed 

printf "faccio solo i 150 perche tutti gli atri rendeni
        in piu conserva i dati per il progetto dual breeding da invaiti\n"

while IFS= read -r i; do
    echo $i
    awk -v breed="$i" '$2 == breed {print "TEST", $1}' ../id_breed.txt | sort | uniq > "breed_${i}.txt"
    echo $(wc -l "breed_${i}.txt")
    #plink --cow --bfile 150K_merged --keep "breed_${i}.txt" --recode 12 --out "breed_${i}"
    plink --cow --file 150K_merged --keep "breed_${i}.txt" --make-bed --out "breed_${i}"
    awk '{ $2 = "a" NR } 1'  "breed_${i}".ped > "breed_${i}"_sas.ped
    cp  "breed_${i}".map  "breed_${i}"_sas.map
done < breed

printf " controlla che i numeri corrispondano"


#------------
#  Parte 4 
#  CONVERTI I FILE CON I VALORI REALI E ELIMINA MATRICOLE IN BASE ALL CALL RATE
#--------------

mkdir prepare_input; cd prepare_input
cp ../breed_10.bed    ../breed_10.bim  ../breed_10.fam . 
cp ../33K_merged.bed ../33K_merged.bim ../33K_merged.fam .
cp ../60K_merged.bed ../60K_merged.bim ../60K_merged.fam .

echo "adesso per ogni file metti la matricola con il trattino basso e coreggi per duplicati"


# Read each line from the file and process
for f in $(ls *.fam)
do
   i="${f%.fam}" # rimuovi fam
  echo "############# $i ##########"
  # Count the lines in the .ped file
  wc -l $i.fam
  echo "totale: $(wc -l < $i.fam)" > n_$(basename $i).txt
  
  # Remove duplicate matriculation numbers
  awk -F ";" '{print $1,$2}' ../../all_ani_info_trattbasso.csv | awk '{print $1,$2}' | sort -k1 > sort.ani
  uniq sort.ani > sort.ani.tmp
  mv sort.ani.tmp sort.ani
  
  plink --cow --bfile $i --missing
  
  awk '{print $2, 1 - $6}' plink.imiss | sort -k1 > sort.miss
  join -1 1 -2 1 sort.miss sort.ani > ani_matr_id.txt
  
  awk '{seen[$3]++; lines[$3] = lines[$3] $0 ORS} END {for (value in seen) if (seen[value] > 1) printf lines[value]}' ani_matr_id.txt > dpl.txt
  sort -k2,2r -k3,3n dpl.txt | sort -k3 -u > massimo_call_rate.txt
  sort dpl.txt | uniq > uguali.txt
  grep -vxFf massimo_call_rate.txt dpl.txt > massimo_errore.txt
  awk '{print "TEST", $1}' massimo_errore.txt > cava0.tmp
  
  cat cava0.tmp > cava.txt
  
  echo "ne perdo: $(wc -l < cava.txt)"
  
  plink --cow --bfile $i --make-bed --remove cava.txt --out $(basename $i)
  plink --cow --bfile $(basename $i) --recode 12 --out $(basename $i)
  
  i=$(basename $i) # for security

  sort -k3 ani_matr_id.txt | uniq -d
  sort -k1,1b $i.ped > sort.ped
  wc -l $i.ped
  join -1 1 -2 2 sort.ani sort.ped > mm.ped.tmp
  wc -l mm.ped.tmp
  awk '{t = $2; $2 = $3; $3 = t; print}' mm.ped.tmp > ok.tmp
  cut -d " " -f2- ok.tmp > $i.ped
  rm *.tmp sort.*
  
  echo "numero dopo la pulizia: $(wc -l < $i.ped)" >> n_$i.txt
  awk '{seen[$2]++; lines[$2] = lines[$2] $0 ORS} END {for (value in seen) if (seen[value] > 1) printf lines[value]}' $(basename $i).ped | wc -l
  
  sort -k1 $i.ped | uniq > mm.ped.tmp
  mv mm.ped.tmp $i"_matr.ped"
  mv $i.map $i"_matr.map"
  
done 

echo "removing basura "
rm  plink* n_breed* plink* uguali.txt 
rm cava*  massimo_call_rate*  

echo "finio"