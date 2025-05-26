
echo  "creo una cartella unica perche anche i 100k"
cd /mnt/user_data/enrico/Genotipi/Rendena/combined

cp ../../Rendena/Neogen_100k/organized_*/neogenPlink.map ../../Rendena/Neogen_100k/organized_*/neogenPlink.ped .
cp  /mnt/user_data/enrico/Genotipi/Rendena/DUALBR/pannel_2025/prepare_input/*.ped  .
cp  /mnt/user_data/enrico/Genotipi/Rendena/DUALBR/pannel_2025/prepare_input/*.map  .


# dato che non lo faccio rimuvo doppioni nel 100k
awk '!seen[$2]++'  neogenPlink.ped > neogenPlink2.ped
mv neogenPlink2.ped  neogenPlink.ped

printf "########################################\n
step 2 rimuovi aniamli in tra i chips \n
la logica che se ho un animale in comune tra due chips prendo quello \n
sul chips piu denso oppure il pannelo da downgradare come pannello piugrande \n
scelgo io l'orodine dei panneli in base all numero di righe ne l bim\n\n"


echo "mettelri in ordine di densitaaaaaaaaaaa"
ls -1 *map | xargs -I {} sh -c 'echo "$(wc -l {})"' | sort -n | awk '{print $2}' > chefile.txt
sed -i 's/.map//g' chefile.txt

# anti grep e poi lo metto sompre cosi 
# i pannello da dowgrade diventa il flucro
#if [ "$todown" != "no" ]; then
#    grep -v "$todown" chefile.txt > tmp
#    echo "${todown}_cln_uniq.bim" >> tmp 
#    mv tmp chefile.txt 
#fi

# greppa solo quelli che mi interssano
#mv pannel_keep.txt chefile.txt

#echo -e ""breed_10_ld33k_merged_matr"\nbreed_10_matr" > chefile.txt
echo " cosa piu conservativa converti in plink bed ofmrat e rimuvoi  anche i cormosomi sesuali"

while IFS= read -r info; do
    echo $info
    plink --cow --allow-extra-chr --chr 1-29 --file $info --make-bed --out $info 
done < chefile.txt


####################################################
#   IN QUESTA PAERTE RIMUOVO I GLI ANIMALI
##################################################

row=2
totr=$(awk 'END{print NR}' chefile.txt)

while IFS= read -r file; do
    if [[ "$row" -gt "$totr" ]]; then
        break
    fi
    row1=$row
    echo "#########"
    echo $file 
    while true; do
        file1=$(echo $file".fam")    
        tmp=$(awk -v r="$row1" 'NR==r{print $1}' chefile.txt)
        file2=$(echo $tmp".fam")
        echo $file1 "vs" $file2
        join -1 2 -2 2 <(sort -k 2 $file1 ) <(sort -k 2 $file2) | awk '{print $2,$1}' > animal.txt 
        filex=$file
        if [[ -s "animal.txt" ]]; then
            echo "has n animal in common " $(wc -l  animal.txt | awk '{print $1}') "remove on chips" $file1
           plink --cow --bfile $filex --make-bed --remove animal.txt --out $filex  > log
        fi
            echo "No animals in common"
            plink --cow --bfile $filex --make-bed --out $filex > log
        ((row1++))
        if ((row1 > totr)); then
            break
        fi
        # Uncomment the line below to execute plink command
        #
    done
    ((row++))
done < chefile.txt | tee duplictae_between_chips_log.txt

rm *~ 

#--------------------------#
#-   Aggiorno le mappe    -# 
#--------------------------#
# copia e rimuvoi se hanno meno di due colonne xk vuold ire che non c'è il chr
map100k=../../mappe_reference/GGP_Bov_100K_HTS_20040701_A1_2.csv
map150k=../../mappe_reference/geneseek-ggp-bovine-150k-manifest-file.csv
map33k=../../mappe_reference/GGP_Bovine_LD_v3_v2.0.csv
map60k=../../mappe_reference/BovineSNP50_v3_A2.csv

#

update_plink_chr_map() {
    #echo $(ls $plinkfile*)
  local refmap="$1"
  local plinkfile="$2"

  awk -F ',' 'NR > 8 {print $2, $10, $11}' "$refmap" | awk 'NF >= 3' > ref.tmp
  
  plink --cow --bfile "$plinkfile" \
        --update-chr ref.tmp 2 1 \
        --update-map ref.tmp 3 1 \
        --make-bed \
        --out "${plinkfile}_updt" | tee "log_upd_"$2
}


update_plink_chr_map $map33k "33K_merged_matr"
update_plink_chr_map $map150k "breed_10_matr"
update_plink_chr_map $map100k "neogenPlink"
update_plink_chr_map $map60k  "60K_merged_matr"



wc -l *updt.fam 


#---------------------------------------
# anche per i nomi degli snp se hanno nomi 
# intanto faccio questo downgrade
# appena so le mappe update map
#---------------------------------------
todown="breed_10_matr.map"


if [ "$todown" != "no" ]; then
echo "si cè una mappa reference ed è: "$todown
        while IFS= read -r file
            do  
                        printf "@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                        
                        echo  "$file" "=>" "$todown"
                        
                        ed=$(echo $file | sed -e 's/.map//g')
                        ed="${ed}_updt.bim"
                        
                        todwn=$(echo $todown | sed -e 's/.map//g')
                        todwn="${todwn}_updt.bim"
                        
                        echo  "$ed" "=>" "$todwn"

                        # faccio il merge per cromosoma e poszione
                        join -1 1 -2 1 <( awk '{print $1"_"$4}' $todwn | sort -k 1 ) \
                                        <(awk '{print $1"_"$4}' $ed | sort -k 1  ) > incommon.txt
        
                         # estri i nmi dalla mappa
                        join -1 1 -2 1 <(awk '{print $1"_"$4,$0}' $ed | sort -k 1b,1 )  <(sort -k 1b,1  incommon.txt) | awk '{print $3}' > keep.txt
  
                        ed2=$(echo $ed | sed -e 's/.bim//g')
                        
                        plink2 --cow --bfile $ed2 --extract keep.txt --make-bed --out $ed2"_rec"  > log
                        
                        
                        tooo=$(awk 'END{print NR}' $ed)
                        frommm=$(awk 'END{print NR}' $ed2"_rec.bim" )
                        echo "from" $tooo "to" $frommm 
                        

        done < chefile.txt
fi


#sed -i 's/60K_merged_matr//g' chefile.txt
sed -i '/^$/d' chefile.txt
echo "create file snp.."

awk '{print $0"_updt_rec"}' chefile.txt > tmp
mv tmp chefile.txt

# rimuvoi snp doppio e in comume

while IFS= read -r i
 do   

    filex=$i
    
    echo "###############################"
    echo $filex
    echo "removing duplicate names..."
    echo $(ls $filex*)
    
    awk '$5 != $6' $filex.bim > good_snps.bim
    
    printf "vedo quanti sono codificati bene \n"
    echo $(wc -l $filex.bim)
    echo $(wc -l good_snps.bim)

    cut -f2 good_snps.bim > keep_snps.txt
    plink2 --cow --bfile  $filex --extract keep_snps.txt --rm-dup \
                        --make-bed --out tmp | tee log_dup
                        
    echo "removing differen snp with same position and same chr.."

    awk '{print $1"_"$4,$2}' "tmp.bim"  > tmp
    awk '{print $1}' tmp | sort -k 1b,1 | uniq -c | awk '$1>1 {print $0}' > dpl.txt
    
    if [[ -s dpl.txt ]]; then
        join -1 2 -2 1 <(sort -k2b,2 dpl.txt) <(sort -k1b,1 tmp) | awk '{print $3}' > rm.tmp
        plink2 --cow --exclude rm.tmp --bfile tmp --make-bed --out "${filex}_ndup" > log
        echo "$(wc -l < dpl.txt) duplicates found."
    else
        echo "No duplicates found — skipping removal step. rinomino lo stesso"
        cp $filex".bed" $filex"_ndup.bed"
        cp $filex".fam" $filex"_ndup.fam"
        cp $filex".bim" $filex"_ndup.bim"
        
        
    fi

done < chefile.txt

 wc -l *_ndup.fam
 
#--- fino q qua 
# ordina e crea il file 
#------------------------------

awk '{print $0"_ndup"}' chefile.txt > tmp
mv tmp chefile.txt


echo "ID chips genotype" > file.snp
row=1
while IFS= read -r i
 do   
    echo "#########"
    filex=$(echo $i | sed -e 's/.map//g')
    filex=$filex"_updt_rec_ndup"
    
    echo "sorting.."
    plink2 --cow --bfile $filex --sort-vars  --make-pgen --out  $filex   > logd
    echo "make raw and pedmap .."
    plink2 --cow --pfile $filex --recode ped 12  --out  $filex  > logd
    plink2 --cow --pfile $filex --recode A --out  $filex   > logd
    echo "adding on file.."
    awk 'NR > 1 { printf "%s%s", $2, OFS; 
        for (i = 7; i <= NF; i++) printf "%s", $i; print "" }' $filex".raw" > tmp
    awk -v IS=$row '{ gsub("NA", "5", $2); print $1,IS,$2 }'  tmp >> file.snp
    ((row++))
done < chefile.txt




while IFS= read -r i
 do   
    filex=$i #"_cln_uniq_rec"
    echo $filex
    awk 'NR==FNR { file1[$2] = FNR; next } 
            $2 in file1 { print $0, file1[$2]; next } { print $0, 0 }' $filex.map s.map >  tmp
    mv tmp s.map
done < chefile.txt


# SNP ID, chromosome number, base pair position, order of SNP for each chip

awk '{printf "%s %s %s ", $2, $1, $4; for (i=5; i<=NF; i++) printf "%s ", $i; print ""}' s.map  > file.map
chip=$(awk 'NR==1 {for (i=5; i<=NF; i++) printf "chip_%d ", i-4; print ""}' s.map)
sed  -i "1i SNP_ID chr pos $chip"   file.map







cp   /home/enrico/aspa2025/data/ana/ana.txt .
awk '{print $1,$2,$3,$5,$4}' ana.txt > ped.txt


printf " es un par sencillo \n despues lo modico"

cat > par.txt  <<EOF
title="Imputation of Project_name";
genotype_file="file.snp";
snp_info_file="file.map";
output_folder="imputazione";
ped_file="ped.txt";
keep_og;
save_genotype; 
parentage_test  ;
ref_chip=4; 
njob=5;
EOF

echo par.txt | /home/enrico/FImpute3 -o | tee log_resul

# vai dentro

awk 'NR > 1 {print $1,$3}' genotypes_imp.txt >  impute_.snp
awk 'NR > 1{print $2 ,$1  ,0 ,$3}' snp_info.txt  > impute.map

#plink --cow --file impute --make-bed --out impute

cp ../ana.txt ana.txt



awk 'length($1) == 16'  impute_.snp > impute.snp

/home/enrico/blupf90_old/seekparentf90 --thr_call_rate .0001 --maxsnp 200000 \
                 --excl_thr_prob 1 \
                 --seeksire_in_ped \
                 --seekdam_in_ped \
                 --seektype 1 \
                 --yob \
                --pedfile ana.txt --snpfile  impute.snp     \
                --assign_thr_prob 0.5 >  seek.log &


# covert in plink 







# now comvere

# crea il file da imputare 
# file snp...

echo "" > merge.txt
while IFS= read -r i
 do   
    echo "#########"
    filex=$i
    plink --cow --bfile $filex --recode 12 --out $filex
    echo $filex".ped" $filex".map"  >> merge.txt
done <  chefile.txt

sed -i '1d' merge.txt
plink --cow --merge-list  merge.txt --recode ped 12 --out s > log_merge

echo "tranforma in blupf90 , questo serve per fare un check.." 
plink --cow --file s --recode A --out s > log_makraw

awk '{printf "%s %s %s ", $2, $1, $4; for (i=5; i<=NF; i++) printf "%s ", $i; print ""}' s.map  > file.map
chip=$(awk 'NR==1 {for (i=5; i<=NF; i++) printf "chip_%d ", i-4; print ""}' s.map)
sed  -i "1i SNP_ID chr pos $chip"   file.map










printf  "TO PLINK TO BLUPF90 \n"

awk 'NR > 1{ 
        printf "%s ", $2; for (i=7; i<=NF; i++) printf "%s", $i;
             printf "\n" }' s.raw | awk '{gsub("NA",5,$2); print $0}'  > rawf90.snp

echo "SNP_ID CHR POS" >  f90.map
awk -F'\t' '{print $2, $1, $4}' s.map >> f90.map

#                       @@@@@@@@@@@@@@@@@@@@
awk 'length($1) == 16' rawf90.snp > rawf901.snp
awk 'length($1) > 16' rawf90.snp > duplicated.snp

cp   /home/enrico/aspa2025/data/ana/ana.txt .
awk '{print $1,$2,$3,$4,$5}' ana.txt  > ana1.txt


/home/enrico/blupf90_old/seekparentf90 \
                --thr_call_rate 0.001 \
                --maxsnp 200000 \
                 --excl_thr_prob 1 \
                --pedfile ana1.txt --snpfile  rawf901.snp | tee  seek_preimpute.log #&




awk '{print $2}'  rawf901.snp | sed -e 's/5/9/g' | awk '{gsub(/[^[:space:]]/, "& "); print}' > tmp

paste -d ' ' <(awk '{print $1 }'  rawf901.snp ) tmp > ai.snp


join -1 1 -2 1 <( awk '{print $1}' ai.snp | sort -k 1  ) <(sort -k 1 ana1.txt )| awk '{print $1,$2,$3}' > ped.txt


source ../../code_and_sftw/venv/bin/activate


if como?=="whole"
python3 ../../code_and_sftw/alphaimpute2/alphaimpute2.py  -genotypes ai.snp  \
                    -pedigree ped.txt \
                    -min_chip 0.3 \
                    -iothreads 1 \
                    -maxthreads 1 \
                    -cycles 1 \
                    -onlykeyed \
                    -final_peeling_threshold 0.99 \
                    -hd_threshold .80 \
                    -out impute_whole  | tee impute.log


else if


#!/bin/bash

# Read the map file to get chromosome lengths
# cre
echo "0" >  N_SNPxCHR.txt
awk -F' ' 'NR>1{print $2}' f90.map | sort -k 1 | uniq -c | awk -F' ' '{print $1}' >> N_SNPxCHR.txt

awk '{ sum += $1 
print sum}' N_SNPxCHR.txt | sed -e '1d' > b.txt

echo "1" >  a.txt
 awk '{print $1=$1+1}' b.txt  >> a.txt 
sed -i '$ d' a.txt

paste a.txt b.txt  # just sto ckeck

awk '{ sum += $1 
print sum}' N_SNPxCHR.txt > a.txt


echo "Start chromosome-by-chromosome imputation"

# Create a directory for imputed chromosomes
mkdir -p impute_by_chr

for ((i=1; i<30; i++)); do
   
    echo "----------------------------------------"
    echo "---------------- Chromosome_$i ----------------"
    
    # Set the start and stop SNPs
    from=$(awk -v r="$i" 'NR==r{print $1}'  a.txt)
    to=$(awk -v r="$i" 'NR==r{print $1}'  b.txt)
    

    
    echo "From" $from "'SNP <-----> To " $to "'s SNP.."

    # Run AlphaImpute2 for each chromosome
    # ulimit -s unlimited
    python3 ../../code_and_sftw/alphaimpute2/alphaimpute2.py \
        -genotypes ai.snp \
        -pedigree ped.txt \
        -min_chip 0.3 \
        -iothreads 2 \
        -maxthreads 2 \
        -cycles 10 \
        -onlykeyed \
        -final_peeling_threshold 0.99 \
        -hd_threshold 0.80 \
        -startsnp "$from" -stopsnp "$to" \
        -phasing_loci_inclusion_threshold 0.9 \
        -out "impute_by_chr/IMPchr_$i"
done

fi


# unisci tutto ora

cd impute_by_chr/
awk '{print $1 }' ../ped.txt > tieni_solo_questi.txt
# creo un unico pannelo sh
chr=$(seq 29)
for i in $chr ;do
    sort -k 1 tieni_solo_questi.txt > Z1
    sort -k 1 IMPchr_$i.genotypes > Z2
    join -1 1 -2 1 Z1  Z2 > chr_$i
    rm Z1 Z2
    echo $i
done
echo "merge finito"


printf "@@@@@@@@@@@@"


join -1 1 -2 1 <(sort -k 1 file  ) <(awk '{print $1}' ../ped.txt | sort -k 1 )  > imp.txt

awk '{ printf "%s ", $1; for (i=2; i<=NF; i++) printf "%s", $i;
             printf "\n" }' imp.txt | awk '{gsub(9,5,$2); print $0}' >  imputeaif90.snp


cp ../f90.map imputeaif90.map


printf "check"



awk 'NR > 1{print $1,$2,$3,$4,$5}' ../ana1.txt > ana.txt

/home/enrico/blupf90_old/seekparentf90 --pedfile ana.txt --snpfile   imputeaif90.snp | tee seek.log
                #--thr_call_rate 0.1 
            #    --maxsnp 200000 \
                # --excl_thr_prob 1 \
                # --seeksire_in_ped \
                # --seekdam_in_ped \
                # --seektype 1 \
#--yob \
                --pedfile ana.txt --snpfile   imputeaif90.snp     \
                --assign_thr_prob 0.5 | tee seek.log #&



# fare file da qua 


while IFS= read -r i
 do   
    echo "#########"
    filex=$(echo $i | sed -e 's/.map//g')
    filex=$filex"_updt_rec_ndup"
    echo $filex
    awk 'NR==FNR { file1[$2] = FNR; next } 
            $2 in file1 { print $0, file1[$2]; next } { print $0, 0 }' $filex.map s.map >  tmp
    mv tmp s.map
done < chefile.txt


# SNP ID, chromosome number, base pair position, order of SNP for each chip



# metterlo nel par magari

FImpute3 par.txt -o


mv file.snp "/"$title
mv file.map "/"$title





# cambia il parmtro per il file machrt
# convert to usable format

cd  impute 
awk 'NR > 1 {print $1,$3}' genotypes_imp.txt >  impute.snp

awk 'NR > 1{print $2 ,$1  ,0 ,$3}' snp_info.txt  > impute.map

#plink --cow --file impute --make-bed --out impute

cp ../ana.txt ana.txt

awk 'NR > 1{print $1,$2,$3,$4,$5}' ../../../../ana/ana_all.txt > ana.txt
nohup  seekparentf90 --thr_call_rate .0001 --maxsnp 200000 \
                 --excl_thr_prob 1 \
                 --seeksire_in_ped \
                 --seekdam_in_ped \
                 --seektype 1 \
                 --yob \
                --pedfile ana.txt --snpfile  impute.snp     \
                --assign_thr_prob 0.5 >  seek.log &



# scirpt per sostituisco
# cpaoire vpnr ho volari diversi a seconda

awk '{print $2}'  rawf90.snp | sed -e 's/5/9/g' | awk '{gsub(/[^[:space:]]/, "& "); print}' > tmp

paste -d ' ' <(awk '{print $1 }'  rawf90.snp ) tmp > ai.snp

join -1 1 -2 1 <( awk '{print $1}' ai.snp | sort -k 1  ) <(sort -k 1 ../../../ana/ana_all.txt )| awk '{print $1,$2,$3}' > ped.txt


source ../../../venv/bin/activate

nohup python3 ../../../alphaimpute2/alphaimpute2.py  -genotypes ai.snp  \
                    -pedigree ped.txt \
                    -min_chip 0.3 \
                    -iothreads 20 \
                    -maxthreads 20 \
                    -cycles 5 \
                    -onlykeyed \
                    -final_peeling_threshold 0.99 \
                    -hd_threshold .80 \
                    -out impute_whole  &

join -1 1 -2 1 <(sort -k 1 impute_whole.genotypes  ) <(awk '{print $1}' ped.txt | sort -k 1 )  > imp.txt

awk '{ printf "%s ", $1; for (i=2; i<=NF; i++) printf "%s", $i;
             printf "\n" }' imp.txt | awk '{gsub(9,5,$2); print $0}' >  imputeaif90.snp


awk 'NR > 1{print $1,$2,$3,$4,$5}' ../../../ana/ana_all.txt > ana.txt

seekparentf90 --thr_call_rate 0.1 --maxsnp 200000 \
                 --excl_thr_prob 1 \
                 --seeksire_in_ped \
                 --seekdam_in_ped \
                 --seektype 1 \
                 --yob \
                --pedfile ana.txt --snpfile   imputeaif90.snp     \
                --assign_thr_prob 0.5 | tee seek.log #&




