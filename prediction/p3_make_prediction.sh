#!/bin/bash

echo "parte uno leggi i dati "

############################################################################
#######        predizione genomica vs normale whole data set      ##########
########               correggi il pedigree                   ##############
#############################################################################
#-----------------  calcola l'h2 del modello con tutti i dati -------------#

printf " step zero usa seek parents to correct ana"

mkdir seekp; cd seekp

ped="data/ana_corretta.txt"
snp="data/impute.snp" # quello completio

awk 'NR > 1{print $1,$2,$3,$4}' ../$ped > anap.txt

export OMP_NUM_THREADS=5
seekparentf90 --pedfile anap.txt --snpfile ../$snp --yob \
  --seeksire_in_ped --seekdam_in_ped --seektype 1 \
   --excl_thr_prob 1 --assign_thr_prob .5 \
     --duplicate  &> seek.log  &


awk '{print $1,$2,$3,$4}' Check_anap.txt > ../data/ana_corretta.txt 

cd ..


printf " setp 1 fai il modello completo\n"


mkdir complete_model 
cd complete_model

cp ../par.txt .
data="data/db.txt"
ped="data/ana_corretta.txt"
prune=5

sed -i "s:DATABASE:../$data:g" par.txt
sed -i "s:PEDIGREE:../$ped:g" par.txt
sed -i "s:PRUNE:$prune:g" par.txt

echo par.txt | renumf90 | tee log_renum
nohup ../../run_blup.sh  &

# make solution  #
python3 ../updata.py postmean
sed -i "s/OPTION method VCE/OPTION sol se/g" renf90.par
echo "OPTION origID " >> renf90.par
echo "OPTION store_accuracy 13" >> renf90.par

cd ..


printf " step 2 fai il modello ridotto\n"

mkdir model_ridotto
cd model_ridotto

cp ../par.txt .
data="data/db_with_test.txt"
ped="data/ana_corretta.txt"
prune=5

sed -i "s:DATABASE:../$data:g" par.txt
sed -i "s:PEDIGREE:../$ped:g" par.txt
sed -i "s:PRUNE:$prune:g" par.txt
echo "OPTION remove_all_missing" >> par.txt

echo par.txt | renumf90 | tee log_renum
nohup ../../run_gibbs.sh  &

printf "renf90.par\n20000\n100\n0" > post 
cat post | postgibbsf90

# remake the analysis with all info

sed -i "s:OPTION remove_all_missing::g" par.txt
echo par.txt | renumf90 | tee log_renum

# make solutions
echo postmean | python3 ../../updata.py 
sed -i "s/OPTION method VCE/OPTION sol se/g" renf90.par
echo "OPTION origID " >> renf90.par
echo "OPTION store_accuracy 11 orig" >> renf90.par

../../run_blup.sh

cd ..

# setp 3 (primo rsultato compara gli idniic)
# calcola l'h2 del modello ridotto

mkdir model_ridotto_genomico
cd model_ridotto_genomico


data="data/db_with_test.txt"
ped="data/ana_corretta.txt"
prune=5
snp="data/pruned.snp"
posg=11

cp ../par_g.txt par.txt

sed -i "s:DATABASE:../$data:g" par.txt
sed -i "s:PEDIGREE:../$ped:g" par.txt
sed -i "s:PRUNE:$prune:g" par.txt
sed -i "s:SNPPE:../$snp:g" par.txt

echo par.txt | renumf90 | tee log_renum

####################################
######  make solutions ssgblup #######
##############################à#########

echo  ../model_ridotto/postmean | python3 ../../updata.py
sed -i "s/OPTION method VCE/OPTION sol se/g" renf90.par
echo "OPTION origID " >> renf90.par
echo "OPTION store_accuracy 11 orig" >> renf90.par
../../run_blup.sh | tee blup.log

mv acc_bf90 acc_geno
mv solutions.orig solutions_geno


######################## io fate direttamente qua ##############################
###### perche se faccio solo pedgree perdo animali, i parenti di quelli genotipizti
#################################################################################

sed -i "s/OPTION SNP_file //g" renf90.par
../../run_blup.sh

mv acc_bf90 acc_no_geno
mv solutions.orig solutions_no_geno

awk '{print $3,$4,$5}' acc_geno    | sort -k 1b,1 > g
awk '{print $3,$4,$5}' acc_no_geno | sort -k 1b,1 > nog

join -1 1 -2 1 g nog > compa.txt 

awk '{print $1,$4,$5}' ../$ped | sort -k 1b,1 > p

join -1 1 -2 1 compa.txt p > compa_with_info.txt

# wirth geno 
join -1 1 -2 1 <(awk '{print $1}' ../$snp | sort -k 1 ) <( sort  -k 1 compa_with_info.txt)  > comp_geno.txt


cp ../par_g.txt par.txt

sed -i "s:DATABASE:../$data:g" par.txt
sed -i "s:PEDIGREE:../$ped:g" par.txt
sed -i "s:PRUNE:$prune:g" par.txt
sed -i "s:SNPPE:../$snp:g" par.txt

echo par.txt | renumf90 | tee log_renum

####################################
######  make solutions ssgblup #######
##############################à#########

echo  ../complete_model/postmean | python3 ../../updata.py
sed -i "s/OPTION method VCE/OPTION sol se/g" renf90.par
echo "OPTION origID " >> renf90.par
echo "OPTION store_accuracy "$posg" orig" >> renf90.par
../../run_blup.sh | tee blup.log

mv acc_bf90 acc_geno
mv solutions.orig solutions_geno


######################## io fate direttamente qua ##############################
###### perche se faccio solo pedgree perdo animali, i parenti di quelli genotipizti
#################################################################################

sed -i "s/OPTION SNP_file //g" renf90.par
../../run_blup.sh

mv acc_bf90 acc_no_geno
mv solutions.orig solutions_no_geno

awk '{print $3,$4,$5}' acc_geno    | sort -k 1b,1 > g
awk '{print $3,$4,$5}' acc_no_geno | sort -k 1b,1 > nog

join -1 1 -2 1 g nog > compa.txt 

awk '{print $1,$4,$5}' ../$ped | sort -k 1b,1 > p

join -1 1 -2 1 compa.txt p > compa_with_info.txt

# wirth geno 
join -1 1 -2 1 <(awk '{print $1}' ../$snp | sort -k 1 ) <( sort  -k 1 compa_with_info.txt)  > comp_geno.txt

cd ..


#######################################################
#            PARTE DUE FAI IL MODELLO LR               #
#######################################################


mkdir model_ridotto_genomico_lr
cd model_ridotto_genomico_lr

# qua rimuvovo i data test perche fanno casino 
# perch devo prevedere i fenotipi #

data="data/db_with_test.txt"
ped="data/ana_corretta.txt"
prune=5
snp="data/pruned.snp"
posg=11
column=10
da_che_anno=2018
a_che_anno=2018

echo "remove young animals"
awk -v nc=$column '$nc!=0{print $0}'  ../$data  > db_no_test.txt
#################################################

echo "remove young animals also form the chips"
join -1 1 -2 1 <(sort -k  1 ../data/train.txt) <(sort -k 1 ../$snp) > snp_test.txt


echo "making the prediction in the whole"

cp ../par_g.txt par.txt
sed -i "s:DATABASE:db_no_test.txt:g" par.txt
sed -i "s:PEDIGREE:../$ped:g" par.txt
sed -i "s:PRUNE:$prune:g" par.txt
sed -i "s:SNPPE:snp_test.txt:g" par.txt

echo "makin whole ssGblup"

echo par.txt | renumf90 | tee log_renum

echo ../model_ridotto/postmean | python3 ../updata.py 

sed -i "s/OPTION method VCE/OPTION sol se/g" renf90.par
echo "OPTION origID " >> renf90.par
echo "OPTION store_accuracy "$posg" orig" >> renf90.par

../../run_blup.sh

####### qua potrei fare il gwas ######

mv acc_bf90 acc_w_geno
mv solutions.orig solutions_w_geno

echo "makin whole Pblup"


sed -i "s/OPTION SNP_file //g" renf90.par
../../run_blup.sh
mv acc_bf90 acc_w_nogeno
mv solutions.orig solutions_w_nogeno

echo "crea il dataset con tutti gli anni"

awk 'NR > 1 {print $1}' db_no_test.txt | sort -k 1 | uniq > anim.txt
wc -l anim.txt
awk '{print $1,$4}' ../$ped | sort -k 1 > tmp
join -1 1 -2 1 anim.txt tmp > to_lr.txt

awk '{print $2}' to_lr.txt | sort -k 1 -n  | uniq -c  > distribution_by_year.txt


#awk -v field="$column" 'NR==FNR{arr[$1]; next} $1 in arr{$field=0} 1' test.txt ../db_no_test.txt | less -S



while [ $da_che_anno -lt $a_che_anno ]; do
    echo "###############################################################"
    echo "estraggo gli anni"
    echo $da_che_anno
    
    mkdir "dir_"$da_che_anno
    cd "dir_"$da_che_anno    
    
    awk -v f=$da_che_anno '$2 >= f {print $0}' ../to_lr.txt | awk '{print $1}' >  test.txt
    
    echo "number of test" $(wc -l test.txt)

      # Replace with the column number you want to set to zero

    awk -v field="$column" 'NR==FNR { arr[$1]; next } 
                    $1 in arr { $field = 0 } 
                        { print }' test.txt ../db_no_test.txt > db.txt

    cp ../../par_g.txt par.txt
    sed -i "s:DATABASE:db.txt:g" par.txt
    sed -i "s:PEDIGREE:../../$ped:g" par.txt
    sed -i "s:PRUNE:$prune:g" par.txt
    sed -i "s:SNPPE:../snp_test.txt:g" par.txt

    echo $(awk C=$column '$C == 0' db.txt | wc -l)

    echo "makin ssGblup"

    echo par.txt | renumf90 > log_renum 

    echo ../../model_ridotto/postmean | python3 ../../../updata.py 

    sed -i "s/OPTION method VCE/OPTION sol se/g" renf90.par
    echo "OPTION origID " >> renf90.par
    echo "OPTION store_accuracy "$posg" orig" >> renf90.par

    ../../../run_blup.sh > loggblup 
    mv acc_bf90 acc_p_geno
    mv solutions.orig solutions_p_geno


    echo "now makin pblup"

    sed -i "s/OPTION SNP_file //g" renf90.par
    ../../../run_blup.sh > lopblup 

    mv acc_bf90 acc_p_nogeno
    mv solutions.orig solutions_p_nogeno

cd ..

    ((da_che_anno++))
done  > log_lr.txt &







