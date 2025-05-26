cp ../../GENOTIPI/GENOTIPI_RENDENA/DATA/data_ready_after_impute/imputed_f901.txt data/impute.snp


# setp 1

mkdir complete_model 
cd complete_model

cp ../par.txt .
data="data/db.txt"
ped="data/ana.txt"
prune=5

sed -i "s:DATABASE:../$data:g" par.txt
sed -i "s:PEDIGREE:../$ped:g" par.txt
sed -i "s:PRUNE:$prune:g" par.txt

echo par.txt | renumf90 | tee log_renum
nohup ../../run_gibbs.sh &

printf "renf90.par\n20000\n100\n0" > post 

cat post | postgibbsf90

cd ..


