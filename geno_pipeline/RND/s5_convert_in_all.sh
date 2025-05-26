# converti con il blup e plink
cd  impute 
awk 'NR > 1 {print $1,$3}' genotypes_imp.txt >  impute.snp
awk 'NR > 1{print $2 ,$1  ,0 ,$3}' snp_info.txt  > impute.map
# converti in plink 

