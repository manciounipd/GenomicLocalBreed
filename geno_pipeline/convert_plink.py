import os
import pandas as pd

def altre_info(value1, genotype):
    sex = "0"
    genotype_str = ' '.join(genotype)
    return f"REND {value1} 0 0 {sex} -9 {genotype_str}\n"

def convert_plink12():
    recode = {'0': '1 1', '2': '2 2', '1': '1 2', '9': '0 0'}
    x=1
    # Read impute.snp and generate imputed.ped
    with open("impute.snp", 'r', encoding='utf-8') as f:
        with open("imputed.ped", 'w') as file_plink:
            for line in f:
                print(x)
                ids = line.split()[0]
                print(ids)
                geno=' '.join(line.split()[1:])
                genotype=[recode.get(geno[x],'') for x in range(0,len(geno))]
                zonta = altre_info(ids, genotype)
                file_plink.write(zonta)
                x += 1
    f.close()

def conv_map():
    t = pd.read_csv("impute.map", sep=" ", header=None)
    t.columns = ['Chr', 'SNPID', 'pos', 'BPPos']
    t["pos"] = "0"
    t.to_csv("imputed.map", sep=" ", index=False, header=False)
    

# Call the function to execute the conversion
convert_plink12()
conv_map()
os.system("plink2 --cow --ped imputed.ped --map imputed.map --make-bed --out imputed")
