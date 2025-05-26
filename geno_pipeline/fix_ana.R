#system("unzip anagrafepermungibilità2.zip ")
ana=data.table::fread("anagrafepermungibilitЕ2.csv",sep=";",quote="" )
ana=ana[,c(1,7,8,5)]
ana=as.data.frame(apply(ana,2,function(x) gsub(" ","_",x)))
names(ana) = c("ID","SIREID","DAMID","SEX")
ana=ana[ana$SEX %in% c("M","F"),]
ana[ana==""]="0"

mistake=which(ana$ID %in% ana$DAMID & ana$ID %in% ana$SIREID)


data.table::fwrite(file="ana.txt",ana,sep=" ",col.names=TRUE)

##############

ana=foreign::read.dbf("ANA0224.DBF")
ana=as.data.frame(apply(ana,2,function(x) gsub(" ","_",x)))
ana[ana=="000000000000000"]="0"
ana=ana[,c(1:3,6)]
names(ana) = c("MATR","PADRE","MADRE","SEX")
data.table::fwrite(file="ana.txt",ana,sep=" ",col.names=TRUE)

# fino a qua 


