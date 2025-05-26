#!/usr/bin/env Rscript

fc=function(tmp,colu) {
    ix=nrow(tmp)
    tm=names(table(tmp[,colu])[table(tmp[,colu]) == 1])
    tmp=tmp[!tmp[,colu] %in% tm,]
    #cat(ix-nrow(tmp),"\n")
    return(tmp)
}

my_renumf90 = function(s,ana,maxgen=1000) {

found_g1 = function(i) {
            tmp=ana[ana$MATR%in%i,]
            return(unique(c(tmp$FID,tmp$FID)))
}

dall = list()
x=0
ntot=s
repeat {
  #  cat(length(s),": ")
    dall[[x+1]] = data.frame(s,x)
    ns=found_g1(s)
    ns=ns[!ns%in% ntot]
    ns=ns[!ns=="0"]
 #   cat(x,"\n")
    ntot = c(ns,ntot)
    s=ns
    x=x+1
    if(x == maxgen | length(s)==0) {break}
     
}

return(list(ntot,do.call("rbind",dall)))

}


fc_loop = function(tmp,eff) {
repeat {
    n=nrow(tmp)
    for (eff_i in eff) {
            tmp = fc(tmp,eff_i)
    }
    if(nrow(tmp)==n) {
        break
    }
}
print(nrow(tmp))
return(tmp)
}


year= 2007
snpfile="data/impute.snp"

cat("anno: ",year,"\nsnpfile: ",snpfile,"\n")

db=data.table::fread("data/db.txt")
ana=data.table::fread("data/ana.txt")

names(ana) = c("MATR","FID","MID","YOB","SEX")

db=merge(ana[,c("MATR","YOB")],db,by.x="MATR",by.y="MATR")

#print(table(db$YOB))

tmp=db[db$YOB > year,]
tmp=as.data.frame(tmp)

cat("rimuovi campi unici \n")

tmp=fc_loop(tmp,eff=c("AZDCNL" , "CLGRAV"))
db_exp=tmp[,-2]

cat("n animals",length(unique(tmp$MATR)),"\n")

geno = system(paste0("awk '{print $1}' ",snpfile), intern = TRUE)

cat(length(geno),"\n")


cat("keep geneotype parents \n")
F=my_renumf90(unique(db$MATR),ana,5)
ge=F[[1]][F[[1]] %in%  geno ]

genopuned_correctid= ge
data.table::fwrite(as.data.frame(ge),"data/train.txt")
cat("test animal ",length(ge),"\n")
cat("now add young annimals > 2016\n")

YU=ana[ana$MATR%in%geno,]
print(summary(YU$YOB))
YU=YU[YU$YOB  > 2016,]
YU=YU[!YU$MATR %in% genopuned_correctid,]
cat("young animal",nrow(YU),"\n")

genopuned_correctid=as.data.frame(unique(c(genopuned_correctid,YU$MATR)))
data.table::fwrite(file=paste0("data/p_",year,".txt"),genopuned_correctid,sep=" ")

data.table::fwrite(YU[,1],"data/test.txt")


cat("now subset file    \n")
system(paste0("sort -k 1  ",snpfile," > tmp1"))
system(paste0("sort -k 1 data/p_",year ,".txt > tmp2"))
system(paste0("join -1 1 -2 1 tmp1 tmp2 > data/pruned.snp"))
system("rm tmp1 tmp2")
system("wc -l data/pruned.snp")


tmp=data.frame(MATR=c(YU[,1]))
names(tmp) = "MATR"
complte=merge(tmp,db_exp,all=TRUE)
complte[is.na(complte)]  = "0"

data.table::fwrite(complte,paste0("data/db_with_test.txt"),sep=" ")


# fare una tabella con due collonne test and triang geneotuye #





if(FALSE) {

############# prova LR #########

db=alldb[alldb$latte!=0, ]
db=merge(db,ana[,c(1,4)])
db[db$YOB > 2017,"latte"] =0

data.table::fwrite(db[,-ncol(db)],paste0("data/LR_db_with_",x,"_",year,".txt"),sep=" ")
#########

dba=db[db$YOB > 2017,"MATR"]
test=unique(dba)

partial=data.table::fread("LRmodel1_red_2010_NOgenomic/solutions.orig")
whole=data.table::fread("model1_red_2010NO_genomic/solutions.orig")

whole=whole[whole$effect %in% 13,c(4,5)]
partial=partial[partial$effect %in% 13,c(4,5)]

tot=merge(partial,whole,by="original_id")

tetsint=tot[tot$original_id %in%  test,]

cor(tetsint[,2:3])
mean(tetsint$solution.x)-mean(tetsint$solution.y)
lm(tetsint$solution.y~tetsint$solution.x)



partial=data.table::fread("LRmodel1_red_2010_genomic/solutions.orig")
whole=data.table::fread("model1_red_2010_genomic/solutions.orig")

whole=whole[whole$effect %in% 13,c(4,5)]
partial=partial[partial$effect %in% 13,c(4,5)]

tot=merge(partial,whole,by="original_id")

tetsint=tot[tot$original_id %in%  test,]

cor(tetsint[,2:3])
mean(tetsint$solution.x)-mean(tetsint$solution.y)
lm(tetsint$solution.y~tetsint$solution.x)

}