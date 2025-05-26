 
getwd()

year=2010

setwd("PROGETTI/acc_rend")

require("tidyverse");require("data.table")

        df_r = foreign::read.dbf("DATA/controlliscs.dbf")
        df_r$MATR = as.character(df_r$MATR)
        df_r$MATR = gsub(" ","_",df_r$MATR)
        head(df_r)
        df_r = df_r %>% select(MATR,AZDCNL,CLGRAV,W0,L1,L2,L3,CLNLETAP,CLNLMP,PRPER )
        head(df_r)
        #df_r$LATTE = scale(df_r$LATTE)
        nrow(df_r);length(unique(paste0(df_r$MATR,df_r$NL)));length(unique(df_r$MATR))
        
        ana=data.table::fread("DATA/ana.txt")[,c(1,4)]
        names(ana) = c("MATR","YOB")
        df_r=merge(ana,df_r) %>% filter(YOB > year )
        
        data.table::fwrite(file="fdb.txt",df_r[,-2],sep= " ",col.names = TRUE,quote = FALSE)

        
            fc=function(tmp,colu) {
                ix=nrow(tmp)
                tm=names(table(tmp[,colu])[table(tmp[,colu]) == 1])
                tmp=tmp[!tmp[,colu] %in% tm,]
                cat(ix-nrow(tmp),"\n")
                return(tmp)
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

            tmp=fc_loop(as.data.frame(df_r),eff=c( "AZDCNL" ))

system("cp /home/enrico/GENOTIPI/GENOTIPI_RENDENA/DATA/data_ready_after_impute/imputed_f901.txt ../snp.txt")
/home/enrico/GENOTIPI/GENOTIPI_RENDENA/DATA/imputation_steps/imp-txt
snpfile="snp.txt"
            #data.table::fwrite(file=paste0("data/db_vpr",year,".txt"),tmp[,-2],sep=" ")

            db_exp=tmp[,-2]
ana=data.table::fread("DATA/ana.txt")
        names(ana) = c("MATR","PS","SD","YOB")


            geno = system(paste0("awk '{print $1}' ",snpfile), intern = TRUE)
            recode=ggroups::renum(ana[,1:3])

            IX=recode$xrf[recode$xrf$ID %in% unique(tmp$MATR),]$newID
            pruned=ggroups::pruneped(recode$newped, IX, mode="strict")

            genoID=recode$xrf[recode$xrf$ID %in% geno,]$newID

            genopuned = pruned[pruned$ID %in% genoID,]$ID

            genopuned_correctid= as.data.frame(recode$xrf[recode$xrf$newID %in% genopuned,]$ID)

        print(getwd())
        system("head fdb.txt")
        system("tail fdb.txt")

head(df_r)


print("add young annimals")
YU=ana[ana$MATR%in%geno,]
print("summary YOB gen ..")
print(summary(YU$YOB))
YU=YU[YU$MATR  %in% db_exp$MATR  & YU$YOB  > 2016,]
genopuned_correctid=as.data.frame(unique(c(genopuned_correctid[,1],YU$MATR)))


data.table::fwrite(file=paste0("pruned_snpanimal_",year,".txt"),genopuned_correctid,sep=" ")


pruned_snpanimal_file <- paste0("pruned_snpanimal_",year,".txt")
output_file <- paste0(gsub(".txt","",snpfile),"_",year,".txt")



# Construct the shell command

system(paste0("sort -k 1  ",snpfile," > tmp1"))
system(paste0("sort -k 1 ",pruned_snpanimal_file ," > tmp2"))
system(paste0("join -1 1 -2 1 tmp1 tmp2 > ",output_file))


# Execute the command

system("rm tmp1 tmp2")


tmp=data.frame(MATR=c(geno))
names(tmp) = "MATR"
alldb=merge(tmp,db_exp,by="MATR",all=TRUE)
cat(nrow(alldb) - nrow(db_exp),"\n")
alldb=as.data.frame(apply(alldb,2,as.character))
#alldb=a
alldb[is.na(alldb)] = "0"
head(alldb)

x=length(strsplit(snpfile,"/",fixed=T)[[1]])
x=strsplit(snpfile,"/",fixed=T)[[1]][x]
x=(gsub(".txt","",x))
print(nrow(alldb)- nrow(db_exp))
data.table::fwrite(alldb,paste0("db_with_",x,"_",year,".txt"),sep=" ",quote =FALSE)


head(alldb)


dir.create("vce")
setwd("vce")
getwd()
system("head ../db_with_snp_2010.txt")


system(paste0("echo 'DATAFILE
      ../db_with_snp_2010.txt
       SKIP_HEADER
        1
       TRAITS
       10
       FIELDS_PASSED TO OUTPUT

       WEIGHT(S)
       0
       RESIDUAL_VARIANCE # chidere A CRI
      10.000   
       EFFECT
       2  cross alpha
       EFFECT
       3   cross alpha     
       EFFECT
       4 cov       
       EFFECT
       5 cov       
       EFFECT
       6 cov       
       EFFECT
       7 cov       
       EFFECT
       8 cross alpha    
       EFFECT
       9 cross alpha        
       EFFECT
       1  cross alpha # animale
       RANDOM
       animal   # pruned ped!
       OPTIONAL
       pe
       FILE
      ../DATA/ana.txt
       FILE_POS
       1 2 3
       PED_DEPTH
       15
       (CO)VARIANCES
  0.99051      
       (CO)VARIANCES_PE
   0.4000       
       OPTION EM-REML 100
       OPTION use_yams
    OPTION conv_crit 1d-10
       OPTION method VCE 
       OPTION remove_all_missing' > par.txt "))
getwd()
       
       system("echo par.txt | renumf90 ")
        

        system("cat renf90.par")



system(paste0("echo 'DATAFILE
 renf90.dat
NUMBER_OF_TRAITS
           1
NUMBER_OF_EFFECTS
          12
OBSERVATION(S)
    1    
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
 2        57573 cross 
3         18 cross
 4         36 cross 9
 5         36 cross 9
 6         36 cross 9
 7         36 cross 9
 4         42 cross 8
 5         42 cross 8
 6         42 cross 8
 7         42 cross 8
10     36983 cross 
10      36983 cross
RANDOM_RESIDUAL VALUES
0.51743E-01
 RANDOM_GROUP
    11
 RANDOM_TYPE
 add_animal   
 FILE
renadd09.ped                                                                    
(CO)VARIANCES
0.25862E-01
 RANDOM_GROUP
    12
 RANDOM_TYPE
 diagonal     
 FILE
                                                                                
(CO)VARIANCES
 0.96405E-02
#OPTION EM-REML 100
OPTION original id
OPTION use_yams
OPTION method VCE 
OPTION sol se
OPTION conv_crit 1d-10' > par_fix.txt "))


system("ulimit -s unlimited; echo par_fix.txt | blupf90+")

#########
setwd("..")


dir.create("blup")
setwd("blup")
getwd()



system(paste0("echo 'DATAFILE
      ../db_with_snp_2010.txt
       SKIP_HEADER
        1
       TRAITS
       10
       FIELDS_PASSED TO OUTPUT

       WEIGHT(S)
       0
       RESIDUAL_VARIANCE # chidere A CRI
      10.000   
       EFFECT
       2  cross alpha
       EFFECT
       3   cross alpha     
       EFFECT
       4 cov       
       EFFECT
       5 cov       
       EFFECT
       6 cov       
       EFFECT
       7 cov       
       EFFECT
       8 cross alpha    
       EFFECT
       9 cross alpha        
       EFFECT
       1  cross alpha # animale
       RANDOM
       animal   # pruned ped!
       OPTIONAL
       pe
       FILE
      ../DATA/ana.txt
       FILE_POS
       1 2 3
       PED_DEPTH
       5
       (CO)VARIANCES
  0.99051      
       (CO)VARIANCES_PE
   0.4000       
       OPTION use_yams' > par.txt "))
getwd()
       
       system("echo par.txt | renumf90 ")
        

        system("cat renf90.par")



system(paste0("echo 'DATAFILE
 renf90.dat
NUMBER_OF_TRAITS
           1
NUMBER_OF_EFFECTS
          12
OBSERVATION(S)
    1    
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
 2        57573 cross 
3         18 cross
 4         36 cross 9
 5         36 cross 9
 6         36 cross 9
 7         36 cross 9
 4         42 cross 8
 5         42 cross 8
 6         42 cross 8
 7         42 cross 8
 10     14334 cross 
 10      14334 cross
RANDOM_RESIDUAL VALUES
0.52081E-01
 RANDOM_GROUP
    11
 RANDOM_TYPE
 add_animal   
 FILE
renadd09.ped                                                                    
(CO)VARIANCES
0.24427E-01
 RANDOM_GROUP
    12
 RANDOM_TYPE
 diagonal     
 FILE
                                                                                
(CO)VARIANCES
 0.13029E-01
#OPTION EM-REML 100
OPTION original id
OPTION use_yams
OPTION sol se
OPTION conv_crit 1d-10' > par_fix.txt "))


system("ulimit -s unlimited; echo par_fix.txt | blupf90+")

###############################
setwd("..")


dir.create("blup_g")
setwd("blup_g")
system("awk '{ if (NR>0) { gsub(/9/, \"5\", $2) } print }' ../snp.txt > ok.txt")



system(paste0("echo 'DATAFILE
      ../db_with_snp_2010.txt
       SKIP_HEADER
        1
       TRAITS
       10
       FIELDS_PASSED TO OUTPUT

       WEIGHT(S)
       0
       RESIDUAL_VARIANCE # chidere A CRI
      10.000   
       EFFECT
       2  cross alpha
       EFFECT
       3   cross alpha     
       EFFECT
       4 cov       
       EFFECT
       5 cov       
       EFFECT
       6 cov       
       EFFECT
       7 cov       
       EFFECT
       8 cross alpha    
       EFFECT
       9 cross alpha        
       EFFECT
       1  cross alpha # animale
       RANDOM
       animal   # pruned ped!
       OPTIONAL
       pe
       FILE
      ../DATA/ana.txt
       FILE_POS
       1 2 3
       SNP_FILE
       ok.txt
       PED_DEPTH
       5
       (CO)VARIANCES
  0.99051      
       (CO)VARIANCES_PE
   0.4000       
       OPTION use_yams' > par.txt "))

getwd()
       
       system("echo par.txt | renumf90 ")
        

        system("cat renf90.par")



system(paste0("echo 'DATAFILE
 renf90.dat
NUMBER_OF_TRAITS
           1
NUMBER_OF_EFFECTS
          12
OBSERVATION(S)
    1    
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
 2        57573 cross 
3         18 cross
 4         36 cross 9
 5         36 cross 9
 6         36 cross 9
 7         36 cross 9
 4         42 cross 8
 5         42 cross 8
 6         42 cross 8
 7         42 cross 8
10     36983 cross 
10      36983 cross
RANDOM_RESIDUAL VALUES
0.52081E-01
 RANDOM_GROUP
    11
 RANDOM_TYPE
 add_animal   
 FILE
renadd09.ped                                                                    
(CO)VARIANCES
0.24427E-01
 RANDOM_GROUP
    12
 RANDOM_TYPE
 diagonal     
 FILE
                                                                                
(CO)VARIANCES
 0.13029E-01
OPTION use_yams
OPTION sol se
OPTION SNP_file ok.txt
OPTION conv_crit 1d-10' > par_fix.txt "))


system("ulimit -s unlimited; echo par_fix.txt | blupf90")


setwd("..")

require("data.table")
require("tidyverse")

getwd()

acc = function(va,se) {return(sqrt(1 - ((se)**2 / va)))}

gef=11
va=0.24427E-01

dat = fread("db_with_snp_2010.txt") %>% filter(PRPER!=0)

db1=fread("blup_g/solutions",skip=1)
db1=db1[db1$V2==gef,c(3,4,5)]
head(db1)
names(db1) = c("id","ebv_g","se_g")
ped=fread("blup_g/renadd09.ped")[,c(1,10)]
names(ped) = c("id","matr")
db1=merge(db1,ped)[,-1]
head(db1)




db2=fread("blup/solutions",skip=1)
db2=db2[db2$V2==gef,c(3,4,5)]
names(db2) = c("id","ebv_p","se_p")
ped=fread("blup/renadd09.ped")[,c(1,10)]
names(ped) = c("id","matr")
db2=merge(db2,ped)[,-1]

head(db1)
head(db2)
summary(db1)
summary(db2)

acc = function(va,se) {return(sqrt(1 - ((se)**2 / va)))}


dbm = merge(db1,db2,by="matr",all=TRUE)


dbm$acc_g = acc(va,dbm$se_g)
dbm$acc_p = acc(va,dbm$se_p)
dbm[is.na(dbm)]=0

nrow(dbm);nrow(db1);nrow(db2)


ana = foreign::read.dbf("DATA/ANA0224.DBF")
mm=ana[,]
head(mm)
names(mm) = c("matr","fid","mid","dn","dna","sex")
mm$matr = gsub(" ","_",mm$matr)
f = merge(mm,dbm)
nrow(f) # perdo 34 giovani animali

geno=system("awk '{print $1 }' blup_g/ok.txt",intern=TRUE)

f$tag="traced-back"
f[f$matr %in% geno,"tag"] = "geno"
f[f$matr %in% dat$MATR,"tag"] = "pheno"
f[f$matr %in% geno  & f$matr %in% dat$MATR,"tag"] = "pheno-geno"

tail(f)
nrow(f)
f %>% filter(is.na(se_g)) %>% filter(yob > 2020)
f$yob=year(as.Date(f$dna,"%Y-%m-%d"))

i=f %>%  filter(sex=="M") %>% pluck("matr") %>% as.data.frame()
names(i) = "matr"
sid=f %>% filter(tag %in% c("pheno","pheno-geno") & fid %in% matr) %>%  count(fid)
head(sid)
sid=merge(i,sid,by.x="matr",by.y="fid",all=TRUE)
head(sid)
sid[is.na(sid)] =0
info=merge(sid,f,by="matr")
head(info)


info %>% filter(tag=="geno" & yob > 2017) %>% #reshape2::melt(id=c("matr","yob","fid","mid","n","dna","sex","dn","ebv_g","ebv_p","se_p","se_g")) %>% 
            group_by(n,yob) %>% summarize(n_anim=n(),mg=mean(acc_g),sd_g = sd(acc_g),mp=mean(acc_p),sd_p = sd(acc_p)) #%>% reshape2::melt(id="n")
              geom_bar( aes(x=n, ,fill=variable,y=value), stat="identity", fill="skyblue", alpha=0.7) 
    geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="orange",   alpha=0.9, size=1.3)

gganimate::animate(zx, height = 500, width = 800, fps = 30, duration = 10,
        end_pause = 60, res = 100)
gganimate::anim_save("ps3 game sales.gif")

jpeg("u.jpg",width = 12, height =7,unit="in",res=300)

#ggpubr::ggarrange(

f %>% filter(yob > 2010 ) %>% filter(sex=="F") %>% 
mutate(dam=ifelse(matr %in% mid,"dam","notdam")) %>% 
ggplot(aes(x=as.factor(yob),fill=tag,group=as.factor(yob),y=(acc_g-acc_p)/acc_p))+
geom_boxplot()+
facet_grid( dam ~ tag,scale="free_y")+
geom_hline(yintercept=0,color="red")+
theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
 xlab("YOB")
dev.off()

jpeg("u22.jpg",width = 12, height =7,unit="in",res=300)

f %>% filter(yob > 2010 & yob < 2022) %>% filter(sex=="F") %>% 
mutate(dam=ifelse(matr %in% mid,"dam","notdam")) %>% 
ggplot(aes(x=as.factor(yob),fill=tag,group=as.factor(yob),y=(acc_g-acc_p)/acc_p))+
geom_boxplot()+
facet_grid( dam ~ tag,scale="free_y")+
geom_hline(yintercept=0,color="red")+
theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
 xlab("YOB")

dev.off()


jpeg("sire.jpg",width = 12, height =7,unit="in",res=300)

#ggpubr::ggarrange(

f %>% filter(yob > 2010 ) %>% filter(sex=="M") %>% 
mutate(dam=ifelse(matr %in% fid,"sire","notsire")) %>% 
ggplot(aes(x=as.factor(yob),fill=tag,group=as.factor(yob),y=(acc_g-acc_p)/acc_p))+
geom_boxplot()+
facet_grid( dam ~ tag,scale="free_y")+
geom_hline(yintercept=0,color="red")+
theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
 xlab("YOB")

dev.off()

jpeg("sire22.jpg",width = 12, height =7,unit="in",res=300)

f %>% filter(yob > 2010 & yob < 2022)   %>% filter(sex=="M") %>% 
mutate(dam=ifelse(matr %in% fid,"sire","notsire")) %>% 
ggplot(aes(x=as.factor(yob),fill=tag,group=as.factor(yob),y=(acc_g-acc_p)/acc_p))+
geom_boxplot()+
facet_grid( dam ~ tag,scale="free_y")+
geom_hline(yintercept=0,color="red")+
theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
 xlab("YOB")

dev.off()

library(ggplot2)
library(ggpmisc)

jpeg("corr.jpg",width = 12, height =7,unit="in",res=300)
ggplot(data = f %>% filter(yob > 2010 ), aes(x = ebv_g, y = ebv_p)) +
geom_point(aes(x = ebv_g, y = ebv_p,fill=yob,color=yob))+
  stat_poly_line(color="red") +
  stat_poly_eq() +
  facet_grid(  ~ sex)+
  theme_bw()
dev.off()

jpeg("corr2.jpg",width = 12, height =7,unit="in",res=300)
ggplot(data = f %>% filter(yob > 2010 & yob < 2019 ), aes(x = ebv_g, y = ebv_p)) +
geom_point(aes(x = ebv_g, y = ebv_p,fill=yob,color=yob))+
  stat_poly_line(color="red") +
  stat_poly_eq() +
  facet_grid(  ~ sex)+
  theme_bw()
dev.off()







xc=mm
xc$type = "not_genotyped"
xc[xc$matr %in% geno,"type"] = "genotyped"
xc$yob=year(as.Date(xc$dna,"%Y-%m-%d"))

jpeg("stren.jpg",width = 13, height =7,unit="in",res=300)
my_colors <- c("red","grey")

ggpubr::ggarrange(

xc %>% filter(matr %in% fid)  %>% group_by(yob,type,sex) %>% 
filter(yob > 1980) %>% count() %>%
        ggplot(aes(yob,n,fill=type))+
        geom_col()+
        theme_bw()+
          scale_fill_manual(values = my_colors),

xc %>% filter(sex=="F")  %>% group_by(yob,type,sex) %>% 
filter(yob > 1980) %>% count() %>%
        ggplot(aes(yob,n,fill=type))+
        geom_col()+
        theme_bw()+
          scale_fill_manual(values = my_colors),
          labels=c("A","B"),
          common.legend = TRUE, legend="bottom"

)
dev.off()

xc %>% group_by(Chip,SEX) %>% filter(Chip!="no_geno") %>%
filter(YOB > 1980) %>% count() %>%  
ggplot(aes(SEX,n,fill=Chip))+
        geom_col()+
        theme_bw()
