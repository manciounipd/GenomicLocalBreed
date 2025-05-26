
#Rscript fai il plot cosi vedo come è struttrata la popolazion

snpfile="data/impute.snp"
geno = system(paste0("awk '{print $1}' ",snpfile), intern = TRUE)

db=data.table::fread("data/db.txt")
ana=data.table::fread("data/ana.txt")
names(ana) = c("anim","fid","mid","yob","sex")

ana$pheno="no"
ana[ana$anim %in% unique(db$MATR),"pheno"] = "yes"

geno=as.data.frame(geno)
names(geno) = "anim"
geno$geno = "yes"
ana=merge(geno,ana,all=TRUE)
ana[is.na(ana$geno),"geno"] ="no"

table(ana$geno)

###### fai i plot ##########

require("tidyverse")
ana %>% filter(geno=="yes") %>% count(sex,yob)
table(ana$pheno)

# torva genitori con il fenotipo cioè quelli utli
ste = ana[ana$pheno =="yes",]

found_g = function(i) {
            
            return(c(ana[ana$anim==i,]$fid,ana[ana$anim==i,]$mid))
}

traceback = function(i) {
    tot=list()
    x=1
    start=i
    repeat  {
        i=found_g(i)
        i=i[!i%in%"0"]
        if(length(i)==0) {
            break
        }
        tot[[x]]=data.frame(anim=i,gen=x,orig=start)
        x=1+x
        
    }
    return(do.call("rbind",tot))
}
# questa utile per singolo animale che mi becca ttt
tracebackf = function(i) {
    tot=list()
    x=1
    start=i
    repeat  {
        i=found_g(i)
        i=i[!i%in%"0"]
        if(length(i)==0) {
            break
        }
        tot[[x]]=i#data.frame(anim=i,gen=x,orig=start)
        x=1+x
        
    }
    return(do.call("rbind",tot))
}

############################################
#    traceback iformation by generation
#############################################


my_renumf90 = function(s) {

found_g1 = function(i) {
            tmp=ana[ana$anim%in%i,]
            return(unique(c(tmp$fid,tmp$mid)))
}

dall = list()
maxgen=100
x=0
ntot=s
repeat {
    cat(length(s),": ")
    dall[[x+1]] = data.frame(s,x)
    ns=found_g1(s)
    ns=ns[!ns%in% ntot]
    ns=ns[!ns=="0"]
    cat(x,"\n")
    ntot = c(ns,ntot)
    s=ns
    x=x+1
    if(x == maxgen | length(s)==0) {break}
     
}

return(list(ntot,do.call("rbind",dall)))

}



anim=ste$anim
f=my_renumf90(anim)
ped=f[[2]]
names(ped)  = c("anim","gen")
ped=merge(ped,ana,all.y=TRUE)
head(ped)
ped[is.na(ped$gen),"gen"] = "not-informative"

table(ped[ped$geno=="yes",]$yob)



png("all_genotype.png")
ped %>% filter( geno=="yes") %>% ggplot(aes(x=yob,fill=sex))+
geom_bar()+theme_light()
dev.off()


ped[(ped$gen) %in% as.character(3:9), "gen"] = "3-9"
# Basic chart

fr=ped  %>% filter(geno=="yes") %>% count(sex,gen,yob)


A=fr %>% ggplot(aes(x=yob,y=n,fill=gen))+
geom_col()+
facet_wrap(~sex, ncol = 1)+
geom_vline(xintercept=2016)+
ylab("Number of animals")+
theme_light()+
theme(legend.position = "bottom")+
scale_fill_discrete(name = "TYPE OF \nANIMAL:")+
xlab("Y.O.B.")

B=ped %>% filter(gen != "not-informative" & yob > 2000) %>% ggplot(aes(x=as.integer(yob),fill=geno))+
geom_bar()+ facet_grid(gen~sex,scale="free_y")+
theme_light()+
theme(legend.position = "bottom")+
scale_fill_discrete(name = "genotyped:")+
xlab("Y.O.B.")



png("informative_complete.png")
ggpubr::ggarrange(A,B)
dev.off()

####################################################

x=1
require(parallel)
r=mclapply2(ste$anim,found_g , mc.cores = 24)
length(unique(do.call("c",r)))


require(parallel)

r=mclapply2(ste$anim,tracebackf , mc.cores = 26)
tot=unique(do.call("c",r))
qanti=tot %>% count(anim)

mindist = tot %>% group_by(anim) %>% filter(gen ==min(gen)) 
nrow(mindist)
length(ste$anim) + nrow(mindist)
length(ste$anim) + length(unique(tot$anim))
head(mindist)
save.image("t.RData")
