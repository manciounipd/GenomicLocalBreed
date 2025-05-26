require("tidyverse")
require("data.table")
require("ggpubr")

  theme_set(theme_bw()+theme(#axis.line = element_line(),
          axis.text.x = element_text(size = 10, face = "italic",angle = 90, vjust = 0.5, hjust=1),
          strip.text.y = element_text(size = 12,face ="italic",angle=.9),
          strip.text.x = element_text(size = 12,face ="italic",angle=.9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          legend.position="none",
          axis.ticks.x =element_line(size = .5, colour = "black"),
          axis.line.x = element_line(size = .5, colour = "black"),
          axis.ticks.y =element_line(size = .5, colour = "black"),
          axis.line.y = element_line(size = .5, colour = "black"),
                axis.text.y = element_text(size = 10, face = "italic"),
          #strip.text.x = element_text(size=24, face="bold"),
          panel.background = element_blank()
  )
  )


cat("parte uno campioni\n")

fx=function(x) length(unique(x))
fxx <- function(x) {
  mean_x <- mean(x)
  ratio <- mean_x / sd(x)
  paste0(round(mean_x, 2), " (", round(ratio, 2), ")")
}


#dir="rend"

fx=function(dir) {
setwd(dir)

db=fread("data/db_with_test.txt")[,1] %>% distinct() 
ana=fread("data/ana.txt") [,c(1,4:5)]
train=fread("data/train.txt")
names(train) = "MATR"
test=fread("data/test.txt")
names(test) = "MATR"

names(ana) = c("MATR","YOB","SEX")
names(db) = c("MATR")

pheno=db[!db$MATR %in% test$MATR] # TOLGO I TEST PERCHE FITIZZIO

db=fread("model_ridotto_genomico/compa_with_info.txt")
names(db) = c("MATR","ebv_g","acc_g","ebv_p","acc_p","YOB")

db=merge(db,ana,by=c("MATR","YOB"))



tmp=db

cat("now add info")

geno = system("awk '{print $1}' data/impute.snp",intern=TRUE)

tmp$genotype="No"
tmp[tmp$MATR %in% geno,'genotype'] = "Yes"

tmp$type=as.character("train")
tmp[tmp$MATR %in% test$MATR,"type"] = "test" 
#tmp[tmp$MATR %in% train$MATR,"type"] = "train" 



tmp$phenotype="Ancestor of \nPhenotyped animal"
tmp[tmp$MATR%in% pheno$MATR,"phenotype"] = "Phenotyped Animal"
tmp[tmp$type=="test","phenotype"] = "Young Animal"

tmp %>% count(phenotype,genotype)

setwd("..")
return(tmp)

}


plotinfo=function(tmp,NOME) {
  
  
tmp1= tmp %>%
  group_by(phenotype,YOB, SEX,genotype) %>% 
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(YOB,SEX, phenotype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

head(tmp1)

tmp1 %>% count(phenotype)
tmp1[tmp1$SEX=="M","SEX"]="Male"
tmp1[tmp1$SEX=="F","SEX"]="Female"
# Create the plot
female=ggplot(tmp1%>% filter(SEX=="Female"),aes(x=YOB,y=count, fill=genotype)) +
  geom_col(alpha=.5) +
  facet_wrap(phenotype ~ ., scales="free") +
  #geom_text( data=tmp1 %>% filter(phenotype=="Phenotyped Animal" & genotype=="Yes"),
  #         aes(x=YOB,y=count, label=sprintf("%.0f%%", percentage)), color="black")  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + # fix
  xlab("")+ylab("")#+ggtitle("Female",element_text(hjust=0.5))


male=ggplot(tmp1%>% filter(SEX=="Male"),aes(x=YOB,y=count, fill=genotype)) +
  geom_col(alpha=.5) +
  facet_wrap(phenotype ~ ., scales="free") +
  #geom_text( data=tmp1 %>% filter(phenotype =="Ancestor of \nPhenotyped animal" & SEX=="M" & genotype=="Yes"),
   #        aes(x=YOB,y=count, label=sprintf("%.0f%%", percentage)), color="black")  +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
 #ggtitle("Male",element_text(hjust=0.5))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + # fix
  xlab("")+ylab("") 

if(NOME) {
p1=annotate_figure(ggarrange(female)+theme_bw() ,
                top = text_grob("Female",face = "italic", size = 18, color = "black"))  #, rot = 90) )
p2=annotate_figure(ggarrange(male)+theme_bw() , 
                top = text_grob("Male",face = "italic", size = 18, color = "black"))#, rot = 90) )

pp=ggarrange(p1,p2,ncol=2,common.legend = TRUE,legend="none")

} else {
   p1=annotate_figure(ggarrange(female)+theme_bw())
    p2=annotate_figure(ggarrange(male)+theme_bw() )
    pp=ggarrange(p1,p2,ncol=2,common.legend = TRUE,legend="none")
}

return(pp)
}




plot_acc=function(tmp,dir) {
setwd(dir)

tmp[tmp$SEX=="M","SEX"]="Male"
tmp[tmp$SEX=="F","SEX"]="Female"

db=fread("model_ridotto_genomico/compa_with_info.txt")
names(db) = c("MATR","ebvg","accg","ebvp","accp","YOB")

whole[whole$acc_p==0,"acc_p"]=0.00000000000000001

db=merge(db,tmp)

whole=db %>% mutate(gain=(accg-accp),gain2=(accg-accp)/accp)


gg=whole %>% group_by(YOB,phenotype,genotype,type,SEX)  %>%  filter(YOB > 2010 & n() > 5) %>%
 summarize(acc_median=median(gain2),accm=mean(gain2),qu=quantile(gain2, probs = 0.95,na.rm=TRUE),ql=quantile(gain2, probs = 0.05,na.rm=TRUE ))


gain1=whole %>% filter(YOB > 2010) %>% filter(genotype=="Yes")%>%
  ggplot(aes(x = as.factor(YOB), y = gain, fill = SEX)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  facet_wrap(phenotype ~ ., scales="free")+xlab("")+ylab("")+#+ylim(c(-0.1,0.5))
  geom_hline(yintercept=0,linetype = "dashed",alpha=.5)#+ ylab(expression("%")) 


gain2=gg %>% filter(genotype=="Yes") %>%
  ggplot(aes(x = as.factor(YOB), y = acc_median, fill = SEX)) +
  geom_col(position = position_dodge(width = 0.75), alpha = 0.5) +
  facet_wrap(phenotype ~ SEX, scales="free",nrow=1)+
  xlab("")+ylab("")+scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  geom_errorbar( aes(x=as.factor(YOB), ymin=ql, ymax=qu), width=0.1, alpha=0.8)
   



gain1a=whole %>% filter(YOB > 2010) %>% filter(genotype=="No")%>%
  ggplot(aes(x = as.factor(YOB), y = gain, fill = SEX)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  facet_wrap(phenotype ~ ., scales="free")+xlab("")+ylab("")+
  geom_hline(yintercept=0,linetype = "dashed",alpha=.5)+ylab(expression("%"))


gain2a=gg %>% filter(genotype=="Yes") %>%
  ggplot(aes(x = as.factor(YOB), y = acc_median, fill = SEX)) +
  geom_col(position = position_dodge(width = 0.75), alpha = 0.5) +
  facet_wrap(phenotype ~ SEX, scales="free",nrow=1)+
  xlab("")+ylab("")+scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  geom_errorbar( aes(x=as.factor(YOB), ymin=ql, ymax=qu), width=0.1, alpha=0.8)
   




setwd("..")
return(list(gain1,gain2,gain1a,gain2a))
}



grey=fx("grey")
vpr=fx("vpr")
vpn=fx("vpn")
rend=fx("rend")

cat("plot due change of accruacy ")


decrittivr=ggarrange(
    annotate_figure(plotinfo(grey,TRUE)+theme_bw(),left = text_grob("ALG",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(plotinfo(vpr,FALSE),left = text_grob("APR",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(plotinfo(vpn,FALSE),left = text_grob("ACN",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(plotinfo(rend,FALSE),left = text_grob("REN",face = "italic", size = 20, color = "black", rot = 90)),ncol=1
)

png("Figure1.png", units="in", width=15, height=10, res=300)
decrittivr
dev.off()

#top = text_grob("Visualizing Tooth Growth", color = "red", face = "bold", size = 14)
#final=annotate_figure(tot,top=text_grob("Female                                              Male", color = "black", face = "italic", size = 14))



gain_g=plot_acc(grey,"grey")#[[1]]
gain_vpr=plot_acc(vpr,"vpr")#[[1]]
gain_vpn=plot_acc(vpn,"vpn")#[[1]]
gain_rend=plot_acc(rend,"rend")#[[1]]

gain_1_geno=ggarrange(
    annotate_figure(gain_g[[1]]+theme(axis.text.x=element_blank()),left = text_grob("ALG",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_vpr[[1]]+theme(axis.text.x=element_blank()),left = text_grob("APR",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_vpn[[1]]+theme(axis.text.x=element_blank()),left = text_grob("ACN",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_rend[[1]],left = text_grob("REN",face = "italic", size = 20, color = "black", rot = 90)),ncol=1
)



png("Figure2.png", units="in", width=15, height=8, res=300)
gain_1_geno
dev.off()



gain_2_geno=ggarrange(
    annotate_figure(gain_g[[3]]+theme(axis.text.x=element_blank()),left = text_grob("ALG",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_vpr[[3]]+theme(axis.text.x=element_blank()),left = text_grob("APR",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_vpn[[3]]+theme(axis.text.x=element_blank()),left = text_grob("ACN",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_rend[[3]],left = text_grob("REN",face = "italic", size = 20, color = "black", rot = 90)),ncol=1
)



png("Figure2_Supllemenatary.png", units="in", width=1, height=10, res=300)
gain_2_geno
dev.off()


gain_vpr[[2]]+ theme( strip.text.x = element_blank() )
theme(strip.background = element_blank(),  strip.text.y = element_blank(),
strip.text = element_blank())


gain_2_geno=ggarrange(
    annotate_figure(gain_g[[2]],left = text_grob("ALG",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_vpr[[2]],left = text_grob("APR",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_vpn[[2]],left = text_grob("ACN",face = "italic", size = 20, color = "black", rot = 90)),
    annotate_figure(gain_rend[[2]],left = text_grob("REN",face = "italic", size = 20, color = "black", rot = 90)),ncol=1
)

png("Figure3.png", units="in", width=15, height=10, res=300)
gain_2_geno
dev.off()



####### now lr ########





get_lr=function(breed,effg) {

setwd(breed)

# metriche lr collexzion la tabelle
setwd("model_ridotto_genomico_lr")

LR=function(TOT) {
        bias=mean(as.double(TOT$solution-TOT$solution_w))
        disp=coef(lm(solution_w~ solution,data=TOT))[["solution"]]
        cor=cor(TOT$solution,TOT$solution_w)
        n=nrow(TOT)
        return(data.frame(y=dir,n=n,cor=cor,disp=disp,bias=bias))
}


wholeg=fread("solutions_w_geno") %>% filter(effect==effg)
wholeg=wholeg[,c("original_id","solution")]
names(wholeg)[2] = "solution_w"

wholep=fread("solutions_w_geno")%>% filter(effect==effg)
wholep=wholep[,c("original_id","solution")]
names(wholep)[2] = "solution_w"


dir=2013
x=1
collect=list()


geno = system("awk '{print $1}' ../data/impute.snp",intern=TRUE)

for(dir in 2013:2017) {

setwd(paste0("dir_",dir))


partial_geno=fread("solutions_p_geno") %>% filter(effect==effg)
partial_geno=partial_geno[,c("original_id","solution")]

partial_nogeno=fread("solutions_p_nogeno")%>% filter(effect==effg)
partial_nogeno=partial_nogeno[,c("original_id","solution")]

G=fread("test.txt",header=FALSE)$V1

F=merge(partial_geno,wholeg,by="original_id")

TOT=F[F$original_id %in% G &F$original_id %in% geno ,]
TOT2=F[F$original_id %in% G & ! F$original_id %in% geno ,]

F=merge(partial_nogeno,wholep,by="original_id")
TOT3=F[F$original_id %in% G &F$original_id %in% geno,]
TOT4=F[F$original_id %in% G & ! F$original_id %in%  geno,]

TMP=list(geno_SS=LR(TOT),peno_SS=LR(TOT2),geno_P=LR(TOT3),peno_P=LR(TOT4))
collect[[x]]=do.call("rbind",TMP)
x=x+1
setwd("..")

}

setwd("../..")
return(do.call("rbind",collect))

}




get_lr("rend",11)
get_lr("vpn",13)
get_lr("vpr",13)
get_lr("grey",15)

cat ("########## results lr  ###################")