
require("tidyverse");require("data.table")

df_r = foreign::read.dbf("data/controlliscs.dbf")
df_r$MATR = as.character(df_r$MATR)
df_r$MATR = gsub(" ","_",df_r$MATR)
df_r = df_r %>% select(MATR,AZDCNL,CLGRAV,W0,L1,L2,L3,CLNLETAP,CLNLMP,LATTE,GRPER,PRPER)
head(df_r)
nrow(df_r);length(unique(paste0(df_r$MATR,df_r$NL)));length(unique(df_r$MATR))
        

ana=foreign::read.dbf("data/ANA0224.DBF")[,c(1,2,3,5,6)]
ana$DN = substr(ana$DN, 1, 4)
ana=as.data.frame(apply(ana,2,function(x) gsub(" ","_",x)))

df_r$MATR=gsub(" "," ",df_r$MATR) 
fwrite(df_r,"data/db.txt",sep=" ",quote=FALSE)


ana[ana=="000000000000000"]="0"
fwrite(ana,"data/ana.txt",sep=" ",quote=FALSE,col.names = FALSE)


