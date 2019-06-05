#Load Data
SD_algae.df<- read.delim("nbc23s.final.subsample.tax.summary")
SD_algae.df<-SD_algae.df[,1:71]

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.us.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
# usage
packages <- c("ggplot2", 
              "plyr", 
              "reshape2", 
              "RColorBrewer", 
              "scales", 
              "grid",
              "VennDiagram",
              "vegan",
              "RAM",
              "tidyr",
              "Rmisc",
              "reshape","reshape2",
              "gapminder","magrittr","dplyr","ggpubr","gridExtra",
              "patternplot","tibble","gplots","broom","data.table","devtools")
ipak(packages)

pd=position_dodge(width=0.4)

#Remove Chloroplast
SD_algae.df=SD_algae.df[SD_algae.df$taxon !="Chloroplast_ge",]

#SHANNON DIVERSITY --> Subset to taxlevel=6
SD_algae.df$taxlevel=as.factor(SD_algae.df$taxlevel)
SD_algae.df.6=subset.data.frame(SD_algae.df, taxlevel==6,drop=TRUE)

SD_alg.shannon.df<-SD_algae.df.6 %>% 
  subset(select=c(-taxlevel,-rankID,-taxon,-daughterlevels,-total)) %>%
  mutate_at(vars(X170818FA:X170915TC),funs(diversity(.,index="shannon"))) %>%
  gather(key="Sample",value="Diversity") %>%
  distinct(Sample,Diversity) %>% 
  separate(Sample, into=c("junk","Date","Substrate","Rep"), sep=c(1,7,-1))

SD_alg.shannon.df$junk=NULL
SD_alg.shannon.df=droplevels(SD_alg.shannon.df[-which(SD_alg.shannon.df$Substrate=='CR'),]) #drop controls
SD_alg.shannon.df$Date=as.factor(SD_alg.shannon.df$Date)
SD_alg.shannon.df$Date=revalue(SD_alg.shannon.df$Date,c("170818"="3","170828"="10","170905"="21","170915"="31"))
SD_alg.shannon.df$Substrate=factor(SD_alg.shannon.df$Substrate,levels=c("H","SO","SH","F","T"),labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))

#Summarize Shannon
SD_shannon.div.summary<-SD_alg.shannon.df %>% 
  group_by(Date, Substrate) %>% 
  summarise(mean=mean(Diversity),se=(sd(Diversity)/sqrt(n())))

#Shannon Figure
SD_shannon.div.alg=ggplot(SD_shannon.div.summary, aes(x=Date,y=mean,group=Substrate,shape=Substrate,fill=Substrate))+
  geom_point(position=pd,size=4)+
  geom_line(aes(linetype=Substrate),position=pd)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.3,position=pd)+
  ylab("Shannon Diversity (H')")+
  theme_classic()+
  theme(legend.position=c(0.25,0.95),
        legend.title=element_blank(),
        legend.direction="horizontal",
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14))+
  scale_shape_manual(name="Substrate",values=c(21,22,23,24,25,26))+
  ggtitle("Algal Communities")+
  scale_y_continuous(breaks=seq(0,5,0.5))
SD_shannon.div.alg

#OTU Richness
SD_algae.df$taxlevel=as.factor(SD_algae.df$taxlevel)
SD_algae.df.6=subset.data.frame(SD_algae.df, taxlevel==6,drop=TRUE)
observationThreshold = 1 #set min number to count as OTU being present
SD_otu=apply(SD_algae.df.6[,6:67]>=observationThreshold, 2, sum,na.rm=T) #check by row to see if each taxa present or not
SD_otu.0=apply(SD_algae.df.6[,6:67],2,sum,na.rm=T)
list_to_df <- function(list_for_df) {
  list_for_df <- as.list(list_for_df)
  nm <- names(list_for_df)
  if (is.null(nm))
    nm <- seq_along(list_for_df)
  df <- data.frame(name = nm, stringsAsFactors = FALSE)
  df$value <- unname(list_for_df)
  df
}
SD_otu.df=list_to_df(SD_otu) #conevrt list to data frame wiht custom function above
SD_otu.df$value=as.numeric(SD_otu.df$value)
SD_otu.df=separate(SD_otu.df,name,into=c("junk","Date","Substrate"),sep=c(1,7))
SD_otu.df=separate(SD_otu.df,Substrate,into=c("Substrate","Rep"),sep=c(-1))
SD_otu.df$junk=NULL
SD_otu.df=droplevels(SD_otu.df[-which(SD_otu.df$Substrate=='CR'),])
SD_otu.df$Date=as.factor(SD_otu.df$Date)
SD_otu.df$Date=revalue(SD_otu.df$Date,c("170818"="3","170828"="10","170905"="21","170915"="31"))
SD_otu.aov=SD_otu.df%>% #tidy pipes!  Run anovs by Date on otu richenss
  group_by(Date)%>%
  do(tidy(aov(value~Substrate,data=.)))
SD_otu.aov
SD_tukey.out=SD_otu.df%>%
  group_by(Date)%>%
  do(multitst=TukeyHSD(aov(value~Substrate,data=.)))
SD_tukey.out %>%tidy(multitst)
SD_otu.df.summary=summarySE(SD_otu.df,measurevar="value",groupvars=c("Date","Substrate"),na.rm=TRUE)
pd=position_dodge(width=0.4)
SD_otu.plot=ggplot(SD_otu.df.summary, aes(x=Date,y=value,group=Substrate,shape=Substrate,fill=Substrate))+
  geom_point(position=pd,size=4)+geom_line(aes(linetype=Substrate),position=pd)+
  geom_errorbar(aes(ymin=value-se,ymax=value+se),width=0.3,position=pd)+
  ylab("Observed OTUs")+
  theme_classic()+
  theme(legend.position=c(0.25,0.95),
        legend.title=element_blank(),
        legend.direction="horizontal",
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14))+
  scale_shape_manual(name="Substrate",values=c(21,22,23,24,25,26))
SD_otu.plot

#Relative Abundance
theme_set(theme_classic(base_size = 18))
SD_algae.df.3=subset.data.frame(SD_algae.df, taxlevel==3,drop=TRUE) #taxlevel 3 is class
sites=list(colnames(SD_algae.df.3[,6:67])) #get list of samples
total=list(colnames(SD_algae.df.3[,1]))
#custome relative abudnace function/  it takes a list of sites and the total column for each taxa and finds the relaive abundance of each taxa in each sample
rabund=function(sites){
  rabund=sites/SD_algae.df.3$total
  return(rabund)
}
SD_rabund.df=rabund(SD_algae.df.3[,6:67])
SD_rabund.df=cbind(SD_algae.df.3$taxon,SD_rabund.df)
#dcast and melt are the modern equivalents of transposition.  they work really well until they don't.
SD_rabund.df=dcast(melt(SD_rabund.df,id.vars="SD_algae.df.3$taxon"),variable~SD_algae.df.3$taxon)
SD_rabund.df=separate(SD_rabund.df,variable,into=c("junk","Date","Substrate"),sep=c(1,7))
SD_rabund.df=separate(SD_rabund.df,Substrate,into=c("Substrate","Rep"),sep=c(-1))
SD_rabund.df$junk=NULL
SD_rabund.df=droplevels(SD_rabund.df[-which(SD_rabund.df$Substrate=='CR'),])
SD_rabund.df$Date=as.factor(SD_rabund.df$Date)
SD_rabund.df$Date=revalue(SD_rabund.df$Date,c("170818"="3","170828"="10","170905"="21","170915"="31"))
SD_rabund.h=aggregate(SD_rabund.df[,4:17],list(SD_rabund.df$Date,SD_rabund.df$Substrate),mean)
SD_taxa=list(colnames(SD_rabund.h[,4:16])) 
SD_rabund.h=melt(SD_rabund.h,id.vars=c("Group.1","Group.2"))

#creates relative abundance subsets according to date
SD_rabund.1=filter(SD_rabund.df, Date==3)
SD_rabund.1$Date=NULL
SD_rabund.1=aggregate(SD_rabund.1[,3:31],list(SD_rabund.1$Substrate),mean)
SD_rabund.1=melt(SD_rabund.1, id="Group.1")
SD_rabund.1$Group.1=as.factor(SD_rabund.1$Group.1)

SD_rabund.2=filter(SD_rabund.df,Date==10)
SD_rabund.2$Date=NULL
SD_rabund.2=aggregate(SD_rabund.2[,3:31],list(SD_rabund.2$Substrate),mean)
SD_rabund.2=melt(SD_rabund.2, id="Group.1")

SD_rabund.3=filter(SD_rabund.df, Date==21)
SD_rabund.3$Date=NULL
SD_rabund.3=aggregate(SD_rabund.3[,3:31],list(SD_rabund.3$Substrate),mean)
SD_rabund.3=melt(SD_rabund.3, id="Group.1")

SD_rabund.4=filter(SD_rabund.df, Date==31)
SD_rabund.4$Date=NULL
SD_rabund.4=aggregate(SD_rabund.4[,3:31],list(SD_rabund.4$Substrate),mean)
SD_rabund.4=melt(SD_rabund.4, id="Group.1")

#select top 40% by substrate
SD_rabund.1.10=SD_rabund.1%>%
  group_by(Group.1)%>%
  arrange(Group.1,desc(value))%>%
  slice(seq(n()*.4)) #here is where you can select your cutoff level
SD_rabund.1.sigma=dcast(SD_rabund.1.10,Group.1~variable)
SD_rabund.1.sigma=t(SD_rabund.1.sigma)
write.table(SD_rabund.1.sigma,file="170818_rabund.txt",sep="\t")
SD_rabund.2.10=SD_rabund.2%>%
  group_by(Group.1)%>%
  arrange(Group.1,desc(value))%>%
  slice(seq(n()*.4))
SD_rabund.2.sigma=dcast(SD_rabund.2.10,Group.1~variable)
SD_rabund.2.sigma=t(SD_rabund.2.sigma)
write.table(SD_rabund.2.sigma,file="170828_rabund.txt",sep="\t")
SD_rabund.3.10=SD_rabund.3%>%
  group_by(Group.1)%>%
  arrange(Group.1,desc(value))%>%
  slice(seq(n()*.4))
SD_rabund.3.sigma=dcast(SD_rabund.3.10,Group.1~variable)
SD_rabund.3.sigma=t(SD_rabund.3.sigma)
write.table(SD_rabund.3.sigma,file="170905_rabund.txt",sep="\t")
SD_rabund.4.10=SD_rabund.4%>%
  group_by(Group.1)%>%
  arrange(Group.1,desc(value))%>%
  slice(seq(n()*.4))
SD_rabund.4.sigma=dcast(SD_rabund.4.10,Group.1~variable)
SD_rabund.4.sigma=t(SD_rabund.4.sigma)
write.table(SD_rabund.4.sigma,file="170915_rabund.txt",sep="\t")

#EnvFit
SD_env.df <- read.csv("SD_env.df.csv") #external data for vector fitting (lab assays, field data,etc)
SD_env.drops=c("8/15/2017",
            "8/22/2017",
            "8/25/2017",
            "9/1/2017",
            "9/8/2017",
            "9/12/2017") 
#beware of date format shifts in excel!!
#I am dropping dates from the env.df because I didn't get community measurements for each date
SD_env.df=droplevels(SD_env.df[!(SD_env.df$Date.Sampled %in% SD_env.drops),],drop=TRUE)
colnames(SD_env.df)=c("Date_Sampled",
                   "Substrate",
                   "NEP",
                   "GPP",
                   "Respiration",
                   "BG-ase",
                   "NAG-ase",
                   "P-ase",
                   "Chl.a",
                   "Biomass",
                   "Day of Incubation",
                   "SRP",
                   "NO3",
                   "NH4",
                   "doy",
                   "width",
                   "xdepth",
                   "xvelocity",
                   "xtemp",
                   "xspc",
                   "Q")
SD_env.c.drops=c("X",  #I don't necessarily want to use all the data here, so I specified some columns to ignore (summary columns, duplicates, etc)
              "Date.Sampled",
              "Sample",
              "time",
              "doy",
              "width",
              "xdepth",
              "xvelocity",
              "xtemp",
              "xspc",
              "Q")
SD_env.df.fit=SD_env.df[ , !(names(SD_env.df) %in% SD_env.c.drops)] #this is the code to remove the specified columns

SD_algae.df.sp=subset.data.frame(SD_algae.df, taxlevel==3,drop=TRUE) %>% 
  select(-taxlevel,-rankID,-daughterlevels,-total)
SD_algae.df.sp<-SD_algae.df.sp  %>% 
  gather(key="Sample",value="obs",-taxon) %>% 
  separate(Sample, into=c("Junk","Date","Substrate","Rep"),sep=c(1,7,-1))
SD_algae.df.sp$Junk=NULL
SD_algae.df.sp=droplevels(SD_algae.df.sp[-which(SD_algae.df.sp$Substrate=='CR'),])
SD_env.df$Date_Sampled=as.factor(SD_env.df$Date_Sampled)
SD_env.df$Date_Sampled=revalue(SD_env.df$Date_Sampled,c("8/18/2017"="3","8/29/2017"="10","9/5/2017"="21","9/15/2017"="31"))
SD_env.df$Substrate=factor(SD_env.df$Substrate,levels=c("H","So","Sh","F","T"),
                        labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))
SD_env.df.names=c("Date_Sampled","Substrate")
SD_env.df.fac=SD_env.df[,names(SD_env.df) %in% SD_env.df.names]
SD_algae.df.sps<-SD_algae.df.sp %>% 
  spread(key=taxon, value=obs) %>% 
  dplyr::select(-Date,-Substrate,-Rep)
#spread take a tall data table and makes it wide.  each "key" value becomes a column populated by the "value"

SD_alg.ord.dist=vegdist(SD_algae.df.sps) #get distances
SD_alg.ord.nmds=metaMDS(SD_alg.ord.dist) #run nmds
SD_alg.ORD=data.frame(MDS1=SD_alg.ord.nmds$points[,1],MDS2=SD_alg.ord.nmds$points[,2]) #exrtact nmds coordinates

#BREAK IS HERE
#SD_alg.ORD=cbind(SD_alg.ORD,SD_env.df.fac) #combine nmds coordinates with env data
#SD_alg.fit=envfit(SD_alg.ord.nmds,SD_env.df.fit,perm=999,na.rm=T,display="sites") #correlate you nmds coordintes with env data  Higher perm number means mroe sensitive p value
#SD_alg.fit.df=as.data.frame(SD_alg.fit$vectors$arrows*sqrt(SD_alg.fit$vectors$r))
#SD_alg.fit.df$species=rownames(SD_alg.fit.df)
#SD_alg.spp.scrs=as.data.frame(scores(SD_alg.fit,display="vectors"))
#SD_alg.spp.scrs=cbind(SD_alg.spp.scrs, Assay=rownames(SD_alg.spp.scrs),SD_alg.fit$vectors$pvals)
#SD_alg.spp.scrs=setnames(SD_alg.spp.scrs, old=("SD_alg.fit$vectors$pvals"), new=c("pvals"))
#SD_alg.spp.scrs=SD_alg.spp.scrs[ which (SD_alg.spp.scrs$pvals<0.05),] #select only significant covariates

