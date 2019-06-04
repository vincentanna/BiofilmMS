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
SD_otu.df=list_to_df(otu) #conevrt list to data frame wiht custom function above
SD_otu.df$value=as.numeric(SD_otu.df$value)
SD_otu.df=separate(SD_otu.df,name,into=c("junk","Date","Substrate"),sep=c(1,7))
SD_otu.df=separate(SD_otu.df,Substrate,into=c("Substrate","Rep"),sep=c(-1))
SD_otu.df$junk=NULL
SD_otu.df=droplevels(SD_otu.df[-which(SD_otu.df$Substrate=='CR'),])
SD_otu.df$Date=as.factor(SD_otu.df$Date)
SD_otu.df$Date=revalue(SD_otu.df$Date,c("170818"="3","170828"="10","170905"="21","170915"="31"))
SD_out.aov=SD_otu.df%>% #tidy pipes!  Run anovs by Date on otu richenss
  group_by(Date)%>%
  do(tidy(aov(value~Substrate,data=.)))
SD_out.aov
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

