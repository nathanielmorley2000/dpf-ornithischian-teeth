# load libraries
library("plyr")
library("dplyr") 
library("ggord")
library("ggplot2")
library("MASS")
library("openxlsx")
library("pcaMethods")
library("purrr")
library("readxl") 
library("tibble")
library("tidyr")
library("vegan")


# install BiocManager and pcaMethods - keep commented out unless needed
#if (!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager") 
#  }
#BiocManager::install("pcaMethods")


# import data
crownAll<-read_excel("ceratopsidIncluded/dpfornithishcians_ceratopsidsIncluded.xlsx", sheet="Raw Data")
crownAll<-crownAll[1:90,1:17]
colnames(crownAll)<-c("Family", "Subfamily", "Genus", "Specimen Number", "Tooth Type", "CBL", "CH","M-CL", "AL",  "NL", "MCL", "DCL", "MDH", "MDL", "DDH", "DDL", "CA")


# alter factor levels to match those used for referencing in the subsequent code
crownAll$Family<-revalue(crownAll$Family, c("Ceratopsidae"="C", "Ankylosauria"="A", "Leptoceratopsidae"="Le", "Pachycephalosauria"="Pa", "Thescelosauridae"="Th"))
crownAll$Subfamily<-revalue(crownAll$Subfamily, c("Centrosaurinae"="Ce", "Chasmosaurinae"="Ch", "Ankylosauridae"="Anky", "Nodosauridae"="No","Leptoceratopsidae"="Le", "Pachycephalosauridae"="Pa", "Orodrominae"="Or"))



# Data Preparation =============================================================

## Define functions ============================================================

### Perform BPCA ===============================================================

### takes data frame of certain dimensions(e.g. cranial_hadro[,8:19]) and imputes, retrieves data
##data frame must have data in columns 8 to 19 while columns 1 to 7 are documentation
BPCA_data <- function(data,range1,range2,range3){
  added_lin<-data[,range1] + 0.001
  log_lin<-log10(added_lin)
  log_all<-cbind(log_lin, data[,range3])
  
  # range represents the log-transformed data columns through numbers in square brackets with your first numerical column and your last numerical column, respectifully. 
  # use correlation matrix for PCA because SSI is a unitless measurement
  est_data<-pcaMethods::pca(log_all,method="bpca",centre=FALSE, cor=TRUE)
  
  est_matrix<-est_data@completeObs # This will retrieve the estimated matrix from the previous function.
  fin<-cbind(data[,range2],est_matrix)
  #est_matrix # This will display your new estimated matrix.  You can  copy and paste it out of R if you want.
  
  return(fin)}


## Define imputation function ==================================================

# impute_classy separates into subsets based on family and then runs BPCA imputation of a subset of a provided dataset of cranial variables (cranial) 
## using the range of the data columns (range1) and the range for the specimen information (range2)
## it then combines the imputed subsets and returns a single imputed dataset
#### NOTE: the data columns will be log-transformed prior to imputation
impute_class<-function(crown,range1,range2,range3){
  # range1 represents the range for the data columns
  # range2 represents the label and other specimen information columns
  # crown is the dataset
  ##################################################
  ### create individual data sets for each subfamily for imputation
  ank<-subset(crown, Family=="A")
  thesc<-crown
  pachy<-subset(crown, Family=="Pa")
  lepto<-subset(crown, Family=="Le")
  cer<-subset(crown, Family=="C")
  ########################################################
  
  ##################################################

  
  ##### run imputation on individual data sets
  finCranAnk<-BPCA_data(ank,range1, range2, range3)
  finCranCer<-BPCA_data(cer,range1, range2, range3)
  finCranLepto<-BPCA_data(lepto,range1, range2, range3)
  finCranThesc<-BPCA_data(thesc,range1, range2, range3)
  finCranPachy<-BPCA_data(pachy,range1, range2, range3)
  finCranThesc <- subset(finCranThesc, Family =="Th")
  
  cranData<-rbind(finCranAnk,finCranCer,finCranLepto,finCranThesc,finCranPachy)### combine the above imputed subsets into a single dataset 
  
  return(cranData) }

###############################################################################################################################################
#################################################### run the imputation #######################################################################
###############################################################################################################################################

# impute missing values and obtain categorical information that will be needed later for analysis; note that all values obtained will be log10 transformed
# note that because ankylosaurs, ceratopsids, leptos, thescelosaurids and pachys are imputed separately, it is not necessary to differentiate between types of analyses at this point
cran<-impute_class(crownAll, 6:16, 1:5, 17)
# note that the returned values are log-transformed

###OLS###
#Set up the function
res <- function(cran, voi, p) {
  out <- data.frame(matrix(ncol = 10, nrow = 90))
  out[,1:5]<-cran[,1:5]
  out[,6]<-cran[,17]
  for (i in voi) {
    testRes<-lm(formula = cran[,i] ~ cran[,p], data = cran)
    out[,i] <- testRes$resid
  }
  colnames(out)<-c("Family", "Subfamily", "Genus", "Specimen Number", "Tooth Type", "CA", "CH","M.CL", "AL",  "NL", "MCL", "DCL", "MDH", "MDL", "DDH", "DDL")
  return(out)
}

#Call the function
crownAdjusted<-res(cran, 7:16, 6)

###Canonical Variates Analysis###
#Create CVA function
ldaData <- function(crownAdjusted, taxonRank) {
  # rename column for fun
  findforme=which(colnames(crownAdjusted)==taxonRank)
  colnames(crownAdjusted)[findforme]="myRank"
  #Create and plot CVA
  linear<-lda(myRank~CH+M.CL+AL+NL+MCL+DCL+MDH+MDL+DDH+DDL, data=crownAdjusted)
  scaling=linear$scaling
  proportion_of_trace=as.vector(linear$svd^2/sum(linear$svd^2))
  myplot=ggord(linear, crownAdjusted$myRank, ellipse = FALSE, hull = TRUE, exp=0.06)
  #Confusion Matrix
  p <- predict(linear, crownAdjusted, method="debiased")$class
  tab <- table(Predicted = p, Actual = crownAdjusted$myRank)
  mytab=sum(diag(tab))/sum(tab)
  #Give ordination and confusion matrix
  giveme=list(myplot, mytab, tab, linear, scaling, proportion_of_trace)
  names(giveme)=c("myplot","perctab", "mytab", "lda", "coefficients", "proportion.of.trace")
  return(giveme)
}

### Call the function ####
# family-level CVA
ldaFamily<-ldaData(crownAdjusted, "Family")
ldaFamily$mytab
ldaFamily$perctab
ldaFamily$lda
ldaFamily$coefficients
ldaFamily$proportion.of.trace


# subfamily-level CVA
ldaSubfamily<-ldaData(crownAdjusted, "Subfamily")
ldaSubfamily$mytab
ldaSubfamily$perctab
ldaSubfamily$lda
ldaSubfamily$coefficients
ldaSubfamily$proportion.of.trace


# store results in a dataframe for export
max_ln <- max(c(length(ldaFamily$proportion.of.trace), length(ldaSubfamily$proportion.of.trace))) 
myLDAproportions <- data.frame(Family = c(ldaFamily$proportion.of.trace,rep(NA, max_ln - length(ldaFamily$proportion.of.trace))), 
                      Subfamily = c(ldaSubfamily$proportion.of.trace,rep(NA, max_ln - length(ldaSubfamily$proportion.of.trace)))) %>%
  t() %>%
  as.data.frame() %>%
  round(., digits = 4) %>%
  rownames_to_column() %>%
  dplyr::rename(CV1 = V1, CV2 = V2, CV3 = V3, CV4 = V4, CV5 = V5, CV6 = V6)
myLDAFam.components=ldaFamily$coefficients%>%as.data.frame() %>% 
  round(., digits = 4) %>%
  rownames_to_column()%>%
  dplyr::rename(CV1 = LD1, CV2 = LD2, CV3 = LD3, CV4 = LD4)
myLDAsubFam.components=ldaSubfamily$coefficients%>%as.data.frame() %>% 
  round(., digits = 4) %>%  
  rownames_to_column() %>%
  dplyr::rename(CV1 = LD1, CV2 = LD2, CV3 = LD3, CV4 = LD4, CV5 = LD5, CV6 = LD6)
myLDAFam.tab=ldaFamily$mytab%>%as.data.frame()%>%spread(.,Actual, Freq)
myLDAsubFam.tab=ldaSubfamily$mytab%>%as.data.frame()%>%spread(.,Actual, Freq)
perctab=data.frame(class=c("Family","Subfamily"), perc=c(ldaFamily$perctab, ldaSubfamily$perctab) )%>% # input percentages into one data.frame for export
  mutate(perc=round(perc,3))# round to three digs after decimal

# export information from data.frame into excel
wrkbk <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb=wrkbk, "crnadj" ) # crownadjusted
openxlsx::writeData(wrkbk, 1, crownAdjusted)

openxlsx::addWorksheet(wb=wrkbk, "TraceProportion" )
openxlsx::writeData(wrkbk, 2, myLDAproportions)

openxlsx::addWorksheet(wb=wrkbk, "FamComponents" )
openxlsx::writeData(wrkbk, 3, myLDAFam.components)

openxlsx::addWorksheet(wb=wrkbk, "SubfamComponents" ) 
openxlsx::writeData(wrkbk, 4, myLDAsubFam.components)

openxlsx::addWorksheet(wb=wrkbk, "FamMatrix" )
openxlsx::writeData(wrkbk, 5, myLDAFam.tab)

openxlsx::addWorksheet(wb=wrkbk, "SubfamMatrix" )
openxlsx::writeData(wrkbk, 6, myLDAsubFam.tab)

openxlsx::addWorksheet(wb=wrkbk, "percTabs" )
openxlsx::writeData(wrkbk, 7, perctab)

openxlsx::saveWorkbook(wrkbk, "ceratopsidIncluded/Results/toothy_LDA_results.xlsx")
############################################################################

########### new function to produce plots ###########
#colscheme is the color scheme desired
#y is the data from the automatic lda plot produced
lda_plottr=function(y, colscheme, shapescheme, groupname){
  ############################################
  # produce convex hulls
  find_hull<-function(df) df[chull(df$one,df$two),]
  combhull=data.frame()
  i=1
  stop=y$Groups%>%unique()%>%length()
  myGroups=y$Groups%>%unique()
  
  while(i<stop+1){
    h1<-find_hull(y%>%subset(Groups==myGroups[[i]])) # repeat for all grouping variables
    combhull=rbind(combhull, h1)
    i=i+1
  }
  myplot=
  ggplot()+
    geom_point(data=y, aes(x=one, y=two, color=Groups, shape=Groups), size=5)+
    geom_polygon(data=combhull, aes(x=one, y=two, color=Groups, fill=Groups), alpha=0.1, show.legend = F)+
    scale_color_manual(name = "Taxon", labels = groupname, values=colscheme)+
    scale_fill_manual(name = "Taxon", labels = groupname, values=colscheme)+
    scale_shape_manual(name = "Taxon", labels = groupname, values=shapescheme)+
    xlab("CV1")+
    ylab("CV2")+
    theme_classic()+
    theme(axis.text=element_text(color="black"),
          axis.title.x=element_text(size=18, vjust=-0.5, color="black"),
          axis.text.x=element_text(size=16, color="black"),
          axis.title.y=element_text(size=18, angle=90, vjust=3, color="black"),
          axis.text.y=element_text(size=16, angle=90, vjust=2, color="black"),
          axis.ticks=element_line(color="black"),
          axis.ticks.length=unit(2, "mm"),
          panel.background = element_blank(),
          panel.border=element_rect(fill=NA, color="black", size=1, linetype="solid"),
          plot.margin = unit(c(3, 2, 4, 4), "mm"),
          legend.title = element_text(size = 24),
          legend.text = element_text(size=18),
          legend.key.size = unit(1.0, 'cm'))
  return(myplot)}

################################################### Family ############################################

# CVA plot
x=ldaFamily$lda
x=ldaFamily$lda%>%.$scaling%>%data.frame()
y=ldaFamily$myplot%>%.$data

colscheme= set_names(x=c("#ff1000", "#00ff36", "#014421", "#0006ff", "#fc00ff"), 
                     nm=c("A", "C", "Le", "Pa", "Th"))
shapescheme=set_names(x=c(15,17,4,8,7), nm=c("A", "C", "Le", "Pa", "Th"))
groupname = set_names(x=c("Ankylosauria", "Ceratopsidae", "Leptoceratopsidae", "Pachycephalosauridae", "Thescelosauridae"), 
                      nm=c("A", "C", "Le", "Pa", "Th") )
lda_plottr(y, colscheme, shapescheme, groupname)

# save CVA plot as .svg file
ggsave(filename= "ceratopsidIncluded/Results/familyLDA.svg", 
  width = 700,
  height = 500,
  units = "mm", 
  dpi= 1000,
  scale = 0.6)

# loadings plot for CV1 and CV2
variables <- myLDAFam.components$rowname
  
loadingsFam <- myLDAFam.components %>%
  subset(., select = c(rowname, CV1, CV2)) %>%
  pivot_longer(!rowname, names_to = "component", values_to = "loading") 

ggplot(loadingsFam, aes(x = factor(rowname, variables), y = loading, fill = component)) +
  geom_bar(stat="identity", width = 0.8, position = position_dodge(width = 0.8)) + 
  labs(y = "Loading",
       x = "Variable") +
  guides(fill=guide_legend(title="Axis")) +
  theme_classic(base_size = 33) +
  theme(axis.title.x=element_text(size=35, vjust=-0.5, color="black"),
        axis.title.y=element_text(size=35, angle=90, vjust=3, color="black"),
        legend.title = element_text(size = 38),
        legend.text = element_text(size=35))

# save loadings plot as .svg file   
ggsave(filename= "ceratopsidIncluded/Results/familyloadings.svg",
       width = 800,
       height = 600,
       units = "mm", 
       dpi= 1000,
       scale = 0.6)

################################################### Subfamily ############################################

# CVA plot
x=ldaSubfamily$lda
x=ldaSubfamily$lda%>%.$scaling%>%data.frame()
y=ldaSubfamily$myplot%>%.$data


colscheme=set_names(x=c("#ff1000", "#ff6c00", "#00ff36", "#964B00", "#014421", "#0006ff", "#fc00ff"), 
                    nm=c("Anky", "No", "Ce", "Ch", "Le", "Pa", "Or"))
shapescheme=set_names(x=c(15,16,17,3,4,8,7), nm=c("Anky", "No", "Ce", "Ch", "Le", "Pa", "Or"))
groupname = set_names(x=c("Ankylosauridae", "Nodosauridae", "Centrosaurinae", "Chasmosaurinae", "Leptoceratopsidae", "Pachycephalosauridae", "Orodrominae"), 
                      nm=c("Anky", "No", "Ce", "Ch", "Le", "Pa", "Or") )
lda_plottr(y, colscheme, shapescheme, groupname)

# save CVA plot as .svg file
ggsave(filename= "ceratopsidIncluded/Results/subfamilyLDA.svg",
        width = 700,
        height = 500,
        units = "mm", 
        dpi= 1000,
        scale = 0.6)

# loadings plot for CV1 and CV2
loadingsSubfam <- myLDAsubFam.components %>%
  subset(., select = c(rowname, CV1, CV2)) %>%
  pivot_longer(!rowname, names_to = "component", values_to = "loading") 

ggplot(loadingsSubfam, aes(x = factor(rowname, variables), y = loading, fill = component)) +
  geom_bar(stat="identity", width = 0.8, position = position_dodge(width = 0.8)) + 
  labs(y = "Loading",
       x = "Variable") +
  guides(fill=guide_legend(title="Axis")) +
  theme_classic(base_size = 33) +
  theme(axis.title.x=element_text(size=35, vjust=-0.5, color="black"),
        axis.title.y=element_text(size=35, angle=90, vjust=3, color="black"),
        legend.title = element_text(size = 38),
        legend.text = element_text(size=35))

# save loadings plot as .svg file   
ggsave(filename= "ceratopsidIncluded/Results/subfamilyloadings.svg" ,
       width = 800,
       height = 600,
       units = "mm", 
       dpi= 1000,
       scale = 0.6)


###############################################################################################################################################
###############################################################################################################
################################## conduct PERMANOVA ###########################################################
###############################################################################################################

THE.STATS<-function(
    crownAdjusted,
    col.to.group.by, # the name of the column containing the grouping information
    lastColName,
    firstColName, # provide the name of the column containing the first measurements
    num.perm, # number of iterations
    return.everything # whether want everything to be returned, significant or not (yes or no)
){
  

  #########################################################################################################
  # preparation of data for analysis
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
  my.data=crownAdjusted%>%subset(Family!='Th') # remove thescelosaurids from dataset (n = 1 won't work with PERMANOVA)
  colnames(my.data)[ which(colnames(my.data)==col.to.group.by) ] <- 'group' # rename column of grouping variable
  my.data=my.data%>%mutate(group = as.factor(group) )
  
  group = my.data$group
 
  # find the column index for the first column containing measurements
  first.val.col= which(colnames(my.data)==firstColName)
  last.val.col= which(colnames(my.data)==lastColName)
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
  
  ################################### Run Test ########################################
  # omnibus PERMANOVA
  Y = cbind(my.data[, first.val.col:last.val.col])
  np.omnibus = adonis2(Y ~ group, permutations = num.perm, method = "mahalanobis", by = NULL, data = my.data)
  
  # create the pairwise function
  pairwise.adonis <- function(group, Y, sim.method, p.adjust.m, data){
    
    # create table
    co = combn(unique(as.character(group)), 2)
    pairs = c()
    F.Model =c()
    R2 = c()
    p.value = c()

    # run pairwise analysis
    for(elem in 1:ncol(co)){
      x1 = vegdist(Y[group %in% c(co[1,elem], co[2,elem]),], method=sim.method)
      
      # fill out table as we go
      ad = adonis2(x1 ~ group[group %in% c(co[1,elem],co[2,elem])]);
      pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
      F.Model =c(as.numeric(F.Model), as.numeric(ad[1,4])) %>% round(digits = 2);
      R2 = c(as.numeric(R2), as.numeric(ad[1,3])) %>% round(digits = 2);
      p.value = c(as.numeric(p.value), as.numeric(ad[1,5]) %>% round(digits = 3))
    }

    #correct p-values  
    p.adjusted = p.adjust(p.value,method=p.adjust.m) %>% round(., digits = 3)
    
    # create dataframe for results
    pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
    return(pairw.res)
  } 
  
  # call the pairwise function
  np.pair.test = pairwise.adonis(group, Y, sim.method = "mahalanobis", p.adjust.m = "holm", data = my.data)
  
  ######################################################################################################
  # prepare the results for return
  if(return.everything=="no"){
    if(np.omnibus$`Pr(>F)`[1]>0.05){np.pair.return="not significant"}
    if(np.omnibus$`Pr(>F)`[1]<=0.05){np.pair.return=np.pair.test}
  } # if only want significant stuff returned
  
  if(return.everything=="yes"){    
    np.pair.return=np.pair.test
  } # if want everything returned
  
  to.return=list("np.omnibus"=np.omnibus, "np.pair.test"=np.pair.return)
  
  return(to.return)  
}

########################################################################################################
########################### Call function and save results #############################################
########################################################################################################

# call PERMANOVA function for family
stres.family=THE.STATS(crownAdjusted,
                       col.to.group.by='Family', 
               lastColName= "DDL", 
               firstColName="CH", 
               num.perm=10000, 
               return.everything="no")

# call PERMANOVa function for subfamily
stres.subfamily=THE.STATS(crownAdjusted,
                       col.to.group.by='Subfamily', 
                       lastColName= "DDL", 
                       firstColName="CH", 
                       num.perm=10000, 
                       return.everything="no")

# open new excel sheets for PERMANOVA results
wrkbk2 <- openxlsx::createWorkbook()

# save results
openxlsx::addWorksheet(wb=wrkbk2, "Family Omnibus" )
openxlsx::writeData(wrkbk2, 1, stres.family$np.omnibus)

openxlsx::addWorksheet(wb=wrkbk2, "Family Pairwise" )
openxlsx::writeData(wrkbk2, 2, stres.family$np.pair.test)

openxlsx::addWorksheet(wb=wrkbk2, "Subfamily Omnibus" )
openxlsx::writeData(wrkbk2, 3, stres.subfamily$np.omnibus)

openxlsx::addWorksheet(wb=wrkbk2, "Subfamily Pairwise" )
openxlsx::writeData(wrkbk2, 4, stres.subfamily$np.pair.test)

# close and save excel sheets for PERMANOVA results
openxlsx::saveWorkbook(wrkbk2, "ceratopsidIncluded/Results/toothy_PERMANOVA_results.xlsx")