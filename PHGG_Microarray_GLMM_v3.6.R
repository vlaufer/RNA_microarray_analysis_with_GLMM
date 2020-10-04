# Licensed under MIT license
# For use in research publications please contact vlaufer@med.umich.edu

#######################################################################
####### Install required packages and troubleshoot installation #######
#######################################################################

#Set or reset libpath if needed:
my_lib_path<-'/data/user/vlaufer/R_libraries'
.libPaths(my_lib_path)
library('BiocManager')

################################################
###### Package Installation, if necessary ######
################################################
# standard
install.packages('BiocManager', lib=my_lib_path)
install.packages("maEndToEnd", lib=my_lib_path) 
install.packages('statmod', lib=my_lib_path)
install.packages('gridExtra', lib=my_lib_path)
install.packages('ggrepel', lib=my_lib_path)
install.packages('dplyr', lib=my_lib_path)   # MAEndtoEnd will fail to install if this isn't installed
install.packages("cli", lib=my_lib_path)     # MAEndtoEnd will fail to install if this isn't installed

# use bioconductor to install:
# may have to run then update then rerun:
BiocManager::install("maEndToEnd", lib=my_lib_path, dependencies=TRUE)
update.packages(lib.loc=my_lib_path)
BiocManager::install("maEndToEnd", lib=my_lib_path, dependencies=TRUE)

# load libraries
update.packages(lib.loc=my_lib_path)
library('ReactomePA')
library('arrayQualityMetrics')
library("huex10sttranscriptcluster.db"); library('hugene10sttranscriptcluster.db')
library('pd.hugene.1.0.st.v1'); library("pd.huex.1.0.st.v2") 
library("maEndToEnd")
library('statmod')
library('gridExtra')
library('ggrepel')
library('fgsea')
library('limma')
library('xlsx')

# Install packages and load libraries. For the maEndtoEnd workflow, required packages are:
package_list<-c("BiocManager", "Biobase", "biomaRt", "oligoClasses", "ArrayExpress", "edgeR", "pd.hugene.1.0.st.v1", "hugene10sttranscriptcluster.db", 
	"oligo", "arrayQualityMetrics", "limma", "topGO", "ReactomePA", "clusterProfiler", "gplots", "ggplot2", "ggrepel", "plotly", "geneplotter", 
	"RColorBrewer", "pheatmap", "dplyr", "tidyr", "stringr", "matrixStats", "genefilter", "openxlsx", "ReactomePA", "arrayQualityMetrics", 
	"huex10sttranscriptcluster.db", "hugene10sttranscriptcluster.db", "pd.hugene.1.0.st.v1", "pd.huex.1.0.st.v2", "maEndToEnd", "statmod", "gridExtra", 
	"ggrepel", "fgsea", "limma")

# Normally, user would run:
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0 ) { install.packages(new_packages) }
lapply(package_list, require, character.only = TRUE)
update.packages()  

install.packages("BiocManager"); library('BiocManager')
BiocManager::install("maEndToEnd"); library("maEndToEnd") 

# But maEndtoEnd throws errors during installation. So, try:
if ( "BiocManager" %in% installed.packages()[,"Package"] == TRUE ) {
	library('BiocManager')
} else { install.packages("BiocManager"); library('BiocManager') }

if ( "maEndToEnd" %in% installed.packages()[,"Package"] == TRUE ) {
	library("maEndToEnd") 
} else { BiocManager::install("maEndToEnd"); update.packages(); library('maEndToEnd') }
library('gdata')

install.packages('statmod')
install.packages("BiocManager"); library('BiocManager')
BiocManager::install("maEndToEnd")

# If this fails, try running the update packages command and re-running:
BiocManager::install("limma"); library('limma')


##################################################################
####### Initial Work to analyze  TMZ affy array data, uses #######
####### maEndToEnd and oligo package for affy array 2.0 ST #######
##################################################################

# Set dirs and initial variables
project_name<-"PHGG_metabolomics"
# base_dir<-"/Users/vincentlaufer/Desktop/Collaborations/PHGG_multiomics/phgg_multiomics"	# PC
base_dir<-"/data/scratch/vlaufer/shared/PHGG_omics_project/shared_code_and_output"			# cluster
script_dir=paste(base_dir, "Code", sep="/"); if (!file.exists(script_dir)) { dir.create(file.path(script_dir)) }
source(paste0(script_dir, "/Microarray_processing_functions.R"))

data_dir=paste(base_dir, "Data", sep="/"); if (!file.exists(data_dir)) { dir.create(file.path(data_dir)) }
meta_data_dir=paste(data_dir, "Meta_Data", sep="/"); if (!file.exists(meta_data_dir)) { dir.create(file.path(meta_data_dir)) }
CEL_dir=paste(data_dir, "CEL_files", sep="/"); if (!file.exists(CEL_dir)) { dir.create(file.path(CEL_dir)) }
R_object_dir<-paste(base_dir, "R_objects", sep="/"); if (!file.exists(R_object_dir)) { dir.create(file.path(R_object_dir)) }
figure_dir<-paste(base_dir, "Figures", sep="/"); if (!file.exists(figure_dir)) { dir.create(file.path(figure_dir)) }
results_dir=paste(base_dir, "Results", sep="/"); if (!file.exists(results_dir)) { dir.create(file.path(results_dir)) }

# QC Metrics Directory
QM_obj_path<-paste( R_object_dir, "QualMetrics", sep='/')
if (!file.exists(QM_obj_path)) { dir.create(file.path(QM_obj_path)) }
QM_figure_path<-paste( figure_dir, "QualMetrics", sep='/')
if (!file.exists(QM_figure_path)) { dir.create(file.path(QM_figure_path)) }

# set file names:
meta_data_fname<-'PHGG_Meta_Data_Master_31_Mar_v6.txt'
raw_data_object_path<-paste0( R_object_dir, "/", project_name, "_raw_expression_data.Rdata")
log2_raw_data_object_path<-paste0( R_object_dir, "/", project_name, "_log2_raw_expression_data.Rdata")
raw_PCA_data_object_path<-paste0( R_object_dir, "/", project_name, "_raw_PC_data.Rdata")
QM_obj_file<-paste(QM_obj_path, "QM_obj.Rdata", sep="/")
RLE_data_object_path<-paste0( R_object_dir, "/", project_name, "_RMA_expression_data.Rdata")
tmz_RLE_data_object_path<-paste0( R_object_dir, "/", project_name, "_tmz_RLE_expression_data.Rdata")
tmz_RMA_data_object_path<-paste0( R_object_dir, "/", project_name, "_tmz_RMA_expression_data.Rdata")
RMA_PCA_data_object_path<-paste0( R_object_dir, "/", project_name, "_RMA_PCA_results.Rdata")
RMA_data_norm_object_path<-paste0( R_object_dir, "/", project_name, "_RMA_norm_expression_data.Rdata")
GSEA_DEG_object_path<-paste0( R_object_dir, "/", project_name, "_GSEA_DEG_data.Rdata")
GLMM_model_object_path<-paste0( R_object_dir, "/", project_name, "_GLMM_Model.Rdata")
FE_model_1_object_path<-paste0( R_object_dir, "/", project_name, "_FE_Model_1.Rdata")
FE_model_2_object_path<-paste0( R_object_dir, "/", project_name, "_FE_Model_2.Rdata")
tmz_final_object_path<-paste0( R_object_dir, "/", project_name, "_tmz_final.Rdata")


# Set initial variables
setwd(base_dir)
# rm(list=ls())
alpha_level<-0.05
set.seed(151845)
options(scipen=999)

###################################################
###### Now, enact the MA end to end workflow ######
###################################################
# more information here: 
# https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html

########################################################################
###### Read in the microarray files and the sample key / metadata ######
########################################################################
SDRF_raw<-read.table(file = paste0(meta_data_dir, "/", meta_data_fname), sep="\t", header=TRUE)

# include only the .CEL files we wish to include from the directory: 
# SDRF_raw<-SDRF_raw[grepl("Pair", SDRF_raw$Patient_Number) == TRUE, ]


# Exclude microarrays found to have excessively high variance during QC (below) during pass 2:
names_to_remove<-c('X59T_1', 'X2609_3', 'X1066_1', 'X1066_2', 'X1066_3', 'X456_6', 'X2587_1', 'X14T_2')

SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X59T_1', ]
SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X2609_3', ]
SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X1066_1', ]
SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X1066_2', ]
SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X1066_3', ]
SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X456_6', ]
SDRF_raw<-SDRF_raw[SDRF_raw$UID != 'X2587_1', ]
SDRF_final<-SDRF_raw[SDRF_raw$UID != 'X14T_2', ]
rownames(SDRF_final) <- SDRF_final$Array
SDRF <- AnnotatedDataFrame(SDRF_final)
SDRF

# Get groupings of sample replicates
Pair_group_list<-list()
for (group_name in levels(SDRF_final$Patient_Number) ) { Pair_group_list[[ group_name ]]<-as.character(SDRF_final[SDRF_final$Patient_Number %in% group_name, ]$Array) }

Treatment_group_list<-list()
for (group_name in levels(SDRF_final$Treatment) ) { 
	SDRF_final$Treatment<-droplevels(SDRF_final$Treatment)
	Treatment_group_list[[ group_name ]]<-as.character(SDRF_final[SDRF_final$Treatment %in% group_name, ]$Array) 
}

# .CEL file names to be included in the analysis:
array.file.list<-paste(CEL_dir, SDRF_final$Array, sep="/")

####################################################################
###### Get uncorrected and RMA transformed expression objects ######
####################################################################
if ( file.exists( raw_data_object_path ) == TRUE ) { 
	load(file=raw_data_object_path); print("cel files previously loaded; retrieving stored data")
}

if ( file.exists( raw_data_object_path ) == FALSE ) { 
	raw_data <- oligo::read.celfiles(filenames = array.file.list, verbose = FALSE, phenoData = SDRF)
	save(raw_data, file = raw_data_object_path)
}

# get log2 normalized data.
if ( file.exists( log2_raw_data_object_path ) == TRUE ) { 
	load(file=log2_raw_data_object_path); print("cel files previously loaded; retrieving stored data")
}

if ( file.exists( log2_raw_data_object_path ) == FALSE ) { 
	log2_raw_data <- log2(Biobase::exprs(raw_data))
	save(log2_raw_data, file = log2_raw_data_object_path)
}

#############################################
###### Conduct PCA on uncorrected data ######
#############################################

if ( file.exists( raw_PCA_data_object_path ) == TRUE ) { 
	load(file=raw_PCA_data_object_path); print("PCA previously run; retrieving stored data")
}

if ( file.exists( raw_PCA_data_object_path ) == FALSE ) { 
	PCA_raw <- prcomp(t(log2_raw_data), scale. = FALSE)
	save(PCA_raw, file = raw_PCA_data_object_path)
}


###############################################
###### Make PC plots of uncorrected data ######
###############################################
vars_to_select<-c('Treatment_2', 'Passage')
variable_class<-c('factor', 'numeric')
raw_data_GG<-make_df_for_microarray_PCA( pData(raw_data), PCA_raw, vars_to_select, variable_class)

# plot uncorrected expression value PCA
main_title="PCA plot of the log2(intensity)"
vars_to_plot<-c('PC1', 'PC2', vars_to_select); raw_PC_plot_12<-make_PC_plot(raw_data_GG, PCA_raw, vars_to_plot, "Uncorrected expression values")
vars_to_plot<-c('PC3', 'PC4', vars_to_select); raw_PC_plot_34<-make_PC_plot(raw_data_GG, PCA_raw, vars_to_plot, "Uncorrected expression values")

# Inspect PC plots visually 
PCA_fname<-paste0( figure_dir, "/", project_name, "_PC_plot_raw.pdf")
PC_plot_list<-list(raw_PC_plot_12, raw_PC_plot_34)
grid.arrange( grobs= PC_plot_list, nrow=2, ncol=1)
ggsave(PCA_fname, width = 30, height = 25, units = "cm", device='pdf', arrangeGrob(grobs = PC_plot_list, ncol = 1) )



#################################################
###### Get QC Metrics for uncorrected data ######
#################################################
if ( file.exists( QM_obj_file ) == TRUE ) { load(file=QM_obj_file); print("Quality Metrics previously generated; retrieving stored data") }

if ( file.exists( QM_obj_file ) == FALSE ) { 
	variables_of_interest<-c("Treatment_2", "Patient_Number")
	QM_obj<-arrayQualityMetrics(expressionset = raw_data, force = TRUE, do.logtransform = TRUE, intgroup = variables_of_interest, outdir = QM_figure_path )
	save(QM_obj, file=QM_obj_file)
}

##################################################################
###### Make boxplots of uncorrected data to compare samples ######
##################################################################

######
## Make simple boxplot of intensities
sample_boxplot_raw<-oligo::boxplot(raw_data, target = "core", main = "Boxplot of log2-intensitites for the raw data",xlabels=raw_data$UID)
raw_boxplot_name<-paste0( figure_dir, "/", project_name, "_intensity_boxplot_raw.pdf")
pdf(file=raw_boxplot_name, width=14, height=12)
oligo::boxplot(raw_data, target = "core", main = "Boxplot of log2-intensitites for the raw data",)
dev.off()


##################################################
###### Make Heatamaps from sample distances ######
##################################################

# make distance HM
distance_HM_raw<-make_sample_distance_HM(raw_data, 'Treatment', 'Patient_Number') 
distance_HM_file_raw<-paste0( figure_dir, "/", project_name, "_distance_HM_raw.pdf")
pdf(distance_HM_file_raw, width=14, height=12)
distance_HM_raw
dev.off()
###### Based on results we can consider dropping 10-26-3 and 10-26-2 and 10-23-1 ######
#######################################################################################




##################################################
###### Make Histograms of intensity values #######
##################################################
# Make histograms of probe intensity to identify possible problems/strategies
lograw_data_medians <- rowVars(log2_raw_data)
hist_res<-hist(lograw_data_medians, 100, col = "cornsilk1", freq = FALSE, main = "Histogram of Median of log2(intensity)", border = "antiquewhite4", xlab = "Median intensities")

intensity_hist_file<-paste0( figure_dir, "/", project_name, "_median_intensity_hist_raw.pdf")
#man_threshold <- 4 #arbitrary manual threshold for low med. intensity, can alter on reanalysis

# make intensity histogram for uncorrected data
man_threshold<-4
pdf(intensity_hist_file)
hist_res
abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()

########################################################################################################################################################################

##############################################################################
###### Apply initial RMA algorithm without normalization to expr values ######
##############################################################################
# This will background-correct, normalize and summarize the data
if ( file.exists( tmz_RLE_data_object_path ) == TRUE ) { 
	load(file=tmz_RLE_data_object_path); print("RMA transform previously applied; retrieving stored data")
}

if ( file.exists( tmz_RLE_data_object_path ) == FALSE ) { 
	tmz_rma_for_RLE <- oligo::rma(raw_data,target="core", normalize=FALSE )
	save(tmz_rma_for_RLE, file = tmz_RLE_data_object_path)
	RLE_data <- Biobase::exprs(tmz_rma_for_RLE)
	save(RLE_data, file = tmz_RLE_data_object_path)
}

##################################################################
###### Make RLE boxplots of RMA data find outlying samples #######
##################################################################
######
## now show the cumulative variance of all the samples in a given sample:
RLE_mat<-as.matrix(Biobase::exprs(tmz_rma_for_RLE))
row_medians_assayData <- Biobase::rowMedians(RLE_mat)
RLE_sweep <- sweep(Biobase::exprs(tmz_rma_for_RLE), 1, row_medians_assayData)
colnames(RLE_sweep)<-pData(tmz_rma_for_RLE)[['UID']]
RLE_data_gathered <- tidyr::gather(as.data.frame(RLE_sweep), UID, log2_expression_deviation) # May need to change "patient_array" and "log2_expression_deviation"

# now make RLE boxplot
RLE_plot<-ggplot(RLE_data_gathered, aes(UID, log2_expression_deviation))
RLE_plot<-RLE_plot + geom_boxplot(outlier.shape = NA)
final_RLE_plot<-RLE_plot + theme(axis.text.x = element_text(colour = "aquamarine4", angle = 60, size = 6.5, hjust = 1 , face = "bold"))
# final_RLE_plot
# dev.off()
a<-oligo::boxplot(tmz_rma_for_RLE, target = "core", main = "Boxplot of log2-intensitites for the raw data", xlabels=tmz_rma_for_RLE$UID)
sample_boxplot_raw<-oligo::boxplot(a, target = "core", main = "Boxplot of log2-intensitites for the raw data", xlabels=a$UID)

RLE_boxplot_name<-paste0( figure_dir, "/", project_name, "_RLE_boxplot_raw_2.pdf")
pdf(file=RLE_boxplot_name)
final_RLE_plot
dev.off()
###### This plot suggests that 10-23-1 and 10-26-2 and 10-26-3 may be problematic #######
#########################################################################################

########################################################
###### Now confirm sample issues using a heatmap #######
########################################################

# make distance HM
distance_HM_raw<-make_sample_distance_HM(tmz_rma_for_RLE, 'Treatment', 'Patient_Number') 
distance_HM_file_raw<-paste0( figure_dir, "/", project_name, "_distance_HM_raw_2.pdf")
pdf(distance_HM_file_raw, width=16, height=14)
distance_HM_raw
dev.off()
###### This plot confirms that 59T-1 and 14T-2 are both problematic. From here we will remove then rerun. #######
#################################################################################################################

if ( file.exists( tmz_RMA_data_object_path ) == TRUE ) { 
	load(file=tmz_RMA_data_object_path); print("RMA transform previously applied; retrieving stored data")
}

if ( file.exists( tmz_RMA_data_object_path ) == FALSE ) { 
	RMA_data <- oligo::rma(raw_data,target="core" )
	save(RMA_data, file = tmz_RMA_data_object_path)
	RMA_data_norm <- Biobase::exprs(RMA_data)
	save(RMA_data_norm, file = RMA_data_norm_object_path)
}


####################################################
###### Make PC plots of RMA transformed  data ######
####################################################

if ( file.exists( RMA_PCA_data_object_path ) == TRUE ) { 
	load(file=RMA_PCA_data_object_path); print("PCA previously run; retrieving stored data")
}

if ( file.exists( RMA_PCA_data_object_path ) == FALSE ) { 
	PCA_RMA <- prcomp(t(RMA_data_norm), scale = FALSE)
	save(PCA_RMA, file = RMA_PCA_data_object_path)
}



############################################################################################################
##### Further process the transformed data to avoid paying multiple testing penalty during final step: #####
############################################################################################################
### Make a high pass filter to remove low counts.
man_threshold <- 4 #arbitrary manual threshold for low med. intensity, can alter on reanalysis
no_of_samples <- table( pData(tmz_rma)[[ 'Treatment' ]] )
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(tmz_rma), 1, function(x){ sum(x > man_threshold) >= samples_cutoff } )
table(idx_man_threshold)
# FALSE  TRUE 
# 2681 19330 
manfiltered <- subset(tmz_rma, idx_man_threshold)


# ### Next, make a low pass filter for within-sample variances
# rma_data<-Biobase::exprs(manfiltered)
# colnames(rma_data)<-colnames(Biobase::exprs(raw_data))
# rma_data_variance_p1s <- rowVars(rma_data[, Treatment_group_list[[ "Treated" ]] ])
# rma_data_variance_p2s <- rowVars(rma_data[, Treatment_group_list[[ "Untreated" ]] ])
# VarSums<-rowSums( cbind (rma_data_variance_p1s, rma_data_variance_p2s), na.rm=TRUE)
# hist(VarSums)
# var_filtered<-VarSums<2.5
# manfiltered<-subset(manfiltered, var_filtered)


##########################################################################################
########### Annotate the RMA transformed data and process into semi-final form ###########
##########################################################################################
anno <- AnnotationDbi::select(huex10sttranscriptcluster.db, keys = (featureNames(manfiltered)), columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
anno <- subset(anno, !is.na(SYMBOL))

# group annotations by gene symbol
anno_grouped <- group_by(anno, PROBEID)
anno_summarized <- dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL)); head(anno_summarized)

# exclude probes that map non-uniquely
anno_filtered <- filter(anno_summarized, no_of_matches > 1); head(anno_filtered)
probe_stats <- anno_filtered; nrow(probe_stats)
ids_to_exlude <- (featureNames(manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)
# FALSE  TRUE
# 17070  2260

# add PCs to semi-final object in case you need to control as covariate



tmz_final <- subset(manfiltered, !ids_to_exlude) # get rid of multiple mappings
fData(tmz_final)$PROBEID <- rownames(fData(tmz_final)) # make an fdata table with probeIDs 
fData(tmz_final) <- left_join(fData(tmz_final), anno) # use anno (like converter) to add SYMBOL
# at this point, still 17070 features, 61 samples 

rownames(fData(tmz_final)) <- fData(tmz_final)$PROBEID 
tmz_final_2<-subset(tmz_final, complete.cases(fData(tmz_final)))
validObject(tmz_final_2)
# now down to 13820 probes ... lose about 4000 genes due to missing data

tmz_final<-tmz_final_2

pData(tmz_final)$PC1<-PCA_RMA$x[,1]
pData(tmz_final)$PC2<-PCA_RMA$x[,2]
pData(tmz_final)$PC3<-PCA_RMA$x[,3]
pData(tmz_final)$PC4<-PCA_RMA$x[,4]


save(tmz_final, file=tmz_final_object_path)

#########################################################################################################
#################### Begin DEG analysis using the normalized and filtered array data ####################
#########################################################################################################
# will use simplest model possible to start, analyzing for DEG between treated and resistant only. 

## get gene signatures of glioma genes ready from the literature:
glioma_signatures<-list()
glioma_signatures[["Diaz"]]<-c("AGT", "EGFR", "CHI3L1", "SOD2", "CCL2", "IGFBPL1", "MBP", "CPE", "OLFM1", "MCF", "PACSIN1") # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4170681/
glioma_signatures[["Arima"]]<-c("TP53", "EGFR", "FAS", "TARP", "CMYC", "NEDD4L", "RB1", "DLK1", "SLITRK6", "PTEN", "PTPRD", "PDGFB", 
	"CCND2", "CDK4", "CDKN2A", "CDKN2B", "DDX1", "KCNRG", "KRAS", "MDM4", "MET", "MYCN", "PDGFRA", "PIK3C2B") # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639942/
glioma_signatures[["Mackay"]]<-c("H3F3A", "BRAF", "TP53", "ATRX", "PDGFRA", "CDKN2A", "CDKN2B", "NF1", "RB1", "PTEN", 
"PIK3CA", "CDK4", "KIT", "BCOR", "KDR", "EGFR", "TERG", "MDM2", "TSC2") # from https://www.cell.com/cancer-cell/pdf/S1535-6108(18)30175-2.pdf

all_glioma_genes<-unique(c(glioma_signatures[["Diaz"]], glioma_signatures[["Arima"]], glioma_signatures[["Mackay"]]))

# prepare data for LM:
Passage <- as.numeric(pData(tmz_final)$Passage)
Treatment <- pData(tmz_final)$Treatment_2
Treatment<-relevel(Treatment, ref='Untreated')
Match_Group <- pData(tmz_final)$Patient_Number
#Match_Group_2 <- paste(pData(tmz_final)$Patient_Number, pData(tmz_final)$Treatment, sep="_")
PC1 <- pData(tmz_final)$PC1; PC2 <- pData(tmz_final)$PC2
PC3 <- pData(tmz_final)$PC3; PC4 <- pData(tmz_final)$PC4

####### For dealing with paired data:
# https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#122_a_linear_model_for_the_data
# this also refers the reader to limma 9.4 (page 43), here: https://bioconductor.org/packages/3.11/bioc/vignettes/limma/inst/doc/usersguide.pdf
# which provides useful info. as well.
# You need to create a blocking variable. For us that will be Match_Group

#####
# So, lets add our blocking variable (to control for pair in paired sample design)
# and lets investigate our PCs to see if we should consider adding them as well:
variable_name<-'TreatmentTreated' # this is the variable that we want to analyze. 
variables_of_interest<-c('Treatment', 'Passage', 'Patient_Number')
design2 <- model.matrix(~Passage + Treatment) # --> with the blocking variable, many are highly significant.
design3 <- model.matrix(~Passage + Match_Group + Treatment) # 

find_PC_correlations(tmz_final, 'PC1', variables_of_interest)
find_PC_correlations(tmz_final, 'PC2', variables_of_interest)
find_PC_correlations(tmz_final, 'PC3', variables_of_interest)
find_PC_correlations(tmz_final, 'PC4', variables_of_interest)
find_PC_correlations(tmz_final, 'Passage', variables_of_interest)

# the correlation between PC and ___ is highest; we will try modeling PC ___ as a covariate like so: 
design4 <- model.matrix(~PC1 + Passage + Match_Group + Treatment)
design5 <- model.matrix(~PC2 + Passage + Match_Group + Treatment)
design6 <- model.matrix(~PC3 + Passage + Match_Group + Treatment) #
design7 <- model.matrix(~PC4 + Passage + Match_Group + Treatment)
design8 <- model.matrix(~PC1 + PC2 + Passage + Match_Group + Treatment)
design9 <- model.matrix(~PC1 + PC3 + Passage + Match_Group + Treatment)
design10<- model.matrix(~PC1 + PC2 + PC3 + Passage + Match_Group + Treatment)
design11 <- model.matrix(~PC3 + Match_Group + Treatment)

#model_1_results<-test_model(tmz_final, design1, glioma_signatures, variable_name)
model_2_results<-test_model(tmz_final, design2, glioma_signatures, variable_name)
model_3_results<-test_model(tmz_final, design3, glioma_signatures, variable_name)
model_4_results<-test_model(tmz_final, design4, glioma_signatures, variable_name)
model_5_results<-test_model(tmz_final, design5, glioma_signatures, variable_name)
model_6_results<-test_model(tmz_final, design6, glioma_signatures, variable_name) # 
model_7_results<-test_model(tmz_final, design7, glioma_signatures, variable_name)
model_8_results<-test_model(tmz_final, design8, glioma_signatures, variable_name)
model_9_results<-test_model(tmz_final, design9, glioma_signatures, variable_name)
model_10_results<-test_model(tmz_final, design10, glioma_signatures, variable_name)
model_11_results<-test_model(tmz_final, design11, glioma_signatures, variable_name)

model_2_results[['results']]
model_3_results[['results']] # 
model_4_results[['results']]
model_5_results[['results']]
model_6_results[['results']]  # ***********************
model_7_results[['results']]
model_8_results[['results']]
model_9_results[['results']]
model_10_results[['results']]
model_11_results[['results']]  # ***********************
#

writeLines(capture.output(sessionInfo()), paste0( R_object_dir, "/", project_name, "_R_sessionInfo.txt") )



save(model_6_results, file=FE_model_1_object_path)
save(model_11_results, file=FE_model_2_object_path)




FE_model_selected<-model_6_results

write.table( FE_model_selected[['fit']]$t, file =paste0(results_dir, "/", project_name, "_FE_model_", 't_statistics',".txt"  ) ) 
write.table( FE_model_selected[['TopHits']], file =paste0(results_dir, "/", project_name, "_FE_model_", 'top_hits',".txt"  ) ) 
write.table( FE_model_selected[['results']], file =paste0(results_dir, "/", project_name, "_FE_model_", 'num_results',".txt"  ) ) 
write.table( FE_model_selected[['design']], file =paste0(results_dir, "/", project_name, "_FE_model_", 'design',".txt"  ) ) 

################################################################################################
#################################### For implementing GLMM  ####################################
################################################################################################
# https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf  section 17.3.6 page 111.
Passage <- as.numeric(pData(tmz_final)$Passage)
Treatment <- pData(tmz_final)$Treatment_2
Treatment<-relevel(Treatment, ref='Untreated')
Match_Group <- pData(tmz_final)$Patient_Number
glmm_design_1 <- model.matrix(~0 + Passage + Treatment)  
glmm_design_2 <- model.matrix(~0 + PC1 + PC2 + PC3 + PC4 + Passage + Treatment)
glmm_design_3 <- model.matrix(~PC3 + Passage + Treatment)                # same as the best model, which was design6 (from FE model).
dupcor <- duplicateCorrelation(tmz_final, glmm_design_1, block= Match_Group )
dupcor$consensus.correlation # degree of w/in sample variation to be absorbed.

# Now we look for differentially expressed genes. We make all possible pairwise comparisons between
# the epithelial cell types, allowing for the correlation within donors:
glmm_res_list1<-list()
fit1 <- lmFit(tmz_final, glmm_design_1, block=Match_Group, correlation=dupcor$consensus.correlation)
glmm_res_list1[["fit"]] <- eBayes(fit1, trend=TRUE, robust=TRUE)
summary(decideTests(fit1, method="global"))
glmm_res_list1[["TopHits"]]<-topTable(glmm_res_list1[["fit"]], coef=variable_name, adjust.method = "BH", number=Inf)
glmm_res_list1[["results"]]<-summary(decideTests(fit1, method="global"))
glmm_res_list1[["design"]]<-fit1$design
glmm_res_list1[["results"]]


glmm_res_list3<-list()
dupcor <- duplicateCorrelation(tmz_final, glmm_design_3, block= Match_Group )
dupcor$consensus.correlation # degree of w/in sample variation to be absorbed.
fit3 <- lmFit(tmz_final, glmm_design_3, block=Match_Group, correlation=dupcor$consensus.correlation)
summary(decideTests(fit3, method="global"))
glmm_res_list3[["fit"]] <- eBayes(fit3, trend=TRUE, robust=TRUE)
glmm_res_list3[["TopHits"]]<-topTable(glmm_res_list3[["fit"]], coef=variable_name, adjust.method = "BH", number=Inf)
glmm_res_list3[["results"]]<-summary(decideTests(fit3, method="global"))
glmm_res_list3[["design"]]<-fit3$design
glmm_res_list3[["results"]]

#
GLMM_selected1<-glmm_res_list1
GLMM_selected3<-glmm_res_list3

#
GLMM_selected<-glmm_res_list3
save(GLMM_selected, file=GLMM_model_object_path)

#
write.table( GLMM_selected[['fit']][['t']], file =paste0(results_dir, "/", project_name, "_GLMM_", 't_statistics',".txt"  ), quote=FALSE, row.names=FALSE, sep='\t' )
write.table( GLMM_selected[['TopHits']], file =paste0(results_dir, "/", project_name, "_GLMM_", 'top_hits',".txt"  ), quote=FALSE, row.names=FALSE, sep='\t' )
write.table( GLMM_selected[['results']], file =paste0(results_dir, "/", project_name, "_GLMM_", 'num_results',".txt"  ), quote=FALSE, row.names=FALSE, sep='\t' )
write.table( GLMM_selected[['design']], file =paste0(results_dir, "/", project_name, "_GLMM_", 'design',".txt"  ), quote=FALSE, row.names=FALSE, sep='\t' )

#
save(GLMM_selected, file=selected_model_object_path)
save(tmz_final, file= )

########### GLMM appears to outcompete the FE model.

# if you want to make specific contrasts (e.g. due to 3 variables.)
# contrasts <- makeContrasts( PC1 + Passage + Resistant - Sensitive, levels=Treatment)
#  fit2 <- contrasts.fit(fit, contrasts)
# fit2 <- eBayes(fit2, trend=TRUE)
#fit_2 <- lmFit(tmz_final, glmm_design_2, block=Match_Group, correlation=dupcor$consensus.correlation)
#summary(decideTests(fit_2, method="global"))


################################################
########### Model Comparison Metrics ###########
################################################
####
# now lets compare our best models: 



###############################################
########### Select a model and Save ###########
###############################################
GLMM_selected1<-glmm_res_list1
GLMM_selected3<-glmm_res_list3
GSEA_model1<-GLMM_selected1
# Chose a model and print DEG results to excel file
DEG_results_name<-paste0( results_dir, "/", project_name, "_DEG_results.txt")
write.table(x = topTable(GSEA_model$fit, number=100, coef=variable_name, adjust.method = "BH"), file = DEG_results_name, 
  sep="\t", quote=F,row.names = T, col.names = T)


