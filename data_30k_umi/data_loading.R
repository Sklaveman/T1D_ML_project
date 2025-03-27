library(data.table)
library(stringr)

clonotype_table_preprocessing <- function(dt, file){
  if('uniqueMoleculeCount' %in% colnames(dt)){
    dt <- dt[,.(aaSeqCDR3, uniqueMoleculeCount)]
    dt <- dt[,sum(uniqueMoleculeCount), by=aaSeqCDR3]
    setnames(dt, c('aaSeqCDR3', 'V1'), c('amino_acid', 'templates'))
  } else {
    dt <- dt[,.(aaSeqCDR3, readCount)]
    dt <- dt[,sum(readCount), by=aaSeqCDR3]
    setnames(dt, c('aaSeqCDR3', 'V1'), c('amino_acid', 'templates'))
  }
  write.table(dt, paste0('~/Diabetes/T1D_ML_project/data_30k_umi/repertoires/',
                         file),
              sep = '\t', quote = F, row.names = F)
}

ngsik <- fread('/home/klupyr/Diabetes/Bulk_repertoires/all_data_from_ngsik.tsv')
files <- list.files('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/T1D_batch1')
ngsik <- ngsik[paste0(ngsik$name, '.tsv') %in% files,]

metadata <- c()

print('old bulk data')
pb <- txtProgressBar(min = 0, max = length(ngsik$name), initial = 0)
step <- 0
for(name_t in ngsik$name){
  dt <- fread(paste0('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/T1D_batch1/',
                     name_t, '.tsv'))
  status <- ngsik[ngsik$name == name_t, Type1Diabetes]
  name_t <- paste0(name_t, '.tsv')
  if(sum(dt$uniqueMoleculeCount) == 30000){
    dt <- clonotype_table_preprocessing(dt, paste0(name_t, '.tsv'))
    metadata <- rbind(metadata, c(paste0(name_t, '.tsv'), status, 'T1D_batch1'))
  }
  step <- step + 1
  setTxtProgressBar(pb, step)
}

############################# new bulk data  ###################################

files <- list.files('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/T1D_batch2/')
files <- files[files != 'NTC_PS_TRB_TRB.tsv']

print('new bulk data')
pb <- txtProgressBar(min = 0, max = length(files), initial = 0)
step <- 0

for(name in files){
  dt <- fread(paste0('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/T1D_batch2/',
                     name))
  print(sum(dt$uniqueMoleculeCount))
  if(sum(dt$uniqueMoleculeCount) == 30000){
    dt <- clonotype_table_preprocessing(dt, name)
    metadata <- rbind(metadata, c(name, 'Yes', 'T1D_batch2'))
  }
  step <- step + 1
  setTxtProgressBar(pb, step)
}

################################### batch 3 ####################################

files <- list.files('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/T1D_batch3/')
files <- files[files != 'NTC_PS_TRB.tsv']

pb <- txtProgressBar(min = 0, max = length(files), initial = 0)
step <- 0

for(name in files){
  dt <- fread(paste0('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/T1D_batch3/',
                     name))
  if(sum(dt$uniqueMoleculeCount) == 30000){
    dt <- clonotype_table_preprocessing(dt, name)
    metadata <- rbind(metadata, c(name, 'Yes', 'T1D_batch3'))
  }
  step <- step + 1
  setTxtProgressBar(pb, step)
}

################################################################################
########################## healthy individuals data ############################
################################################################################

########################## rosati new alignment  ###############################

files <- list.files('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/rosati/')

print('rosati')
pb <- txtProgressBar(min = 0, max = length(files), initial = 0)
step <- 0

for(file in files){
  dt <- fread(paste0('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/rosati/', file))
  if(sum(dt$readCount) == 30000){
    dt <- clonotype_table_preprocessing(dt, file)
    metadata <- rbind(metadata, c(file, 'No', 'rosati'))
  }
  step <- step + 1
  setTxtProgressBar(pb, step)
}

############################### aging dataset ##################################

files <- list.files('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/aging/')

print('aging')
pb <- txtProgressBar(min = 0, max = length(files), initial = 0)
step <- 0

for(file in files){
  dt <- fread(paste0('/projects/cdr3_common/user/klupyr/T1D_and_healthy_downsample_30k/tsv_files/aging/', file))
  if(sum(dt$uniqueMoleculeCount) == 30000){
    dt <- clonotype_table_preprocessing(dt, file)
    metadata <- rbind(metadata, c(file, 'No', 'aging'))
  }
  step <- step + 1
  setTxtProgressBar(pb, step)
}

files <- list.files('~/Diabetes/T1D_ML_project/data_30k_umi/repertoires/')
metadata <- as.data.frame(metadata)
metadata <-metadata[595:nrow(metadata),]

files[!(files %in% metadata$V1)]
setnames(netadata, c('V1', 'V2', 'V3'), c('name', 'Type1Diabetes', 'batch'))
write.table(metadata, file = '~/Diabetes/T1D_ML_project/data_30k_umi/metadata.tsv',
            sep = '\t', quote = F, row.names = )
setwd('~/Diabetes/T1D_ML_project/')
for(file in files){
  system(paste0('git add ~/Diabetes/T1D_ML_project/data_30k_umi/repertoires/', file))
}