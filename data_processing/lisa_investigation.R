load_bed_data<-function(sample_dta,target_chromosomes){
  
  data = lapply(with(sampleTable, setNames(path, name)),
                function(x){
                  if(endsWith(x = x,'.gz')){
                    tmp = fread(cmd=paste("gunzip -c", x), showProgress=FALSE,
                                col.names=c("seqnames", "start", "end", "score"),
                                colClasses=c("character", "numeric", "numeric", "numeric"))
                  }else{
                    tmp = fread(x, showProgress=FALSE,
                                col.names=c("seqnames", "start", "end",'name', "score",'strand'),
                                colClasses=c("character", "numeric", "numeric",'character', "numeric",'character'))
                    tmp<-tmp[,c("seqnames", "start", "end", "score")]
                  }
                  
                  #chr23 and chr24 map to X and Y
                  tmp[seqnames == "chr23", seqnames := "chrX"]
                  tmp[seqnames == "chr24", seqnames := "chrY"]
                  
                  #Add 'chr' if it isn't there
                  if(startsWith(x=tmp$seqnames[1],'chr')==FALSE){
                    tmp$seqnames<-paste0('chr',tmp$seqnames)
                  }
                  setkeyv(tmp, c("seqnames", "start", "end"))
                  
                  tmp = tmp[seqnames%in%chromosomes,]
                  return(as.data.frame(tmp))
                })
  
  return(data)
}



library(data.table)
sampleTable = data.table(name = c("P1c-SMUG1","P1a-SMUG1","P1b-SMUG1","P2a-SMUG1","P2b-SMUG1","P2c-SMUG1",
                                  "KO1a-RNPII","KO1b-RNPII","KO2a-RNPII","KO2b-RNPII",
                                  "P1a","P1b","P1c","P2a","P2b","P2c",
                                  "KO1a","KO1b","KO1c","KO2a","KO2b","KO2c",
                                  # "P1a","P1b","P1c","P2a","P2b","P2c",
                                  # "KO1a","KO1b","KO1c","KO2a","KO2b","KO2c",
                                  "P1a","P1b","P1c","P2a","P2b","P2c","P3a","P3b","P3c",
                                  "KO1a","KO1b","KO1c","KO2a","KO2b","KO2c"
),
path_loc=c(rep('/media/yohanl/Expansion/seq-processing/CutTag_SMUG1/results/',10),
            rep('/media/yohanl/Expansion/seq-processing/Lisa_5hmu/results/',12),
            # rep('/media/yohanl/Expansion/seq-processing/Lisa_ATAC/ATACseq_first_run/results/',12),
            rep('/media/yohanl/Expansion/seq-processing/Lisa_ATAC/ATACseq_second_run/results/',15))
)



path_vect<-c()
peak_vect<-c()
for(i in 1:nrow(sampleTable)){
  name<-sampleTable$name[i]
  path_loc<-sampleTable$path_loc[i]
  new_path<-paste0(path_loc,name,'/peaks_bw_bed/',name,'_sorted.bed')
  path_vect<-c(path_vect,new_path)
  if(any(grepl("\\.narrowPeak$", list.files(paste0(path_loc,name,'/peaks_bw_bed/'))))==TRUE){
    new_path<-paste0(path_loc,name,'/peaks_bw_bed/',name,'_sorted_peaks.narrowPeak')
  }else{
    new_path<-''
  }
  peak_vect<-c(peak_vect,new_path)
}
sampleTable$path<-path_vect
sampleTable$peaks_path<-peak_vect



#Filter for standard chromosomes
chromosomes = paste0('chr',c(1:22, "X",'Y'))
data=load_bed_data(sample_dta=sampleTable,target_chromosomes=chromosomes)