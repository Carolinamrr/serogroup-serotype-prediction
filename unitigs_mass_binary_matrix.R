#script creating binary matrix with unitigs as predictor

path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#loading libraries
library(stringr)

#Load data
id_links <- read.delim(paste0(path, "/datasets/strain_id_links.txt"))
meta <- read.csv(paste0(path, "/datasets/mass_metadata.csv"))
data_path <- paste0(path, "/datasets/mass_unitigs.txt.gz.")
zz <- gzfile(data_path, "rt")
data <- read.csv(zz, header = F)

#Change # to _ in id_links
id_links$Fastq <- gsub("#", "_", id_links$Fastq)

#merge Fastq with Accession + Serotype
id_sero <- merge(id_links, meta[, c("Serotype", "Sequence.Reads.Accession")], 
                 by.x = "Accession", by.y = "Sequence.Reads.Accession", 
                 all = TRUE)
id_sero <- id_sero[(id_sero$Fastq %in% unique(data$Fastq)),]
id_sero <- id_sero[order(id_sero$Fastq),]
rownames(id_sero) <- c(1:length(id_sero$Serotype))

#Separate data into unitigs and string of fastq ids
data <- str_split_fixed(data[,1], " ", 3)
data <- data[, -2]
data <- as.data.frame(data)
colnames(data) <- c("unitigs", "Fastq")
data$Fastq <- gsub(" ", "", data$Fastq)
data$Fastq <- gsub(":([0-9]+).*", "", data$Fastq)
data <- data[order(data$Fastq),]

gc() #free memory

# Create row and column indix for sparse matrix
rows <- c()
columns <- c()

for(i in 1:dim(data)[1]){
  fastq_sequences <- unique(unlist(strsplit(as.character(data[i, "Fastq"]),
                                            split = " ")))
  row <- rownames(id_sero)[(id_sero$Fastq %in% fastq_sequences)]
  rows <- c(rows, row)
  columns <- c(columns, rep(i, length(row_Fastq)))
}

saveRDS(rows, paste0(path, "/results/mass_rows.rds"))
saveRDS(columns, paste0(path, "/results/mass_columns.rds"))


