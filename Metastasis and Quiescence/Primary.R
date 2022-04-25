# data from primary samples with no metastasis
prim_no_met <- read.csv('prim_no_met.csv')

# data from primary samples with metastasis
prim_met <- read.csv('prim_met.csv')

# clincial data to work around with
cplay <- read.csv('clinical_data.csv')

# removing rows with frequency values less than 0.1
prim_met <- prim_met[prim_met$Freq >= 0.1, ]
prim_no_met <- prim_met[prim_no_met$Freq >= 0.1, ]

# removing percent signs from frequency values
prim_met$Freq<-gsub("%$","",prim_met$Freq)
prim_no_met$Freq<-gsub("%$","",prim_no_met$Freq)

# converting frequency values to numeric
prim_met$Freq <- as.numeric(prim_met$Freq)
prim_no_met$Freq <- as.numeric(prim_no_met$Freq)

# joining both sample data sets into one, .x is for no metastasis, .y is for metastasis
prim <- dplyr::inner_join(prim_no_met, prim_met, by = 'Gene')

# showing the genes with greated difference between the two sets
head(prim[order((prim$Freq.x - prim$Freq.y) / prim$Freq.x), 'Gene'])
tail(prim[order((prim$Freq.x - prim$Freq.y) / prim$Freq.x), 'Gene'])
prim[prim$Gene %in% head(prim[order((prim$Freq.x - prim$Freq.y) / prim$Freq.x), 'Gene']), c('Gene', "Freq.x", 'Freq.y')]
prim[prim$Gene %in% tail(prim[order((prim$Freq.x - prim$Freq.y) / prim$Freq.x), 'Gene']), c('Gene', "Freq.x", 'Freq.y')]

#genes with frequency greater than 5 in no met or met
freq_greter_than_five <- prim[prim$Freq.x > 5 | prim$Freq.y > 5, c('Gene', "Freq.x", "Freq.y")]

#genes with at least 50% more freq in no met than met
freq_greter_than_five[((freq_greter_than_five$Freq.x - freq_greter_than_five$Freq.y) / freq_greter_than_five$Freq.x > 0.5),]

#genes with at least 50% more in met than no met
freq_greter_than_five[((freq_greter_than_five$Freq.x - freq_greter_than_five$Freq.y) / freq_greter_than_five$Freq.x < -0.5),]

#finding graph of tumor purity to decide where to cut it off in future research
tumor_purity <- clinical$Tumor.Purity[!is.na(clinical$Tumor.Purity)]
hist(tumor_purity, breaks = 10)
mean(tumor_purity)
median(tumor_purity)


