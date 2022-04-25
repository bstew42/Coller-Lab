# primary metastatic samples with tumor purity greater than 40%
pure_prim_met <- read.csv('pure_prim_met.csv')

# primary non-metastatic samples with tumor purity greater than 40%
pure_prim_no_met <- read.csv('pure_no_met_prim.csv')

# data set containing both the primary metastatic and primary-non metastatic samples joint by the gene column. Columns with .x are metastatic, .y is non-metastatic
pure_prim <- dplyr::inner_join(pure_prim_no_met, pure_prim_met, by = 'Gene')

# removing percent signs from frequency column
pure_prim$Freq.x <- gsub('%', '', pure_prim$Freq.x)
pure_prim$Freq.y <- gsub('%', '', pure_prim$Freq.y)

# removing frequencies labeled as >0.1.
pure_prim <- pure_prim[pure_prim$Freq.x >= 0.1 & pure_prim$Freq.y >= 0.1,]

# converting frequencies to numeric to work with as numeric values instead of strings
pure_prim$Freq.x <- as.numeric(pure_prim$Freq.x)
pure_prim$Freq.y <- as.numeric(pure_prim$Freq.y)

# run this to find list of genes with the frequency of metastatic genes greter than 50% difference in frequency than non-metastatic
pure_prim[(pure_prim$Freq.x - pure_prim$Freq.y) / pure_prim$Freq.x > 0.5, c('Gene', 'Freq.x', 'Freq.y')]

# run this to find list of genes with the frequency of metastatic genes less than 50% difference in frequency than non-metastatic
pure_prim[(pure_prim$Freq.x - pure_prim$Freq.y) / pure_prim$Freq.x < -0.5, c('Gene', 'Freq.x', 'Freq.y')]

# just for interest, shows genes with greater than 5% frequency in both the metastatic and non-metastatic
pure_greater_than_five <- pure_prim[pure_prim$Freq.x > 5 & pure_prim$Freq.y > 5, ]

# same as line 22 but with more meaningful values
pure_greater_than_five[(pure_greater_than_five$Freq.x - pure_greater_than_five$Freq.y) / pure_greater_than_five$Freq.x > 0.5, c('Gene', 'Freq.x', 'Freq.y')]

# same as line 25 but more meaningful values
pure_greater_than_five[(pure_greater_than_five$Freq.x - pure_greater_than_five$Freq.y) / pure_greater_than_five$Freq.x < -0.5, c('Gene', 'Freq.x', 'Freq.y')]

#histogram showing is frequency differences between genes that a metastatic vs not, doesn't give much info
hist(pure_prim$Freq.x, breaks = 200, col = "blue", density = 20)
hist(pure_prim$Freq.y, breaks = 200, col = "red", density = 30, add = TRUE)



