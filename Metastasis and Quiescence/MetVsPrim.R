clinical <- data.frame(read.csv('clinical_data.csv'))
clinical

met_mut <- read.csv('Metastatic_Mutated.csv')
met_struc <- read.csv('Metastatic_Structural.csv')
met_cna <- read.csv('Metastatic_CNA.csv')

prim_mut <- read.csv('Primary_Mutated.csv')
prim_struc <- read.csv('Primary_Structural.csv')
prim_cna <- read.csv('Primary_CNA.csv')

quiescent <- read.csv('Q_Data.csv')


sum(table(clinical$Patient.ID) > 1)

library(UpSetR)


length(met_mut$Gene)
sum(met_mut$Gene %in% prim_mut$Gene)

length(prim_mut$Gene)
sum(prim_mut$Gene %in% met_mut$Gene)

length(met_struc$Gene)
sum(met_struc$Gene %in% prim_struc$Gene)

length(prim_struc$Gene)
sum(prim_struc$Gene %in% met_struc$Gene)

met_mut$Gene
quiescent$Gene

sum(met_mut$Gene %in% quiescent$Gene) / length(met_mut$Gene)
sum(prim_mut$Gene %in% quiescent$Gene) / length(prim_mut$Gene)

sum(met_struc$Gene %in% quiescent$Gene) / length(met_struc$Gene)
sum(prim_struc$Gene %in% quiescent$Gene) / length(prim_struc$Gene)

sum(met_cna$Gene %in% quiescent$Gene) / length(met_cna$Gene)
sum(prim_cna$Gene %in% quiescent$Gene) / length(prim_cna$Gene)


#List containing quiescnet, metastatic(mutation and structural), and primary(mutaation and structural)
L1 <- list( q = quiescent$Gene,
      mm = met_mut$Gene,
      pm = prim_mut$Gene,
      ms = met_struc$Gene,
      ps = prim_struc$Gene
)
#upset plot from L1
upset(fromList(L1), order.by = 'freq')

#List containing quiescent, metastatic mutation, and primary mutation
L2 <- list(
  q = quiescent$Gene,
  mm = met_mut$Gene,
  pm = prim_mut$Gene
)
#upset from L2
upset(fromList(L2), order.by = 'freq')

#list containing quiescent, metastatic structural,and primary structural
L3 <- list(
  q = quiescent$Gene,
  ms = met_struc$Gene,
  ps = prim_struc$Gene
)
#upset from L3
upset(fromList(L3), order.by = 'freq')

#list containing quiescent, metastatic cna, and primary cna
L4 <- list(
  q = quiescent$Gene,
  mc = met_cna$Gene,
  pc = prim_cna$Gene
)
#upset from L4
upset(fromList(L4), order.by = 'freq')


#Now lets try to find mutations/structural/CNA that have a high frequency in metastatic and not primary or vice versa

both_mut <- met_mut$Gene[met_mut$Gene %in% prim_mut$Gene]

met_mut[met_mut$Gene == 'TP53', 'Freq']
prim_mut[prim_mut$Gene == 'TP53', 'Freq']

m <- numeric(0)
p <- numeric(0)

for (i in both_mut){
  f <- c(f, met_mut[met_mut$Gene == i, 'Freq'])
  p <- c(p, prim_mut[prim_mut$Gene == i, 'Freq'])
}
f
p


met_mut[met_mut[, 'Freq'] == '11.80%',]
prim_mut[prim_mut[, 'Freq'] == '8.90%',]

names(clinical)
sum(clinical$Number.of.Samples.Per.Patient > 1)

both_mut <- dplyr::inner_join(met_mut, prim_mut, by = 'Gene')
dif_freq_mut <- both_mut$Freq.x - both_mut$Freq.y
