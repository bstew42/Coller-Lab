Dormancy Analysis (CNA)
================

### Defining Data Sets

Chloe’s Dormancy Data

``` r
dormant <- read.csv("Dormant_Data.csv")
```

cna Genes from
<https://www.cbioportal.org/study/summary?id=msk_met_2021>

``` r
cna <- read.csv("CNA.csv")
```

Changing the column name for gene to “Gene” in Chloe’s dormancy data and
innerjoining the two sets by the Gene column

``` r
colnames(dormant)[4] <- "Gene"
data <- dplyr::inner_join(dormant, cna, by = "Gene")
```

### Quick numbers regarding genes the two sets share

Number of genes in Chloe’s data

``` r
length(dormant$Gene)
```

    ## [1] 5878

Number of genes in cbioportal data

``` r
length(cna$Gene)
```

    ## [1] 928

Number of Genes shared between the two

``` r
length(data$Gene)
```

    ## [1] 366

data will now have genes from the innerjoined set that have data from at
least 4 studies in Chloe’s data

``` r
data <- data[data$cell_line_counts >= 4, ]
```

## Upregulated Genes

Upregulated genes will be defined as genes with at least 4 cell line
counts and 1 or less downregulated count.

``` r
upreg <- data[data$downreg_counts <= 1, ]
```

The amount of upregulated genes fitting our description above

``` r
length(upreg[, 'Gene'])
```

    ## [1] 14

The upregulated genes

``` r
upreg$Gene
```

    ##  [1] "CDKN1A" "CDKN1A" "HLA-A"  "HLA-A"  "BCL2L1" "BCL2L1" "MAPK3"  "MAPK3" 
    ##  [9] "SESN2"  "SESN2"  "PNRC1"  "PNRC1"  "SMAD3"  "SMAD3"

Frequencies of upregulated genes from cbioportal data

``` r
rbind(upreg$Gene, upreg$Freq)
```

    ##      [,1]     [,2]     [,3]    [,4]    [,5]     [,6]     [,7]    [,8]   
    ## [1,] "CDKN1A" "CDKN1A" "HLA-A" "HLA-A" "BCL2L1" "BCL2L1" "MAPK3" "MAPK3"
    ## [2,] "0.20%"  "<0.1%"  "<0.1%" "0.10%" "1.50%"  "<0.1%"  "<0.1%" "0.30%"
    ##      [,9]    [,10]   [,11]   [,12]   [,13]   [,14]  
    ## [1,] "SESN2" "SESN2" "PNRC1" "PNRC1" "SMAD3" "SMAD3"
    ## [2,] "<0.1%" "0.20%" "0.20%" "<0.1%" "<0.1%" "0.20%"

BCL2L1 and MAPK3 have the highest frequencies of CNA, so I will explore
those further

### BCL2L1

``` r
BCL2L1 <- upreg[upreg$Gene == "BCL2L1" & upreg$Freq == "1.50%", ]
```

Type of CNA

``` r
BCL2L1$CNA
```

    ## [1] "AMP"

Mean of BCL2L1 cell line counts

``` r
mean(as.numeric(BCL2L1[which(!is.na(BCL2L1[5:18])) + 4]))
```

    ## [1] 0.7583011

Median of BCL2L1 cell line counts

``` r
median(as.numeric(BCL2L1[which(!is.na(BCL2L1[5:18])) + 4]))
```

    ## [1] 0.5866225

### MAPK3

``` r
MAPK3 <- upreg[upreg$Gene == "MAPK3" & upreg$Freq == "0.30%", ]
```

Type of CNA

``` r
MAPK3$CNA
```

    ## [1] "AMP"

Mean of MAPK3 cell line counts

``` r
mean(as.numeric(MAPK3[which(!is.na(MAPK3[5:18])) + 4]))
```

    ## [1] 0.3859559

Median of MAPK3 cell line counts

``` r
median(as.numeric(MAPK3[which(!is.na(MAPK3[5:18])) + 4]))
```

    ## [1] 0.5434011

## Downregulated Genes

Downregulated genes will be defined as genes with at least 4 cell line
counts and 1 or less upregulated count.

``` r
downreg <- data[data$upreg_counts <= 1, ]
```

Amount of downregulated genes fitting our description above

``` r
length(downreg[, 'Gene'])
```

    ## [1] 46

The downregulated genes

``` r
downreg$Gene
```

    ##  [1] "EIF1AX" "EIF1AX" "NPM1"   "NPM1"   "EZH2"   "EZH2"   "DNMT1"  "DNMT1" 
    ##  [9] "CDK4"   "CDK4"   "NUF2"   "CHEK1"  "CHEK1"  "RAD51D" "RAD51D" "BLM"   
    ## [17] "BLM"    "TOP1"   "TOP1"   "BRCA1"  "BRCA1"  "POLD1"  "POLD1"  "RHOA"  
    ## [25] "RHOA"   "XPO1"   "XPO1"   "MSH2"   "MSH2"   "BRCA2"  "BRCA2"  "RB1"   
    ## [33] "RB1"    "EIF4E"  "EIF4E"  "ELOC"   "SRSF2"  "SRSF2"  "PBRM1"  "PBRM1" 
    ## [41] "RAD21"  "RAD21"  "FANCA"  "FANCA"  "XRCC2"  "XRCC2"

Frequencies of downregulated genes

``` r
rbind(downreg$Gene, downreg$Freq)
```

    ##      [,1]     [,2]     [,3]    [,4]    [,5]    [,6]    [,7]    [,8]    [,9]   
    ## [1,] "EIF1AX" "EIF1AX" "NPM1"  "NPM1"  "EZH2"  "EZH2"  "DNMT1" "DNMT1" "CDK4" 
    ## [2,] "0.20%"  "0.50%"  "0.10%" "<0.1%" "<0.1%" "0.20%" "0.40%" "<0.1%" "2.10%"
    ##      [,10]   [,11]   [,12]   [,13]   [,14]    [,15]    [,16]   [,17]   [,18]  
    ## [1,] "CDK4"  "NUF2"  "CHEK1" "CHEK1" "RAD51D" "RAD51D" "BLM"   "BLM"   "TOP1" 
    ## [2,] "<0.1%" "0.50%" "<0.1%" "0.10%" "0.10%"  "0.10%"  "<0.1%" "0.20%" "<0.1%"
    ##      [,19]   [,20]   [,21]   [,22]   [,23]   [,24]   [,25]   [,26]   [,27]  
    ## [1,] "TOP1"  "BRCA1" "BRCA1" "POLD1" "POLD1" "RHOA"  "RHOA"  "XPO1"  "XPO1" 
    ## [2,] "0.70%" "<0.1%" "<0.1%" "0.10%" "0.10%" "<0.1%" "<0.1%" "0.10%" "<0.1%"
    ##      [,28]   [,29]   [,30]   [,31]   [,32]   [,33]   [,34]   [,35]   [,36]  
    ## [1,] "MSH2"  "MSH2"  "BRCA2" "BRCA2" "RB1"   "RB1"   "EIF4E" "EIF4E" "ELOC" 
    ## [2,] "<0.1%" "<0.1%" "0.40%" "0.30%" "0.20%" "1.80%" "<0.1%" "<0.1%" "1.40%"
    ##      [,37]   [,38]   [,39]   [,40]   [,41]   [,42]   [,43]   [,44]   [,45]  
    ## [1,] "SRSF2" "SRSF2" "PBRM1" "PBRM1" "RAD21" "RAD21" "FANCA" "FANCA" "XRCC2"
    ## [2,] "0.20%" "<0.1%" "0.20%" "<0.1%" "2.30%" "<0.1%" "<0.1%" "0.40%" "<0.1%"
    ##      [,46]  
    ## [1,] "XRCC2"
    ## [2,] "<0.1%"

RAD21, ELOC, RB1, and CDK4 have the highest CNA so I will examine those
closer.

### RAD21

``` r
RAD21 <- downreg[downreg$Gene == "RAD21" & downreg$Freq == "2.30%", ]
```

Type of CNA

``` r
RAD21$CNA
```

    ## [1] "AMP"

Mean of RAD21 cell line counts

``` r
mean(as.numeric(RAD21[which(!is.na(RAD21[5:18])) + 4]))
```

    ## [1] -0.2423539

Median of RAD21 cell line counts

``` r
median(as.numeric(RAD21[which(!is.na(RAD21[5:18])) + 4]))
```

    ## [1] -0.4356688

### ELOC

``` r
ELOC <- downreg[downreg$Gene == "ELOC" & downreg$Freq == "1.40%", ]
```

Type of CNA

``` r
ELOC$CNA
```

    ## [1] "AMP"

Mean of ELOC cell line counts

``` r
mean(as.numeric(ELOC[which(!is.na(ELOC[5:18])) + 4]))
```

    ## [1] -0.602545

Median of ELOC cell line counts

``` r
median(as.numeric(ELOC[which(!is.na(ELOC[5:18])) + 4]))
```

    ## [1] -0.9628831

### RB1

``` r
RB1 <- downreg[downreg$Gene == "RB1" & downreg$Freq == "1.80%", ]
```

Type of CNA

``` r
RB1$CNA
```

    ## [1] "HOMDEL"

Mean of RB1 cell line counts

``` r
mean(as.numeric(RB1[which(!is.na(RB1[5:18])) + 4]))
```

    ## [1] -0.5633385

Median of RB1 cell line counts

``` r
median(as.numeric(RB1[which(!is.na(RB1[5:18])) + 4]))
```

    ## [1] -0.9174057

### CDK4

``` r
CDK4 <- downreg[downreg$Gene == "CDK4" & downreg$Freq == "2.10%", ]
```

Type of CNA

``` r
CDK4$CNA
```

    ## [1] "AMP"

Mean of CDK4 cell line counts

``` r
mean(as.numeric(CDK4[which(!is.na(CDK4[5:18])) + 4]))
```

    ## [1] -0.8808549

Median of CDK4 cell line counts

``` r
median(as.numeric(CDK4[which(!is.na(CDK4[5:18])) + 4]))
```

    ## [1] -0.7170509

### Primary with metastasis vs without

Defining data sets from
<https://www.cbioportal.org/study/summary?id=msk_met_2021>. prim\_met
has all of the mutated genes data from samples prom primary tumors with
metastasis, while prim\_no\_met has the mutated genes data from primary
tumors with no metastasis

``` r
prim_met <- read.csv("prim_met.csv")
prim_no_met <- read.csv("prim_no_met.csv")
```

### BCL2L1

Frequency of BCL2L1 in primary with metastasis

``` r
prim_met[prim_met$Gene == "BCL2L1", "Freq"]
```

    ## [1] "0.20%"

Frequency of BCL2L1 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "BCL2L1", "Freq"]
```

    ## [1] "0.40%"

### ELOC

Frequency of ELOC in primary with metastasis

``` r
prim_met[prim_met$Gene == "ELOC", "Freq"]
```

    ## [1] "0.10%"

Frequency of ELOC in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "ELOC", "Freq"]
```

    ## [1] "<0.1%"

### RB1

Frequency of RB1 in primary with metastasis

``` r
prim_met[prim_met$Gene == "RB1", "Freq"]
```

    ## [1] "4.90%"

Frequency of RB1 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "RB1", "Freq"]
```

    ## [1] "4.50%"

### CDK4

Frequency of CDK4 in primary with metastasis

``` r
prim_met[prim_met$Gene == "CDK4", "Freq"]
```

    ## [1] "0.40%"

Frequency of CDK4 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "CDK4", "Freq"]
```

    ## [1] "0.50%"
