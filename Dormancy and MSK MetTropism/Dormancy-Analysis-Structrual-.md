Dormancy Analysis(Structural)
================

### Defining Data Sets

Chloe’s Dormancy Data

``` r
dormant <- read.csv("Dormant_Data.csv")
```

structural Genes from
<https://www.cbioportal.org/study/summary?id=msk_met_2021>

``` r
structural <- read.csv("Structural Variant.csv")
```

Changing the column name for gene to “Gene” in Chloe’s dormancy data and
innerjoining the two sets by the Gene column

``` r
colnames(dormant)[4] <- "Gene"
data <- dplyr::inner_join(dormant, structural, by = "Gene")
```

### Quick numbers regarding genes the two sets share

Number of genes in Chloe’s data

``` r
length(dormant$Gene)
```

    ## [1] 5878

Number of genes in cbioportal data

``` r
length(structural$Gene)
```

    ## [1] 2521

Number of Genes shared between the two

``` r
length(data$Gene)
```

    ## [1] 812

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

    ## [1] 32

The upregulated genes

``` r
upreg$Gene
```

    ##  [1] "ZMAT3"  "CDKN1A" "LRP1"   "CTSB"   "GLG1"   "PEPD"   "ITM2B"  "FURIN" 
    ##  [9] "PLCD3"  "BCL2L1" "PLPP1"  "CREM"   "TRIOBP" "MAPK3"  "ASCC3"  "UXS1"  
    ## [17] "RAP1A"  "PTK2B"  "SESN2"  "LAMC1"  "MYOF"   "TPM1"   "PNRC1"  "VOPP1" 
    ## [25] "ZFAND3" "SQSTM1" "IGSF8"  "CRYL1"  "SMAD3"  "ACSF2"  "SLC3A2" "HLA-E"

Frequencies of upregulated genes from cbioportal data

``` r
rbind(upreg$Gene, upreg$Freq)
```

    ##      [,1]    [,2]     [,3]    [,4]    [,5]    [,6]    [,7]    [,8]    [,9]   
    ## [1,] "ZMAT3" "CDKN1A" "LRP1"  "CTSB"  "GLG1"  "PEPD"  "ITM2B" "FURIN" "PLCD3"
    ## [2,] "<0.1%" "<0.1%"  "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%"
    ##      [,10]    [,11]   [,12]   [,13]    [,14]   [,15]   [,16]   [,17]   [,18]  
    ## [1,] "BCL2L1" "PLPP1" "CREM"  "TRIOBP" "MAPK3" "ASCC3" "UXS1"  "RAP1A" "PTK2B"
    ## [2,] "<0.1%"  "<0.1%" "<0.1%" "<0.1%"  "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%"
    ##      [,19]   [,20]   [,21]   [,22]   [,23]   [,24]   [,25]    [,26]    [,27]  
    ## [1,] "SESN2" "LAMC1" "MYOF"  "TPM1"  "PNRC1" "VOPP1" "ZFAND3" "SQSTM1" "IGSF8"
    ## [2,] "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%"  "<0.1%"  "<0.1%"
    ##      [,28]   [,29]   [,30]   [,31]    [,32]  
    ## [1,] "CRYL1" "SMAD3" "ACSF2" "SLC3A2" "HLA-E"
    ## [2,] "<0.1%" "<0.1%" "<0.1%" "<0.1%"  "<0.1%"

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

    ## [1] 68

The downregulated genes

``` r
downreg$Gene
```

    ##  [1] "DIAPH3"  "HNRNPA3" "HNRNPC"  "HNRNPM"  "RFC3"    "EIF1AX"  "RMI2"   
    ##  [8] "NPM1"    "SNRPE"   "MAPRE1"  "EZH2"    "YWHAE"   "CNTRL"   "CDK2"   
    ## [15] "TUBA1B"  "DNMT1"   "GINS2"   "CDK4"    "COPS3"   "NUF2"    "TPM3"   
    ## [22] "CHEK1"   "CCT2"    "CENPV"   "RAD51D"  "BLM"     "TOP1"    "PRC1"   
    ## [29] "BRCA1"   "TACC3"   "POLD1"   "RHOA"    "XPO1"    "XRN2"    "PHGDH"  
    ## [36] "CLSPN"   "MSH2"    "PABPN1"  "NUP107"  "PPP2CA"  "PSMD14"  "MRPS15" 
    ## [43] "EXOSC8"  "TUBGCP3" "LDLR"    "EIF5A"   "GNL2"    "CENPE"   "BRCA2"  
    ## [50] "RB1"     "EIF4A3"  "SAE1"    "ENSA"    "PRIM2"   "TWF1"    "EIF4E"  
    ## [57] "ELOC"    "ABHD3"   "SRSF2"   "NUP35"   "PBRM1"   "RAD21"   "PPIP5K1"
    ## [64] "FASN"    "THOP1"   "MRPS22"  "FANCA"   "XRCC2"

Frequencies of downregulated genes

``` r
rbind(downreg$Gene, downreg$Freq)
```

    ##      [,1]     [,2]      [,3]     [,4]     [,5]    [,6]     [,7]    [,8]   
    ## [1,] "DIAPH3" "HNRNPA3" "HNRNPC" "HNRNPM" "RFC3"  "EIF1AX" "RMI2"  "NPM1" 
    ## [2,] "<0.1%"  "<0.1%"   "<0.1%"  "<0.1%"  "<0.1%" "<0.1%"  "<0.1%" "<0.1%"
    ##      [,9]    [,10]    [,11]   [,12]   [,13]   [,14]   [,15]    [,16]   [,17]  
    ## [1,] "SNRPE" "MAPRE1" "EZH2"  "YWHAE" "CNTRL" "CDK2"  "TUBA1B" "DNMT1" "GINS2"
    ## [2,] "<0.1%" "<0.1%"  "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%"  "0.10%" "<0.1%"
    ##      [,18]   [,19]   [,20]   [,21]   [,22]   [,23]   [,24]   [,25]    [,26]  
    ## [1,] "CDK4"  "COPS3" "NUF2"  "TPM3"  "CHEK1" "CCT2"  "CENPV" "RAD51D" "BLM"  
    ## [2,] "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%"  "<0.1%"
    ##      [,27]   [,28]   [,29]   [,30]   [,31]   [,32]   [,33]   [,34]   [,35]  
    ## [1,] "TOP1"  "PRC1"  "BRCA1" "TACC3" "POLD1" "RHOA"  "XPO1"  "XRN2"  "PHGDH"
    ## [2,] "<0.1%" "<0.1%" "0.10%" "0.20%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%"
    ##      [,36]   [,37]   [,38]    [,39]    [,40]    [,41]    [,42]    [,43]   
    ## [1,] "CLSPN" "MSH2"  "PABPN1" "NUP107" "PPP2CA" "PSMD14" "MRPS15" "EXOSC8"
    ## [2,] "<0.1%" "<0.1%" "<0.1%"  "<0.1%"  "<0.1%"  "<0.1%"  "<0.1%"  "<0.1%" 
    ##      [,44]     [,45]   [,46]   [,47]   [,48]   [,49]   [,50]   [,51]    [,52]  
    ## [1,] "TUBGCP3" "LDLR"  "EIF5A" "GNL2"  "CENPE" "BRCA2" "RB1"   "EIF4A3" "SAE1" 
    ## [2,] "<0.1%"   "<0.1%" "<0.1%" "<0.1%" "<0.1%" "0.10%" "0.30%" "<0.1%"  "<0.1%"
    ##      [,53]   [,54]   [,55]   [,56]   [,57]   [,58]   [,59]   [,60]   [,61]  
    ## [1,] "ENSA"  "PRIM2" "TWF1"  "EIF4E" "ELOC"  "ABHD3" "SRSF2" "NUP35" "PBRM1"
    ## [2,] "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "<0.1%" "0.20%"
    ##      [,62]   [,63]     [,64]   [,65]   [,66]    [,67]   [,68]  
    ## [1,] "RAD21" "PPIP5K1" "FASN"  "THOP1" "MRPS22" "FANCA" "XRCC2"
    ## [2,] "0.10%" "<0.1%"   "<0.1%" "<0.1%" "<0.1%"  "<0.1%" "<0.1%"

RB1, PBRM1, RAD21, TACC3, BRCA2, and DNMT1 all have frequencies of 0.1%
or greater

### RB1

``` r
RB1 <- downreg[downreg$Gene == "RB1", ]
```

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

### PBRM1

``` r
PBRM1 <- downreg[downreg$Gene == "PBRM1", ]
```

Mean of PBRM1 cell line counts

``` r
mean(as.numeric(PBRM1[which(!is.na(PBRM1[5:18])) + 4]))
```

    ## [1] -0.4485714

Median of PBRM1 cell line counts

``` r
median(as.numeric(PBRM1[which(!is.na(PBRM1[5:18])) + 4]))
```

    ## [1] -0.893644

### RAD21

``` r
RAD21 <- downreg[downreg$Gene == "RAD21", ]
```

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

### BRAC2

``` r
BRAC2 <- downreg[downreg$Gene == "BRCA2", ]
```

Mean of BRAC2 cell line counts

``` r
mean(as.numeric(BRAC2[which(!is.na(BRAC2[5:18])) + 4]))
```

    ## [1] -1.123837

Median of BRAC2 cell line counts

``` r
median(as.numeric(BRAC2[which(!is.na(BRAC2[5:18])) + 4]))
```

    ## [1] -1.069951

### TACC3

``` r
TACC3 <- downreg[downreg$Gene == "TACC3", ]
```

Mean of TACC3 cell line counts

``` r
mean(as.numeric(TACC3[which(!is.na(TACC3[5:18])) + 4]))
```

    ## [1] -0.5826649

Median of TACC3 cell line counts

``` r
median(as.numeric(TACC3[which(!is.na(TACC3[5:18])) + 4]))
```

    ## [1] -1.015571

### DNMT1

``` r
DNMT1 <- downreg[downreg$Gene == "DNMT1", ]
```

Mean of DNMT1 cell line counts

``` r
mean(as.numeric(DNMT1[which(!is.na(DNMT1[5:18])) + 4]))
```

    ## [1] -0.8225876

Median of DNMT1 cell line counts

``` r
median(as.numeric(DNMT1[which(!is.na(DNMT1[5:18])) + 4]))
```

    ## [1] -0.960106

## Primary with metastasis vs without

Defining data sets from
<https://www.cbioportal.org/study/summary?id=msk_met_2021>. prim\_met
has all of the mutated genes data from samples from primary tumors with
metastasis, while prim\_no\_met has the mutated genes data from primary
tumors with no metastasis

``` r
prim_met <- read.csv("prim_met.csv")
prim_no_met <- read.csv("prim_no_met.csv")
```

### PBRM1

Frequency of PBRM1 in primary with metastasis

``` r
prim_met[prim_met$Gene == "PBRM1", "Freq"]
```

    ## [1] "3.60%"

Frequency of PBRM1 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "PBRM1", "Freq"]
```

    ## [1] "3.60%"

### RAD21

Frequency of RAD21 in primary with metastasis

``` r
prim_met[prim_met$Gene == "RAD21", "Freq"]
```

    ## [1] "0.90%"

Frequency of RAD21 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "RAD21", "Freq"]
```

    ## [1] "1.00%"

### TACC3

Frequency of TACC3 in primary with metastasis

``` r
prim_met[prim_met$Gene == "TACC3", "Freq"]
```

    ## character(0)

Frequency of TACC3 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "TACC3", "Freq"]
```

    ## character(0)

### DNMT1

Frequency of DNMT1 in primary with metastasis

``` r
prim_met[prim_met$Gene == "DNMT1", "Freq"]
```

    ## [1] "1.70%"

Frequency of DNMT1 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "DNMT1", "Freq"]
```

    ## [1] "2.00%"
