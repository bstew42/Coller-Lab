Dormancy\_Analysis
================

### Defining Data Sets

Chloe’s Dormancy Data

``` r
dormant <- read.csv("Dormant_Data.csv")
```

Mutated Genes from
<https://www.cbioportal.org/study/summary?id=msk_met_2021>

``` r
mutated <- read.csv("Mutated_Genes.csv")
```

Changing the column name for gene to “Gene” in Chloe’s dormancy data and
innerjoining the two sets by the Gene column

``` r
colnames(dormant)[4] <- "Gene"
data <- dplyr::inner_join(dormant, mutated, by = "Gene")
```

### Quick numbers regarding genes the two sets share

Number of genes in Chloe’s data

``` r
length(dormant$Gene)
```

    ## [1] 5878

Number of genes in cbioportal data

``` r
length(mutated$Gene)
```

    ## [1] 477

Number of Genes shared between the two

``` r
length(data$Gene)
```

    ## [1] 188

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

    ## [1] 7

The upregulated genes

``` r
upreg$Gene
```

    ## [1] "CDKN1A" "HLA-A"  "BCL2L1" "MAPK3"  "SESN2"  "PNRC1"  "SMAD3"

Frequencies of upregulated genes from cbioportal data

``` r
rbind(upreg$Gene, upreg$Freq)
```

    ##      [,1]     [,2]    [,3]     [,4]    [,5]    [,6]    [,7]   
    ## [1,] "CDKN1A" "HLA-A" "BCL2L1" "MAPK3" "SESN2" "PNRC1" "SMAD3"
    ## [2,] "0.90%"  "0.90%" "0.20%"  "0.50%" "0.70%" "0.40%" "1.30%"

CDKN1A and SMAD3 have the highest frequencies in the upregulated genes.
I will examine those two a little closer.

### CDKN1A

``` r
CDKN1A <- upreg[upreg$Gene == "CDKN1A", ]
```

Mean of CDKN1A cell line counts

``` r
mean(as.numeric(CDKN1A[which(!is.na(CDKN1A[5:18])) + 4]))
```

    ## [1] 1.513284

Median of CDKN1A cell line counts

``` r
median(as.numeric(CDKN1A[which(!is.na(CDKN1A[5:18])) + 4]))
```

    ## [1] 1.55555

### SMAD3

``` r
SMAD3 <- upreg[upreg$Gene == "SMAD3", ]
```

Mean of SMAD3 cell line counts

``` r
mean(as.numeric(SMAD3[which(!is.na(SMAD3[5:18])) + 4]))
```

    ## [1] 0.6396924

Median of SMAD3 cell line counts

``` r
median(as.numeric(SMAD3[which(!is.na(SMAD3[5:18])) + 4]))
```

    ## [1] 0.5689799

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

    ## [1] 24

The downregulated genes

``` r
downreg$Gene
```

    ##  [1] "EIF1AX" "NPM1"   "EZH2"   "DNMT1"  "CDK4"   "NUF2"   "CHEK1"  "RAD51D"
    ##  [9] "BLM"    "TOP1"   "BRCA1"  "POLD1"  "RHOA"   "XPO1"   "MSH2"   "BRCA2" 
    ## [17] "RB1"    "EIF4E"  "ELOC"   "SRSF2"  "PBRM1"  "RAD21"  "FANCA"  "XRCC2"

Frequencies of downregulated genes

``` r
rbind(downreg$Gene, downreg$Freq)
```

    ##      [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]     [,9]   
    ## [1,] "EIF1AX" "NPM1"  "EZH2"  "DNMT1" "CDK4"  "NUF2"  "CHEK1" "RAD51D" "BLM"  
    ## [2,] "0.60%"  "0.30%" "1.00%" "1.70%" "0.40%" "0.80%" "0.60%" "0.30%"  "1.60%"
    ##      [,10]   [,11]   [,12]   [,13]   [,14]   [,15]   [,16]   [,17]   [,18]  
    ## [1,] "TOP1"  "BRCA1" "POLD1" "RHOA"  "XPO1"  "MSH2"  "BRCA2" "RB1"   "EIF4E"
    ## [2,] "1.00%" "2.30%" "2.00%" "0.80%" "1.10%" "1.40%" "4.00%" "4.90%" "0.20%"
    ##      [,19]   [,20]   [,21]   [,22]   [,23]   [,24]  
    ## [1,] "ELOC"  "SRSF2" "PBRM1" "RAD21" "FANCA" "XRCC2"
    ## [2,] "0.10%" "0.50%" "3.40%" "0.90%" "1.70%" "0.40%"

BRCA2, RB1, and PRBM1 have the highest frequency of mutation in down
regulated genes. I will examine those three a little closer

### BRCA2

``` r
BRCA2 <- downreg[downreg$Gene == "BRCA2", ]
```

Mean of BRCA2 cell line counts

``` r
mean(as.numeric(BRCA2[which(!is.na(BRCA2[5:18])) + 4]))
```

    ## [1] -1.123837

Median of BRCA2 cell line counts

``` r
median(as.numeric(BRCA2[which(!is.na(BRCA2[5:18])) + 4]))
```

    ## [1] -1.069951

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

## Extremely up or down regulated

Now we will look into the genes that are either very upregulated or very
donwregulated. We will define these genes by having at least 6 cell line
counts with either all upregulated or all downregulated.

``` r
temp_data <-dplyr::inner_join(dormant, mutated, by = "Gene")

extr_data <- temp_data[temp_data$cell_line_counts >= 6, ]
extr_up <- extr_data[extr_data$downreg_counts == 0, ]
extr_down <- extr_data[extr_data$upreg_counts == 0, ]
```

### Extreme Up

Frequencies of the very upregulated genes

``` r
rbind(extr_up$Gene, extr_up$Freq)
```

    ##      [,1]    
    ## [1,] "CDKN1A"
    ## [2,] "0.90%"

### Extreme Down

Frequencies of the very downregulated genes

``` r
rbind(extr_down$Gene, extr_down$Freq)
```

    ##     
    ## [1,]
    ## [2,]

## Genes of Most Interest

### Up: CDKN1A

### Down: BRCA2

## Primary with metastasis vs without

Defining data sets from
<https://www.cbioportal.org/study/summary?id=msk_met_2021>. prim\_met
has all of the mutated genes data from samples prom primary tumors with
metastasis, while prim\_no\_met has the mutated genes data from primary
tumors with no metastasis

``` r
prim_met <- read.csv("prim_met.csv")
prim_no_met <- read.csv("prim_no_met.csv")
```

### CDKN1A

Frequency of CDKN1A in primary with metastasis

``` r
prim_met[prim_met$Gene == "CDKN1A", "Freq"]
```

    ## [1] "1.10%"

Frequency of CDKN1A in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "CDKN1A", "Freq"]
```

    ## [1] "1.20%"

### BRCA2

Frequency of BRCA2 in primary with metastasis

``` r
prim_met[prim_met$Gene == "BRCA2", "Freq"]
```

    ## [1] "3.80%"

Frequency of BRCA2 in primary without metastasis

``` r
prim_no_met[prim_no_met$Gene == "BRCA2", "Freq"]
```

    ## [1] "4.40%"

Interesting that BRCA2 is mutated more in samples without metastasis
