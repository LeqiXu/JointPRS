# JointPRS
JointPRS is a cross-population PRS model that only requires GWAS summary statistics and LD reference panel from multiple populations. It has two versions: 
- **JointPRS-auto**: **no need** for a validation dataset.
- **JointPRS**: requires a validation dataset for tuning parameters. 

# Getting Started
In this section, we will offer step-by-step guidance on JointPRS implementation.

## 1. JointPRS Installation

## 2. LD Reference Panel Download
We use reference panels from [PRScsx](https://github.com/getian107/PRScsx#getting-started) and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory

- **1KG**: This subfolder should contain LD reference panels constructed using the 1000 Genomes Project phase 3 samples and the corresponding SNP information file.
- **UKBB**: This subfolder should contain LD reference panels constructed using the UK Biobank data and the corresponding SNP information file.

Then place the downloaded LD reference panels and the SNP information file into their corresponding subfolders.

## 3. Summary Statistics Preparation
We require the following format for summary statistics input (including the header line):

```
SNP         A1  A2  BETA           P
rs3131972   G   A   -0.000422623   0.869
rs3131969   G   A   0.000148544    0.955
rs1048488   T   C   0.000678053    0.798
...
```
Here
- `SNP`: the rs ID.
-  `A1`: the effect allele.
-  `A2`: the alternative allele.
-  `BETA`: the effect of the A1 allele, which is only used to determine the direction of an association.
-  `P`: the p-value of the effect, which is used to calculate the standardized effect size.

## 4. JointPRS Implementation
In this section, we assume we need to model four populations jointly. Please modify the code to reflect the actual number of populations present in your data.

### 4.1 Preparation
```
conda activate JointPRS

JointPRS_path=
reference_path= ;type=
bim_path= ;bim_prefix=
outcome_path=
param_phi=

pop1= ;pop2=; pop3= ;pop4=
r1= ;r2= ;r3= ;r4=
sst1= ;sst2= ;sst3= ;sst4= 
sample_size1= ;sample_size2= ;sample_size3= ;sample_size4= 
```

### 4.2 JointPRS-auto


### 4.3 JointPRS

```

```
