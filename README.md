# JointPRS
JointPRS is a cross-population PRS model that only requires GWAS summary statistics and LD reference panel from multiple populations. It has two versions: 
- **JointPRS-auto**: **no need** for a validation dataset.
- **JointPRS**: requires a validation dataset for tuning parameters. 

# Getting Started
In this section, we will offer step-by-step guidance on JointPRS implementation.

## 1. JointPRS Installation

## 2. LD Reference Panel Download
We use reference panels from [PRScsx](https://github.com/getian107/PRScsx#getting-started) and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory

- **1KG**: This subfolder should contain LD reference panels constructed using the 1000 Genomes Project phase 3 samples.
- **UKBB**: This subfolder should contain LD reference panels constructed using the UK Biobank data.

Place the LD reference panels and the SNP information file into their corresponding subfolders.

## 3. Summary Statistics Preparation
We require the following summary statistics input

```
SNP A1 A2 BETA P
rs3131972 G A -0.000422623 0.869
rs3131969 G A 0.000148544 0.955
rs1048488 T C 0.000678053 0.798
...
```

## 4. JointPRS Implementation
