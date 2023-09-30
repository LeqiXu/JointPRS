# JointPRS
JointPRS is a cross-population PRS model that only requires GWAS summary statistics and LD reference panel from multiple populations. It has two versions: 
- **JointPRS-auto**: **no need** for a validation dataset.
- **JointPRS**: requires a validation dataset for tuning parameters. 

# Getting Started
In this section, we will offer step-by-step guidance on JointPRS implementation.

## 1. JointPRS Installation

## 2. LD Reference Panel Download
We use the reference panel from [PRScsx](https://github.com/getian107/PRScsx#getting-started) and you can download it following their guidance. We highly recommend you to cosntruct two folders under the reference folder.

- **1KG**: LD reference panels constructed using the 1000 Genomes Project phase 3 samples.
- **UKBB**: LD reference panels constructed using the UK Biobank data.

And put the corresponding LD reference panels as well as the SNP information file in the folder.

## 3. Summary Statistics Preparation

## 4. JointPRS Implementation
