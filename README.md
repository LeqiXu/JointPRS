# JointPRS
JointPRS is a cross-population PRS model that only requires GWAS summary statistics and LD reference panel from multiple populations. It has two versions: 
- **JointPRS-auto**: **no need** for a validation dataset.
- **JointPRS**: requires a validation dataset for tuning parameters. 

## Version History

## Getting Started
In this section, we will offer step-by-step guidance on JointPRS implementation.

### 1. JointPRS Installation
For the first time, you need to use the following code to install JointPRS:
```
git clone https://github.com/LeqiXu/JointPRS.git
cd JointPRS
conda env create -f environment.yml
conda activate JointPRS
python setup.py build_ext --inplace
```
After this, you only need to use
```
conda activate JointPRS
```

### 2. LD Reference Panel Download
We use reference panels from [PRScsx](https://github.com/getian107/PRScsx#getting-started) and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory

- **1KG**: This subfolder should contain LD reference panels constructed using the 1000 Genomes Project phase 3 samples and the corresponding SNP information file.
- **UKBB**: This subfolder should contain LD reference panels constructed using the UK Biobank data and the corresponding SNP information file.

Then place the downloaded LD reference panels and the SNP information file into their corresponding subfolders.

### 3. Summary Statistics Preparation
We require the following format for summary statistics input (including the header line):
```
SNP         A1  A2  BETA         P
rs3934834	  C	  T	  0.0063086	   0.00512
rs3766192	  T	  C	  0.00761278	  5.14e-06
rs9442372	  G	  A	  0.00690567	  2.66e-05
...
```

Here
- `SNP`: the rs ID.
-  `A1`: the effect allele.
-  `A2`: the alternative allele.
-  `BETA`: the effect of the A1 allele, which is only used to determine the direction of an association.
-  `P`: the p-value of the effect, which is used to calculate the standardized effect size.

In addition, you need to obtain the sample size for the summary statistics, and take the median value if the sample size is different across SNPs.

### 4. JointPRS Implementation
In this section, we assume we need to model four populations jointly. Please modify the code to reflect the actual number of populations present in your data.

#### 4.1 Preparation
```
conda activate JointPRS

JointPRS_path=
reference_path= ; type=
bim_path= ; bim_prefix=
outcome_path=
param_phi=
chr=

pop1= ; pop2= ; pop3= ; pop4=
r1= ; r2= ; r3= ; r4=
sst1= ; sst2= ; sst3= ; sst4= 
sample_size1= ; sample_size2= ; sample_size3= ; sample_size4= 
```

- `${JointPRS_path}`: full path to the JointPRS software directory.
- `${reference_path}`: full path to the reference directory; `${type}`: **1KG** or **UKBB** if you construct two subfolders as recommended.
- `${bim_path}`: full path to the bim file for the target dataset; `${bim_prefix}`: prefix of the bim file for the target dataset.
- `${outcome_path}`: full path to the outcome directory.
- `${param_phi}`: **remove** this line if you use **JointPRS-auto**; use set **{1e-06, 1e-04, 1e-02, 1e+00}** if you use **JointPRS**.
- `${chr}`: the chromosome we want to consider (1-22) and we recommend estimate 22 chromosomes in parallel.

- `${pop1},${pop2},${pop3},${pop4}`: population name from set **{EUR,EAS,AFR,SAS,AMR}**.
- `${r1},${r2},${r3},${r4}`: upper bound for the correlation pairs, `${r1} * ${r2}` represents the upper bound for the correlation between `${pop1}` and `${pop2}`. And we recommand the following setting:
  * `r1=1,r2=1,r3=1,r4=1`: if you use **JointPRS-auto** or use **JointPRS** and `${param_phi}` comes from set **{1e-04, 1e-02, 1e+00}**. This setting represents the positive correlation assumption.
  * `r1=0,r2=0,r3=0,r4=0`: if you use **JointPRS** and `${param_phi}` comes from set **{1e-06}**. This setting represents no correlation assumption and is equivalent to [the PRScsx model](https://github.com/getian107/PRScsx/blob/master/README.md#prs-csx).
- `${sst1},${sst2},${sst3},${sst4}`: full path and full name of the summary statistics for the corresponding poopulation.
- `${sample_size1},${sample_size2},${sample_size3},${sample_size4}`: sample size for the corresponding summary statistics, and take the median value if the sample size is different across SNPs.

#### 4.2 JointPRS-auto
```
python ${JointPRS_path}/JointPRS.py \
--ref_dir=${refernce_path}/${type} \
--bim_prefix=${bim_path}/${bim_prefix} \
--pop=${pop1},${pop2},${pop3},${pop4} \
--rho_cons=${r1},${r2},${r3},${r4} \
--sst_file=${sst1},${sst2},${sst3},${sst4} \
--n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
--chrom=${chr} \
--out_dir=${outcome_path} \
--out_name=JointPRS_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
```

#### 4.3 JointPRS
```
python ${JointPRS_path}/JointPRS.py \
--ref_dir=${refernce_path}/${type} \
--bim_prefix=${bim_path}/${bim_prefix} \
--pop=${pop1},${pop2},${pop3},${pop4} \
--rho_cons=${r1},${r2},${r3},${r4} \
--sst_file=${sst1},${sst2},${sst3},${sst4} \
--n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
--chrom=${chr} \
--phi=${param_phi} \
--out_dir=${outcome_path} \
--out_name=JointPRS_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
```

## Example
The example contains EUR, EAS, AFR and GWAS summary statistics and a bim file for 500 SNPs on chromosome 1 for HDL.
- The summary statistics comes from [GLGC](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/).
- The bim file comes from the [1000 Genome Project](https://www.internationalgenome.org/data).

The following code is a demo to use the example data:
```
conda activate JointPRS

JointPRS_path=
reference_path= ; type=1KG
bim_path=${JointPRS_path}/example_data; bim_prefix=example
outcome_path=${JointPRS_path}/example_data
param_phi=1e-04

pop1=EUR; pop2=EAS; pop3=AFR; pop4=SAS
r1=1; r2=1; r3=1; r4=1
sst1=${JointPRS_path}/example_data/EUR_sumstat.txt; sst2=${JointPRS_path}/example_data/EAS_sumstat.txt
sst3=${JointPRS_path}/example_data/AFR_sumstat.txt; sst4=${JointPRS_path}/example_data/SAS_sumstat.txt 
sample_size1=885546; sample_size2=116404; sample_size3=90804; sample_size4=33953

python ${JointPRS_path}/JointPRS.py \
--ref_dir=${refernce_path}/${type} \
--bim_prefix=${bim_path}/${bim_prefix} \
--pop=${pop1},${pop2},${pop3},${pop4} \
--rho_cons=${r1},${r2},${r3},${r4} \
--sst_file=${sst1},${sst2},${sst3},${sst4} \
--n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
--chrom=${chr} \
--phi=${param_phi} \
--out_dir=${outcome_path} \
--out_name=JointPRS_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
```
Here you still need to specify `JointPRS_path` and `reference_path` by yourself.

## Acknowledgment

## Support

## Citation


