# JointPRS
JointPRS is a multi-population PRS model that only requires GWAS summary statistics and LD reference panel from multiple populations. When individual-level tuning data is available, it adpots a data-adaptive apporach that combines meta-analysis and tuning strategies. Based on the availability of tuning data, JointPRS has two implementations: 
- **JointPRS-auto**: **No** tuning data scenario.
  - Compute the auto version of JointPRS.
- **JointPRS**: **Exist** tuning data scenario.
  - Compute the meta version and the tune version of JointPRS.
  - Utilize a data-adaptive approach to select the optimal version of JointPRS.

![JointPRS_workflow](https://github.com/user-attachments/assets/15ca8ae6-c786-473e-8526-88e5e699f964)

## Update History
**Sept 30, 2023**: Repository made public.

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
We use reference panels from [PRS-CSx](https://github.com/getian107/PRScsx#getting-started) and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory:

- **1KG**: This subfolder should contain LD reference panels constructed using the 1000 Genomes Project phase 3 samples and the corresponding SNP information file.
- **UKBB**: This subfolder should contain LD reference panels constructed using the UK Biobank data and the corresponding SNP information file.

Then place the downloaded LD reference panels and the SNP information file into their corresponding subfolders.

### 3. Summary Statistics Preparation
We require the following format for summary statistics input (including the header line):
```
SNP               A1      A2      BETA            P
rs3934834	  C	  T	  0.0063086	  0.00512
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

### 4. JointPRS Model Implementation
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
- `${param_phi}`: the global shrinkage prior that will be discussed later based on different versions of JointPRS.
- `${chr}`: the chromosome we want to consider (1-22) and we recommend estimate 22 chromosomes in parallel.

- `${pop1},${pop2},${pop3},${pop4}`: population name from set **{EUR,EAS,AFR,SAS,AMR}**.
- `${r1},${r2},${r3},${r4}`: upper bound for the correlation pairs, `${r1} * ${r2}` represents the upper bound for the correlation between `${pop1}` and `${pop2}`. And we recommand the following setting:
  * `r1=1,r2=1,r3=1,r4=1`: if `${param_phi}` does **not** come from set **{1e-06}**. This setting represents the positive correlation assumption.
  * `r1=0,r2=0,r3=0,r4=0`: if `${param_phi}` comes from set **{1e-06}**. This setting represents no correlation assumption and is equivalent to [the PRS-CSx model](https://github.com/getian107/PRScsx/blob/master/README.md#prs-csx).
- `${sst1},${sst2},${sst3},${sst4}`: full path and full name of the summary statistics for the corresponding poopulation.
- `${sample_size1},${sample_size2},${sample_size3},${sample_size4}`: sample size for the corresponding summary statistics, and take the median value if the sample size is different across SNPs.

#### 4.2 No Tuning Data Scenario: JointPRS-auto
When there is no tuning data, we compute the auto version of JointPRS.

- **Auto Version**:
  - Auto version utilizes the original GWAS summary statistics as input data and does not require any parameter tuning.
  - Code:
    ```
    python ${JointPRS_path}/JointPRS.py \
    --ref_dir=${reference_path}/${type} \
    --bim_prefix=${bim_path}/${bim_prefix} \
    --pop=${pop1},${pop2},${pop3},${pop4} \
    --rho_cons=${r1},${r2},${r3},${r4} \
    --sst_file=${sst1},${sst2},${sst3},${sst4} \
    --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
    --chrom=${chr} \
    --out_dir=${outcome_path} \
    --out_name=JointPRS_auto_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
    ```

#### 4.3 Exist Tuning Data Scenario: JointPRS
When there exist tuning data, we need to compute the meta version and the tune version of JointPRS.

- **Meta Version**:
  - Meta version leverages meta-analysis to integrate the orginal GWAS summary statistics with tuning datasets, utilizing [the METAL software](https://csg.sph.umich.edu/abecasis/metal/).
  - Example Code:
    ```
    /gpfs/gibbs/pi/zhao/lx94/Software/generic-metal/metal
    SCHEME   STDERR
    MARKER   SNP
    WEIGHT   N
    ALLELE   A1 A2
    FREQ     MAF
    EFFECT   BETA
    STDERR   SE
    PVAL     P

    PROCESS /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/HDL_AFR_inter_clean.txt
    PROCESS /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/UKB_glm/same_cohort/HDL_AFR_UKB_val_1.txt
    OUTFILE /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/with_UKB_meta/same_cohort/HDL_AFR_inter_UKB_val_1_meta .tbl
    ANALYZE
    QUIT
    ```
  - Meta version utilizes the updated GWAS summary statistics as input data and does not required any parameter tuning.
  - Code:
    ```
    python ${JointPRS_path}/JointPRS.py \
    --ref_dir=${reference_path}/${type} \
    --bim_prefix=${bim_path}/${bim_prefix} \
    --pop=${pop1},${pop2},${pop3},${pop4} \
    --rho_cons=${r1},${r2},${r3},${r4} \
    --sst_file=${sst1},${sst2},${sst3},${sst4} \
    --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
    --chrom=${chr} \
    --out_dir=${outcome_path} \
    --out_name=JointPRS_meta_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
    ```
    
- **Tune Version**:
  - Tune version utilizes the original GWAS summary statistics as input data and choose the global shrinkage parameter from a broader range **{1e-06,1e-04,1e-02,1e+00,auto}**.
  - Code for **{auto}** is equivalent to the implementation of JointPRS-auto:
    ```
    python ${JointPRS_path}/JointPRS.py \
    --ref_dir=${reference_path}/${type} \
    --bim_prefix=${bim_path}/${bim_prefix} \
    --pop=${pop1},${pop2},${pop3},${pop4} \
    --rho_cons=${r1},${r2},${r3},${r4} \
    --sst_file=${sst1},${sst2},${sst3},${sst4} \
    --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
    --chrom=${chr} \
    --out_dir=${outcome_path} \
    --out_name=JointPRS_auto_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
    ```
  - Code for **{1e-06,1e-04,1e-02,1e+00}**:
    ```
    python ${JointPRS_path}/JointPRS.py \
    --ref_dir=${reference_path}/${type} \
    --bim_prefix=${bim_path}/${bim_prefix} \
    --pop=${pop1},${pop2},${pop3},${pop4} \
    --rho_cons=${r1},${r2},${r3},${r4} \
    --sst_file=${sst1},${sst2},${sst3},${sst4} \
    --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
    --chrom=${chr} \
    --phi=${param_phi} \
    --out_dir=${outcome_path} \
    --out_name=JointPRS_tune_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
    ```

#### 4.4 Output
**JointPRS has two types of outputs:**
- Posterior SNP effect size estimates for each chromosome for each population (pst_eff):
```
1 rs3934834 1005806 T C -3.715860e-05 6.496727e-03
1 rs3766192 1017197 T C -1.608738e-05 1.829326e-03
1 rs9442372 1018704 G A 7.444141e-05  1.959886e-02
...
```
 Seven columns represent chromosome, rsID, base position, A1, A2, posterior effect size estimates and local shrinkage estimates
- Correlation matrix estimates for each chromosome (pst_corr):
```
9.999999999999896749e-01 4.676351646602236456e-01 7.278592525977222172e-01 5.202891308281465399e-01
4.676351646602236456e-01 9.999999999999896749e-01 6.543938271240943294e-01 3.986037958722946639e-01
7.278592525977222172e-01 6.543938271240943294e-01 9.999999999999896749e-01 4.418160360865302505e-01
5.202891308281465399e-01 3.986037958722946639e-01 4.418160360865302505e-01 9.999999999999896749e-01
```

**For the polygenic risk score calculation**, we recommend utilizing the PLINK/2 --score command, which requires three columns: rsID, A1, and posterior effect size estimates. Detailed documentation for the --score command can be found on the [PLINK/2 website](https://www.cog-genomics.org/plink/2.0/score) 

**We recommend integrating the beta files from all 22 chromosomes into a single beta file to be used as input in PLINK/2.**
- Example Code:
  ```
  target_geno=
  target_pop=
  trait=

  plink2 --bfile ${target_geno} \
  --double-id \
  --threads 1 \
  --score JointPRS_tune_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}_${target_pop}_pst_eff_a1_b0.5_phi${param_phi}.txt 2 4 6 cols=+scoresums,-scoreavgs \
  --out ${trait}_${pop1}_${pop2}_${pop3}_${pop4}_${type}_JointPRS_tune_${target_pop}_phi${param_phi}
  ```
**Note:** The parameter `cols=+scoresums,-scoreavgs` is unnecessary if a single integrated beta file for all chromosomes is used as input to compute the score. However, it is **required** if the scores are computed separately for each of the 22 chromosomes.

**For further score combination**, it depends on the availability of the tuning dataset:
- If there is **no** tuning dataset, we use **JointPRS-auto**, and the genome-wide score of the target population from the previous step will be the final score for complex traits prediction.
- If there **exists** a tuning dataset, we further use a data-adaptive approach to select between the meta version and the tune version, as detailed in section **5. JointPRS Data-Adaptive Approach**.

#### 5. JointPRS Data-Adaptive Approach
JointPRS data-adaptive approach is performed only when the tuning dataset is available, aiming to select the optimal PRS between the meta version and tune version.

- **Meta Version**:
  - Obtain JointPRS_meta_PRS: Obtain the PRS for the meta version of JointPRS for the target population, and name it JointPRS_meta_PRS.
- **Tune Version**:
  - Obtain JointPRS_tune_PRS_${param_phi}_${pop}: Obtain the PRS for the tune version of JointPRS for each global shrinkage parameter `${param_phi}` from set **{1e-06, 1e-04, 1e-02, 1e+00, auto}**, and calculate the polygenic score of the target population using each of the discovery population posterior SNP effect size estimates (four discovery populations in our example and four PRS scores in this step).
  - Obtain JointPRS_tune_PRS_${param_phi}_linear: Perform a linear combination of the polygenic scores across discovery populations for each global shrinkage parameter `${param_phi}` from set **{1e-06, 1e-04, 1e-02, 1e+00, auto}**.
  - Obtain JointPRS_tune_PRS_optimal_linear: Select the optimal PRS across global shrinkage parameter choices in the tuning dataset.
  - Example Code:
    ```
    library(data.table)
    library(stringr)
    library(dplyr)

    param_list = c("auto","1e-06","1e-04","1e-02","1e+00")
    type="1KG"

    for(s in c(1:5)){
    for(trait in c("HDL")){
    for(pop in c("EAS","AFR","SAS")){
    
    Trait_JointPRS_phiauto_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_",type,"_JointPRS_tune_",pop,"_phiauto.sscore"))
    Trait_JointPRS_phi1e_06_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_",type,"_JointPRS_tune_",pop,"_phi1e-06.sscore"))
    Trait_JointPRS_phi1e_04_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_",type,"_JointPRS_tune_",pop,"_phi1e-04.sscore"))
    Trait_JointPRS_phi1e_02_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_",type,"_JointPRS_tune_",pop,"_phi1e-02.sscore"))
    Trait_JointPRS_phi1e_00_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_",type,"_JointPRS_tune_",pop,"_phi1e+00.sscore"))

    scale_pheno = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/split/",trait,"_scale_",pop,"_val_",s,"_doubleid.tsv"))
    scale_pheno =scale_pheno[,c(1,3)]
    colnames(scale_pheno) = c("eid","pheno")

    ## validation
    ## JointPRS
    Trait_JointPRS_phiauto_val = Trait_JointPRS_phiauto_val[,c(1,5,6,7,8)]
    colnames(Trait_JointPRS_phiauto_val) = c("eid","EUR","EAS","AFR","SAS")
    Trait_JointPRS_phiauto_val = scale_pheno[Trait_JointPRS_phiauto_val, on = .(eid = eid)]
    Trait_JointPRS_phiauto_val = Trait_JointPRS_phiauto_val[,c(3:6) := lapply(.SD,scale),.SDcols = c(3:6)]
    Trait_JointPRS_phiauto_val = na.omit(Trait_JointPRS_phiauto_val)
    
    Trait_JointPRS_phi1e_06_val = Trait_JointPRS_phi1e_06_val[,c(1,5,6,7,8)]
    colnames(Trait_JointPRS_phi1e_06_val) = c("eid","EUR","EAS","AFR","SAS")
    Trait_JointPRS_phi1e_06_val = scale_pheno[Trait_JointPRS_phi1e_06_val, on = .(eid = eid)]
    Trait_JointPRS_phi1e_06_val = Trait_JointPRS_phi1e_06_val[,c(3:6) := lapply(.SD,scale),.SDcols = c(3:6)]
    Trait_JointPRS_phi1e_06_val = na.omit(Trait_JointPRS_phi1e_06_val)

    Trait_JointPRS_phi1e_04_val = Trait_JointPRS_phi1e_04_val[,c(1,5,6,7,8)]
    colnames(Trait_JointPRS_phi1e_04_val) = c("eid","EUR","EAS","AFR","SAS")
    Trait_JointPRS_phi1e_04_val = scale_pheno[Trait_JointPRS_phi1e_04_val, on = .(eid = eid)]
    Trait_JointPRS_phi1e_04_val = Trait_JointPRS_phi1e_04_val[,c(3:6) := lapply(.SD,scale),.SDcols = c(3:6)]
    Trait_JointPRS_phi1e_04_val = na.omit(Trait_JointPRS_phi1e_04_val)
    
    Trait_JointPRS_phi1e_02_val = Trait_JointPRS_phi1e_02_val[,c(1,5,6,7,8)]
    colnames(Trait_JointPRS_phi1e_02_val) = c("eid","EUR","EAS","AFR","SAS")
    Trait_JointPRS_phi1e_02_val = scale_pheno[Trait_JointPRS_phi1e_02_val, on = .(eid = eid)]
    Trait_JointPRS_phi1e_02_val = Trait_JointPRS_phi1e_02_val[,c(3:6) := lapply(.SD,scale),.SDcols = c(3:6)]
    Trait_JointPRS_phi1e_02_val = na.omit(Trait_JointPRS_phi1e_02_val)

    Trait_JointPRS_phi1e_00_val = Trait_JointPRS_phi1e_00_val[,c(1,5,6,7,8)]
    colnames(Trait_JointPRS_phi1e_00_val) = c("eid","EUR","EAS","AFR","SAS")
    Trait_JointPRS_phi1e_00_val = scale_pheno[Trait_JointPRS_phi1e_00_val, on = .(eid = eid)]
    Trait_JointPRS_phi1e_00_val = Trait_JointPRS_phi1e_00_val[,c(3:6) := lapply(.SD,scale),.SDcols = c(3:6)]
    Trait_JointPRS_phi1e_00_val = na.omit(Trait_JointPRS_phi1e_00_val)

    # JointPRS validation data select the best performed parameter
    lm_JointPRS_phiauto_val = lm(pheno ~ . + 0, data = Trait_JointPRS_phiauto_val[,c(2:6)])
    lm_JointPRS_phi1e_06_val = lm(pheno ~ . + 0, data = Trait_JointPRS_phi1e_06_val[,c(2:6)])
    lm_JointPRS_phi1e_04_val = lm(pheno ~ . + 0, data = Trait_JointPRS_phi1e_04_val[,c(2:6)])
    lm_JointPRS_phi1e_02_val = lm(pheno ~ . + 0, data = Trait_JointPRS_phi1e_02_val[,c(2:6)])
    lm_JointPRS_phi1e_00_val = lm(pheno ~ . + 0, data = Trait_JointPRS_phi1e_00_val[,c(2:6)])

    Trait_JointPRS_val_R2 = data.table(JointPRS_phiauto = summary(lm_JointPRS_phiauto_val)$`r.squared`,
                                    JointPRS_phi1e_06 = summary(lm_JointPRS_phi1e_06_val)$`r.squared`, 
                                    JointPRS_phi1e_04 = summary(lm_JointPRS_phi1e_04_val)$`r.squared`,
                                    JointPRS_phi1e_02 = summary(lm_JointPRS_phi1e_02_val)$`r.squared`, 
                                    JointPRS_phi1e_00 = summary(lm_JointPRS_phi1e_00_val)$`r.squared`)
    JointPRS_val_weight = data.table(rbind(lm_JointPRS_phiauto_val$coefficient,lm_JointPRS_phi1e_06_val$coefficient,lm_JointPRS_phi1e_04_val$coefficient,lm_JointPRS_phi1e_02_val$coefficient,lm_JointPRS_phi1e_00_val$coefficient))
                         
    ## best index
    JointPRS_index = which.max(Trait_JointPRS_val_R2)
    best_param = param_list[JointPRS_index]
    print(JointPRS_index)

    Trait_JointPRS_optimal_weight = JointPRS_val_weight[JointPRS_index,]
    Trait_JointPRS_optimal_phi = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_JointPRS_phi",best_param,"_beta.txt"))
      write.table(Trait_JointPRS_optimal_weight,paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_tune/",trait,"_JointPRS_linear_val_",s,"_EUR_EAS_AFR_SAS_weight_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    write.table(Trait_JointPRS_optimal_phi,paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_tune/",trait,"_JointPRS_linear_val_",s,"_EUR_EAS_AFR_SAS_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    }
    }
    }
    ``` 
- **Select Between Meta Version and Tune Version**:
  - Perform model selection using F-test for continuous traits and $\chi^2$-test for binary traits.
  - Example Code for continuous traits:
    ```
    library(data.table)
    library(stringr)
    library(dplyr)

    param_list = c("auto","1e-06","1e-04","1e-02","1e+00")
    type="1KG"

    for(s in c(1:5)){
    for(trait in c("HDL")){
    for(pop in c("EAS","AFR","SAS")){
    
    Trait_JointPRS_phiauto_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_AFR_SAS_",type,"_JointPRS_tune_",pop,"_phiauto.sscore"))
    
    scale_pheno = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/split/",trait,"_scale_",pop,"_val_",s,"_doubleid.tsv"))
    scale_pheno =scale_pheno[,c(1,3)]
    colnames(scale_pheno) = c("eid","pheno")

    ## validation
    ## JointPRS
    Trait_JointPRS_phiauto_val = Trait_JointPRS_phiauto_val[,c(1,5,6,7,8)]
    colnames(Trait_JointPRS_phiauto_val) = c("eid","EUR","EAS","AFR","SAS")
    Trait_JointPRS_phiauto_val = scale_pheno[Trait_JointPRS_phiauto_val, on = .(eid = eid)]
    Trait_JointPRS_phiauto_val = Trait_JointPRS_phiauto_val[,c(3:6) := lapply(.SD,scale),.SDcols = c(3:6)]
    Trait_JointPRS_phiauto_val = na.omit(Trait_JointPRS_phiauto_val)
    col_idx = c(2,which(colnames(Trait_JointPRS_phiauto_val) == pop))
    
    # JointPRS validation data select the best performed parameter
    lm_JointPRS_phiauto_val_sub = lm(pheno ~ . + 0, data = Trait_JointPRS_phiauto_val[,..col_idx])
    lm_JointPRS_phiauto_val_full = lm(pheno ~ . + 0, data = Trait_JointPRS_phiauto_val[,c(2:6)])

    R2_sub = data.table(JointPRS_phiauto = summary(lm_JointPRS_phiauto_val_sub)$`r.squared`)
    R2_full = data.table(JointPRS_phiauto = summary(lm_JointPRS_phiauto_val_full)$`r.squared`)

    p_value = anova(lm_JointPRS_phiauto_val_sub,lm_JointPRS_phiauto_val_full)$`Pr(>F)`
    p_value = na.omit(p_value)

    write.table(R2_sub,paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_meta/",trait,"_JointPRS_meta_val_",s,"_EUR_EAS_AFR_SAS_r2_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    write.table(R2_full,paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_tune/",trait,"_JointPRS_linear_val_",s,"_EUR_EAS_AFR_SAS_r2_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    write.table(p_value, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_meta/",trait,"_JointPRS_meta_val_",s,"_EUR_EAS_AFR_SAS_pvalue_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    }
    }
    }
    ```
  - Example Code for binary traits:
    ```
    library(data.table)
    library(stringr)
    library(dplyr)
    library(pROC)

    param_list = c("auto","1e-06","1e-04","1e-02","1e+00")
    type="1KG"

    for(s in c(1:5)){
    for(trait in c("CAD")){
    for(pop in c("EAS")){
    
    Trait_JointPRS_phiauto_val = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/same_cohort/JointPRS/",trait,"_EUR_EAS_",type,"_JointPRS_tune_",pop,"_phiauto.sscore"))
    
    pheno = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/split/",trait,"_",pop,"_val_",s,"_doubleid.tsv"))
    pheno =pheno[,c(1,3)]
    colnames(pheno) = c("eid","pheno")

    ## validation
    ## JointPRS
    Trait_JointPRS_phiauto_val = Trait_JointPRS_phiauto_val[,c(1,5,6)]
    colnames(Trait_JointPRS_phiauto_val) = c("eid","EUR","EAS")
    Trait_JointPRS_phiauto_val = pheno[Trait_JointPRS_phiauto_val, on = .(eid = eid)]
    Trait_JointPRS_phiauto_val = Trait_JointPRS_phiauto_val[,c(3:4) := lapply(.SD,scale),.SDcols = c(3:4)]
    Trait_JointPRS_phiauto_val = na.omit(Trait_JointPRS_phiauto_val)
    col_idx = c(2,which(colnames(Trait_JointPRS_phiauto_val) == pop))
    
    # JointPRS validation data select the best performed parameter
    glm_JointPRS_phiauto_val_sub = glm(pheno ~ . + 0, data = Trait_JointPRS_phiauto_val[,..col_idx],family=binomial(link="logit"))
    glm_JointPRS_phiauto_val_full = glm(pheno ~ . + 0, data = Trait_JointPRS_phiauto_val[,c(2:4)],family=binomial(link="logit"))

    AUC_sub = data.table(JointPRS_phiauto = roc(Trait_JointPRS_phiauto_val$pheno, predict(glm_JointPRS_phiauto_val_sub, type="response"), quiet=T, plot=F)$auc)
    AUC_full = data.table(JointPRS_phiauto = roc(Trait_JointPRS_phiauto_val$pheno, predict(glm_JointPRS_phiauto_val_full, type="response"), quiet=T, plot=F)$auc)

    p_value = anova(glm_JointPRS_phiauto_val_sub,glm_JointPRS_phiauto_val_full,test="Chisq")$`Pr(>Chi)`
    p_value = na.omit(p_value)

    write.table(AUC_sub,paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_meta/",trait,"_JointPRS_meta_val_",s,"_EUR_EAS_auc_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    write.table(AUC_full,paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_tune/",trait,"_JointPRS_linear_val_",s,"_EUR_EAS_auc_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    write.table(p_value, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/same_cohort/JointPRS_meta/",trait,"_JointPRS_meta_val_",s,"_EUR_EAS_pvalue_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
    }
    }
    }
    ```
  - Example Code for final model selection:
    - Continuous traits:
      ```
      if (R2_sub > 0.01 && (R2_full - R2_sub > 0) && pvalue < 0.05){
      JointPRS_PRS = JointPRS_tune_PRS_optimal_linear
      } else {
      JointPRS_PRS = JointPRS_meta_PRS
      }
      ```
    - Binary traits:
      ```
      if (AUC_sub > 0.501 && (AUC_full - AUC_sub > 0) && pvalue < 0.05){
      JointPRS_PRS = JointPRS_tune_PRS_optimal_linear
      } else {
      JointPRS_PRS = JointPRS_meta_PRS
      }
      ```

**Note: We recommend standardize the polygenic scores (i.e., converting the scores to zero mean and unit variance) in both tuning and testing datasets before linear combination.**


## Example
The example contains EUR, EAS, AFR and GWAS summary statistics and a bim file for 500 SNPs on chromosome 1 for HDL.
- The summary statistics comes from [GLGC](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/).
- The bim file comes from the [1000 Genome Project](https://www.internationalgenome.org/data).

The following code is a demo to use the example data for section **4. JointPRS Model Implementation**:
```
conda activate JointPRS

JointPRS_path=
reference_path= ; type=1KG
bim_path=${JointPRS_path}/example_data; bim_prefix=example
outcome_path=${JointPRS_path}/example_data
param_phi=1e-04
chr=1

pop1=EUR; pop2=EAS; pop3=AFR; pop4=SAS
r1=1; r2=1; r3=1; r4=1
sst1=${JointPRS_path}/example_data/EUR_sumstat.txt; sst2=${JointPRS_path}/example_data/EAS_sumstat.txt; sst3=${JointPRS_path}/example_data/AFR_sumstat.txt; sst4=${JointPRS_path}/example_data/SAS_sumstat.txt 
sample_size1=885546; sample_size2=116404; sample_size3=90804; sample_size4=33953

python ${JointPRS_path}/JointPRS.py \
--ref_dir=${reference_path}/${type} \
--bim_prefix=${bim_path}/${bim_prefix} \
--pop=${pop1},${pop2},${pop3},${pop4} \
--rho_cons=${r1},${r2},${r3},${r4} \
--sst_file=${sst1},${sst2},${sst3},${sst4} \
--n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4} \
--chrom=${chr} \
--phi=${param_phi} \
--out_dir=${outcome_path} \
--out_name=JointPRS_tune_${pop1}_${pop2}_${pop3}_${pop4}_${r1}${r2}${r3}${r4}_${type}
```
Here you still need to specify `JointPRS_path` and `reference_path` by yourself. For further implementation of the JointPRS data-adaptive approach when the individual-level tuning data is available, please use your own tuning data and follow the pipeline described in section **5. JointPRS Data-Adaptive Approach**.

## Acknowledgment
Part of the code is adapted from [PRS-CSx](https://github.com/getian107/PRScsx/tree/master). We thank Dr. Tian Ge for sharing his code and LD reference panels.

## Support
Please direct any problems or questions to Leqi Xu (leqi.xu@yale.edu).

## Citation
Xu L, Zhou G, Jiang W, Zhang H, Dong Y, Guan L, Zhao H. JointPRS: A Data-Adaptive Framework for Multi-Population Genetic Risk Prediction Incorporating Genetic Correlation. bioRxiv. 2024 Sept 1:2023-10. (https://www.biorxiv.org/content/10.1101/2023.10.29.564615v4)

