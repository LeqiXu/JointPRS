#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for JointPRS.

"""

import math
import numpy as np
from scipy import linalg
from numpy import random

import psi_update

def mcmc(a, b, phi, rho_cons, snp_dict, beta_sumstat, frq_dict, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, chrom, out_dir, out_name, meta, seed):
    print('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    n_pst = (n_iter - n_burnin) / thin  # number of times for updating parameters
    n_pop = len(pop)  # number of ethnic groups
    n_rho = int(0.5 * n_pop * (n_pop - 1))
    p_tot = len(snp_dict['SNP'])  # total number of SNPs in all ethnic groups (union)

    p = np.zeros(n_pop, dtype=int)  # number of SNPs in each ethnic group
    n_blk = np.zeros(n_pop, dtype=int)  # number of blocks in each ethnic group
    het = {}  # variances of each SNP in each group
    for pp in range(n_pop):
        p[pp] = len(beta_sumstat[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = np.sqrt(2.0 * frq_dict[pp] * (1.0 - frq_dict[pp]))

    idx_grp = np.zeros((p_tot,n_pop),dtype=int) # index matrix indicating whether this SNP is available in this ethnic group (1) or not (0)
    nonzero_grp = {} # ethnic groups index that have this SNP
    n_grp = np.zeros(p_tot,dtype=int)  # number of ethnic groups have this SNP
    for jj in range(p_tot):
        for pp in range(n_pop):
            if jj in idx_dict[pp]:
                idx_grp[jj,pp] = 1
                n_grp[jj] += 1
        nonzero_grp[jj] = np.array(np.nonzero(idx_grp[jj,:])).flatten().astype('int')

    uni_idx_grp = np.unique(idx_grp,axis=0) # unique number of index patterns
    n_uni_type = np.shape(uni_idx_grp)[0]
    uni_nonzero_grp = {} # unique non zero groups
    n_uni_grp = np.zeros(n_uni_type,dtype=int)
    for gg in range(n_uni_type):
        uni_nonzero_grp[gg] = np.array(np.nonzero(uni_idx_grp[gg,:])).flatten().astype('int')
        n_uni_grp[gg] = len(uni_nonzero_grp[gg])

    map_nonzero_grp = np.zeros(p_tot, dtype=int) # index of this SNP pattern in the unique nonzero group
    for jj in range(p_tot):
        for key, val in uni_nonzero_grp.items():
            if np.array_equal(val,nonzero_grp[jj]):
                map_nonzero_grp[jj] = key

    # initialization
    n = np.array(n)
    n_sqrt = np.array(np.sqrt(n))  # squared root sample size in all ethnic groups

    beta = {} # true effect dict for all ethnic groups
    beta_all = np.zeros((p_tot,n_pop)) # enlarge the beta to a union array
    for pp in range(n_pop):
        beta[pp] = np.zeros(p[pp]) # true effect for SNPs in each ethnic group

    sigma = np.ones(n_pop) # variance vector in all ethnic groups

    map_rho = np.zeros((n_rho,2), dtype=int)
    rho_row_idx = 0
    for pp1 in range(n_pop-1):
        for pp2 in range(pp1+1,n_pop):
            map_rho[rho_row_idx,0] = pp1
            map_rho[rho_row_idx,1] = pp2
            rho_row_idx += 1

    rho_share_idx = {}
    rho_share_len = np.zeros(n_rho, dtype=int)
    for rr in range(n_rho):
        pp1 = map_rho[rr,0]
        pp2 = map_rho[rr,1]
        set_idx1 = set(idx_dict[pp1]); set_idx2 = set(idx_dict[pp2])
        rho_share_idx[rr] = list(set_idx1.intersection(set_idx2))
        rho_share_len[rr] = len(rho_share_idx[rr])

    cov = {}  # covariance matrix dict for all types of non zero groups
    inv_cov = {} # inverse of the covariance matrix dict for all types of non zero groups
    m_inv_cov = {} # sample size multiplies inverse of the covariance matrix dict for all types of non zero groups
    inv_cov_all = {} # inverse of the covariance for each SNP with all groups and 0 represents missing part
    for gg in range(n_uni_type):
        cov[gg] = np.identity(n_uni_grp[gg])
        inv_cov[gg] = np.identity(n_uni_grp[gg])
        group = uni_nonzero_grp[gg]
        m_inv_cov[gg] = np.diag(n[group])
    for pp in range(n_pop):
        inv_cov_all[pp] = np.zeros((p_tot,n_pop))
        inv_cov_all[pp][:,pp] = 1
    rho = np.identity(n_pop)

    psi = np.ones(p_tot)  # local marker specific parameter

    if phi == None:  # global shared scaling parameter
        phi = 1.0
        phi_updt = True
    else:
        phi_updt = False

    # space allocation
    beta_est = {}  # the average of n_pst estimated SNP effects
    beta_sq_est = {}  # the average of n_pst estimated SNP effects variances
    sigma_est = np.zeros(n_pop) # the average of n_pst estimated variances
    for pp in range(n_pop):
        beta_est[pp] = np.zeros(p[pp])
        beta_sq_est[pp] = np.zeros(p[pp])
        
    psi_est = np.zeros(p_tot)  # the average of n_pst estimated local marker specific parameters
    phi_est = 0.0  # the average of n_pst estimated global shared scaling parameters
    rho_est = np.zeros((n_pop,n_pop))

    # MCMC
    quad = np.zeros(n_pop)  # for an ethnic group, its beta^T * LD matrix * its beta
    err = np.zeros(n_pop) # for an ethnic group, its sample size n / 2.0 * (1.0 - 2.0 * sum(beta * beta_sumstat) + quad)

    accept_rho = np.zeros(n_rho)
    rho_step = 0.02 * np.ones(n_rho)

    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('--- iter-' + str(itr) + ' ---')
        
        for pp in range(n_pop):
            mm = 0; quad[pp] = 0.0
            psi_pp = psi[idx_dict[pp]]
            beta_all_pp = beta_all[idx_dict[pp],:]
            inv_cov_all_pp = inv_cov_all[pp][idx_dict[pp],:]
            n_sqrt_pp = np.delete(n_sqrt,pp)
            for kk in range(n_blk[pp]):  # number of LD blocks in that population
                if blk_size[pp][kk] == 0:
                    continue
                else:                        
                    idx_blk = range(mm,mm+blk_size[pp][kk])
                    inv_cov_psi_blk = inv_cov_all_pp[idx_blk,:] / psi_pp[idx_blk][:,None]
                    dinvt = ld_blk[pp][kk] + sigma[pp] * np.diag(inv_cov_psi_blk[:,pp])
                    dinvt_chol = linalg.cholesky(dinvt)

                    a_matrix = np.multiply(np.delete(inv_cov_psi_blk, pp, 1), np.delete(beta_all_pp[idx_blk,:],pp,1))
                    beta_calculate = beta_sumstat[pp][idx_blk] \
                        - sigma[pp] / n_sqrt[pp] * a_matrix @ n_sqrt_pp    
                    beta_tmp = linalg.solve_triangular(dinvt_chol, beta_calculate, trans='T') \
                               + np.sqrt(sigma[pp]) / n_sqrt[pp] * random.randn(len(idx_blk))
                    beta[pp][idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')

                    quad[pp] += beta[pp][idx_blk].T @ ld_blk[pp][kk] @ beta[pp][idx_blk]
                    mm += blk_size[pp][kk]

            beta_all[idx_dict[pp],pp] = beta[pp]        
            err[pp] = max(n[pp] / 2.0 * (1.0 - 2.0 * sum(beta[pp] * beta_sumstat[pp]) + quad[pp] + sum(beta[pp] ** 2/psi_pp)), n[pp] / 2.0 * sum(beta[pp] ** 2/psi_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err[pp])

        for rr in range(n_rho):
            pp1 = map_rho[rr,0]
            pp2 = map_rho[rr,1]
            rho_star = np.random.uniform(rho[pp1,pp2] - rho_step[rr], rho[pp1,pp2] + rho_step[rr], size=1)[0]
            rho_constraint = rho_cons[pp1] * rho_cons[pp2]
            if rho_star >= 0 and rho_star <= min(rho_constraint,0.99):
                rho_new = 1.0 * rho
                rho_new[pp1,pp2] = 1.0 * rho_star
                rho_new[pp2,pp1] = 1.0 * rho_star

                if np.all(linalg.eigvals(rho_new) > 0):
                    psigma_star = rho_star * np.sqrt(sigma[pp1] * sigma[pp2])
                    psigma_original = rho[pp1,pp2] * np.sqrt(sigma[pp1] * sigma[pp2])
                    beta_share1 = beta_all[rho_share_idx[rr],pp1]
                    beta_share2 = beta_all[rho_share_idx[rr],pp2]
                    psi_share = psi[rho_share_idx[rr]]

                    log_t_star = - 0.5 * rho_share_len[rr] * math.log(sigma[pp1] * sigma[pp2] - psigma_star ** 2) \
                        - 0.5 * n[pp1] * n[pp2] / (sigma[pp1] * sigma[pp2] - psigma_star ** 2) \
                            * (sum(beta_share1 ** 2 / psi_share) * sigma[pp2] / n[pp2] \
                            - 2 * sum(beta_share1 * beta_share2 / psi_share) * psigma_star / (n_sqrt[pp1] * n_sqrt[pp2]) \
                                + sum(beta_share2 ** 2 / psi_share) * sigma[pp1] / n[pp1])
                    log_t = - 0.5 * rho_share_len[rr] * math.log(sigma[pp1] * sigma[pp2] - psigma_original ** 2) \
                        - 0.5 * n[pp1] * n[pp2] / (sigma[pp1] * sigma[pp2] - psigma_original ** 2) \
                            * (sum(beta_share1 ** 2 / psi_share) * sigma[pp2] / n[pp2] \
                            - 2 * sum(beta_share1 * beta_share2 / psi_share) * psigma_original / (n_sqrt[pp1] * n_sqrt[pp2]) \
                                + sum(beta_share2 ** 2 / psi_share) * sigma[pp1] / n[pp1])
                                                    
                    log_ratio_t = log_t_star - log_t

                    if log_ratio_t >= 0:
                        accept_rho[rr] += 1
                        rho = rho_new
                    elif math.exp(log_ratio_t) >= np.random.uniform(low=0.0, high=1.0, size=1)[0]:
                        accept_rho[rr] += 1
                        rho = rho_new
                    else: 
                        accept_rho[rr] += 0

        for gg in range(n_uni_type):
            for gg1 in range(n_uni_grp[gg]):
                cov[gg][gg1,gg1] = sigma[uni_nonzero_grp[gg][gg1]]
                for gg2 in range(gg1+1,n_uni_grp[gg]):
                    cov[gg][gg1,gg2] = rho[uni_nonzero_grp[gg][gg1],uni_nonzero_grp[gg][gg2]] * np.sqrt(sigma[uni_nonzero_grp[gg][gg1]] * sigma[uni_nonzero_grp[gg][gg2]])
                    cov[gg][gg2,gg1] = cov[gg][gg1,gg2]
            inv_cov[gg] = np.linalg.inv(cov[gg])
            group = uni_nonzero_grp[gg]
            m_inv_cov[gg] = np.diag(n_sqrt[group]) @ inv_cov[gg] @ np.diag(n_sqrt[group])

        for jj in range(p_tot):
            gg = map_nonzero_grp[jj]
            inv_cov_jj = inv_cov[gg]
            for ss in range(n_uni_grp[gg]):
                for tt in range(n_uni_grp[gg]):
                    inv_cov_all[uni_nonzero_grp[gg][ss]][jj,uni_nonzero_grp[gg][tt]] = inv_cov_jj[ss,tt]
       
        delta = random.gamma(a+b, 1.0/(psi+phi))
        xx = np.zeros(p_tot)

        psi = psi_update.psi_update(p_tot, map_nonzero_grp, uni_nonzero_grp, m_inv_cov, beta_all, 
                                    xx, psi, a, n_grp, delta)    
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p_tot*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst

            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst
            rho_est = rho_est + rho/n_pst

    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    # meta
    if meta == 'TRUE':
        vv = np.zeros(p_tot)
        zz = np.zeros(p_tot)
        for pp in range(n_pop):
            vv[idx_dict[pp]] += 1.0/(beta_sq_est[pp]-beta_est[pp]**2)
            zz[idx_dict[pp]] += 1.0/(beta_sq_est[pp]-beta_est[pp]**2)*beta_est[pp]
        mu = zz/vv

    # write posterior effect sizes
    for pp in range(n_pop):
        if phi_updt == True:
            eff_file = out_dir + '/' + '%s_%s_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (out_name, pop[pp], a, b, chrom)
        else:
            eff_file = out_dir + '/' + '%s_%s_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (out_name, pop[pp], a, b, phi, chrom)

        snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
        bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
        a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
        a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]
        psi_est_pp = [psi_est[ii] for ii in idx_dict[pp]]

        with open(eff_file, 'w') as ff:
            for snp, bp, a1, a2, beta, psi in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp], psi_est_pp):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\t%.6e\n' % (chrom, snp, bp, a1, a2, beta, psi))

    # write correlation matrix
    if phi_updt == True:
        eff_file = out_dir + '/' + '%s_pst_corr_a%d_b%.1f_phiauto_chr%d.txt' % (out_name, a, b, chrom)
    else:
        eff_file = out_dir + '/' + '%s_pst_corr_a%d_b%.1f_phi%1.0e_chr%d.txt' % (out_name, a, b, phi, chrom)
    np.savetxt(eff_file,rho_est)

    if meta == 'TRUE':
        if phi_updt == True:
            eff_file = out_dir + '/' + '%s_META_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (out_name, a, b, chrom)
        else:
            eff_file = out_dir + '/' + '%s_META_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (out_name, a, b, phi, chrom)

        with open(eff_file, 'w') as ff:
            for snp, bp, a1, a2, beta in zip(snp_dict['SNP'], snp_dict['BP'], snp_dict['A1'], snp_dict['A2'], mu):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    print('... Estimated Correlation Matrix for chrom %d is : ...' % chrom)
    print(rho_est)
    print('... Done ...')