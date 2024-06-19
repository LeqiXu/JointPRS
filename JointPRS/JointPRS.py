#!/usr/bin/env python
import os
import sys
import getopt
import parse_genet
import mcmc_gtb

def parse_param():
    long_opts_list = ['ref_dir=', 'bim_prefix=', 'sst_file=', 'rho_cons=', 'a=', 'b=', 'phi=', 'n_gwas=', 'pop=', 
                      'n_iter=', 'n_burnin=', 'thin=', 'out_dir=', 'out_name=', 'chrom=', 'meta=', 'seed=', 'help']

    param_dict = {'ref_dir': None, 'bim_prefix': None, 'sst_file': None, 'rho_cons': None, 'a': 1, 'b': 0.5, 'phi': None, 'n_gwas': None, 'pop': None,
                  'n_iter': None, 'n_burnin': None, 'thin': 5, 'out_dir': None, 'out_name': None, 'chrom': range(1,23), 'meta': 'FALSE', 'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except:
            print('* Option not recognized.')
            print('* Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": param_dict['ref_dir'] = arg
            elif opt == "--bim_prefix": param_dict['bim_prefix'] = arg
            elif opt == "--sst_file": param_dict['sst_file'] = arg.split(',')
            elif opt == "--rho_cons": param_dict['rho_cons'] = list(map(float,arg.split(',')))
            elif opt == "--a": param_dict['a'] = float(arg)
            elif opt == "--b": param_dict['b'] = float(arg)
            elif opt == "--phi": param_dict['phi'] = float(arg)
            elif opt == "--n_gwas": param_dict['n_gwas'] = list(map(int,arg.split(',')))
            elif opt == "--pop": param_dict['pop'] = arg.split(',')
            elif opt == "--n_iter": param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": param_dict['thin'] = int(arg)
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
            elif opt == "--chrom": param_dict['chrom'] = arg.split(',')
            elif opt == "--meta": param_dict['meta'] = arg.upper()
            elif opt == "--seed": param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['ref_dir'] == None:
        print('* Please specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['bim_prefix'] == None:
        print('* Please specify the directory and prefix of the bim file for the target (validation/testing) dataset using --bim_prefix\n')
        sys.exit(2)
    elif param_dict['sst_file'] == None:
        print('* Please provide at least one summary statistics file using --sst_file\n')
        sys.exit(2)
    elif param_dict['rho_cons'] == None:
        print('* Please provide the constraint for each pop --rho_constraint\n')
    elif param_dict['n_gwas'] == None:
        print('* Please provide the sample size of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['pop'] == None:
        print('* Please specify the population of the GWAS sample using --pop\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] == None:
        print('* Please specify the prefix of the output file using --out_name\n')
        sys.exit(2)
    elif (len(param_dict['sst_file']) != len(param_dict['n_gwas']) or 
          len(param_dict['sst_file']) != len(param_dict['pop'])):
        print('* Length of sst_file, n_gwas and pop does not match\n')
        sys.exit(2)

    n_pop = len(param_dict['pop'])
    if param_dict['n_iter'] == None or param_dict['n_burnin'] == None:
        param_dict['n_iter'] = n_pop*1000
        param_dict['n_burnin'] = n_pop*500

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()
    n_pop = len(param_dict['pop'])

    if n_pop == 1:
        print('*** only %d discovery population detected ***\n' % n_pop)
    else:
        print('*** %d discovery populations detected ***\n' % n_pop)


    for chrom in param_dict['chrom']:
        print('##### process chromosome %d #####' % int(chrom))

        if os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3'):
            ref = '1kg'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3', int(chrom), ref)
        elif os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3'):
            ref = 'ukbb'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3', int(chrom), ref)

        vld_dict = parse_genet.parse_bim(param_dict['bim_prefix'], int(chrom))

        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = parse_genet.parse_sumstats(ref_dict, vld_dict, param_dict['sst_file'][pp], param_dict['pop'][pp], param_dict['n_gwas'][pp])

        ld_blk = {}
        blk_size = {}
        for pp in range(n_pop):
            ld_blk[pp], blk_size[pp] = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict[pp], param_dict['pop'][pp], int(chrom), ref)

        snp_dict, beta_dict, frq_dict, idx_dict = parse_genet.align_ldblk(ref_dict, vld_dict, sst_dict, n_pop, int(chrom))

        mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], param_dict['rho_cons'], snp_dict, beta_dict, frq_dict, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
                      param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], int(chrom),
                      param_dict['out_dir'], param_dict['out_name'], param_dict['meta'], param_dict['seed'])

        print('\n')


if __name__ == '__main__':
    main()
