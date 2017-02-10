from __future__ import division

import sys
sys.path.insert(0, "..")

import os
from glob import glob
from collections import OrderedDict

import numpy as np
import pandas as pd

from special_utils.os_utils import smart_makedirs

from config.config import config, read_config
from config.parse_input_args import parse_args

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies


def read_models(model_dir):
    models_path = glob(os.path.join(model_dir, "*.csv"))
    models = OrderedDict()

    for model_path in models_path:
        model_filename = os.path.basename(model_path)
        model_filename = model_filename.split('.')[0]
        strategy, chain = model_filename.split('_')
        strategy, chain = MutationStrategies[strategy], Chains[chain]
        models[(strategy, chain)] = pd.read_csv(model_path)

    return models


def coverage_kmers_utils_all_models(models, chains=None, strategies=None, verbose=False):
    if chains is None:
        chains = Chains
    if strategies is None:
        strategies = MutationStrategies

    from pprint import pprint
    res = OrderedDict()

    for strategy in strategies:
        for chain in chains:
            res[(strategy, chain)] = \
                coverage_kmers_utils(models[strategy][chain], strategy, chain)
            if verbose:
                pprint("Strategy: %s, Chain: %s" %(strategy.name, chain.name))
                pprint(res[strategy][chain])
                print("")
    return pd.DataFrame(res, index=res[(MutationStrategies.NoKNeighbours, Chains.IGH)].keys())


def coverage_kmers_utils(model, strategy, chain):
    from kmer_utilities.kmer_utilities import kmer_names
    kmer_names = np.array(kmer_names())

    def nonnan_ind(model, column_names, kmer_names=kmer_names):
        columns = [np.array(model[x]) for x in column_names]
        nonnan_ind = [~np.isnan(x)[:,np.newaxis] for x in columns]
        nonnan_ind = np.concatenate(nonnan_ind, axis=1)
        nonnan_ind = np.prod(nonnan_ind, axis=1, dtype=np.bool)
        return set(kmer_names[nonnan_ind])

    if chain == Chains.IG:
        return

    from genomic_kmers.genomic_kmers import get_genomic_kmers

    good_sp_kmers_fr = nonnan_ind(model, ["start_point_beta_FR_shape1",
                                          "start_point_beta_FR_shape2"])


    good_sp_kmers_cdr = nonnan_ind(model, ["start_point_beta_CDR_shape1",
                                           "start_point_beta_CDR_shape2"])


    good_sp_kmers_subst = nonnan_ind(model, ["start_point_dir_shape1",
                                             "start_point_dir_shape2",
                                             "start_point_dir_shape3"])

    gen_kmers_fr, gen_kmers_cdr = get_genomic_kmers(chain)
    gen_kmers_fr, gen_kmers_cdr = set(gen_kmers_fr), set(gen_kmers_cdr)
    gen_kmers_all = gen_kmers_fr | gen_kmers_cdr


    est_fr = nonnan_ind(model, ["beta_FR_shape1",
                                "beta_FR_shape2"])

    est_cdr = nonnan_ind(model, ["beta_CDR_shape1",
                                 "beta_CDR_shape2"])

    est_subst = nonnan_ind(model, ["dir_shape1",
                                   "dir_shape2",
                                   "dir_shape3"])
    so_fr = np.array(model.success_optim_beta_FR, dtype=np.bool)
    so_cdr = np.array(model.success_optim_beta_CDR, dtype=np.bool)
    so_subst = np.array(model.success_optim_dir, dtype=np.bool)
    so_fr, so_cdr, so_subst = kmer_names[so_fr], kmer_names[so_cdr], kmer_names[so_subst]
    so_fr, so_cdr, so_subst = set(so_fr), set(so_cdr), set(so_subst)

    assert good_sp_kmers_fr <= gen_kmers_fr
    assert good_sp_kmers_cdr <= gen_kmers_cdr
    assert good_sp_kmers_subst <= gen_kmers_all

    assert est_fr == good_sp_kmers_fr
    assert est_cdr == good_sp_kmers_cdr
    assert est_subst <= good_sp_kmers_subst
    return OrderedDict([# '# est FR / # good start point FR': len(est_fr) / len(good_sp_kmers_fr),
                # '# est CDR / # good start point CDR': len(est_cdr) / len(good_sp_kmers_cdr),
                # '# est subst / # good start point subst': len(est_subst) / len(good_sp_kmers_subst),
                # '# (good sp FR) - (est FR)': len(good_sp_kmers_fr - est_fr),
                # '# (good sp CDR) - (est CDR)': len(good_sp_kmers_cdr - est_cdr),
                # '# (good sp subst) - (est subst)': len(good_sp_kmers_subst - est_subst),
                ('# genomic fr', len(gen_kmers_fr)),
                ('# genomic cdr', len(gen_kmers_cdr)),
                ('# genomic subst (fr + cdr)', len(gen_kmers_all)),
                ('# well estimated kmers fr', len(est_fr & so_fr)),
                ('# well estimated kmers cdr', len(est_cdr & so_cdr)),
                ('# well estimated kmers fr or cdr', len((est_fr & so_fr) | (est_cdr & so_cdr))),
                ('# well estimated kmers fr and cdr', len((est_fr & so_fr) & (est_cdr & so_cdr))),
                ('# well estimated kmers subst', len(est_subst & so_subst)),
                ('# genomic - well estimated kmers fr', len(gen_kmers_fr - (est_fr & so_fr))),
                ('# genomic - well estimated kmers cdr', len(gen_kmers_cdr - (est_cdr & so_cdr))),
                ('# genomic - well estimated kmers subst', len(gen_kmers_all - (est_subst & so_subst))),
                ('# genomic - good start point fr', len(gen_kmers_fr - good_sp_kmers_fr)),
                ('# genomic - good start point cdr', len(gen_kmers_cdr - good_sp_kmers_cdr)),
                ('# genomic - good start point subst', len(gen_kmers_all - good_sp_kmers_subst)),
                # ('# good start point FR / # gen FR', len(good_sp_kmers_fr) / len(gen_kmers_fr)),
                # ('# good start point CDR / # gen CDR', len(good_sp_kmers_cdr) / len(gen_kmers_cdr)),
                # ('# good start point subst / #gen subst', len(good_sp_kmers_subst) / len(gen_kmers_all)),
                # ('# kmers in (gen FR) - (good sp FR)', len(gen_kmers_fr - good_sp_kmers_fr)),
                # ('# kmers in (gen CDR) - (good sp CDR)', len(gen_kmers_cdr - good_sp_kmers_cdr)),
                # ('# kmers in (gen subst) - (good sp subst)', len(gen_kmers_all - good_sp_kmers_subst)),
                # ('# kmers optim not successful that est, FR', len(est_fr - so_fr)),
                # ('# kmers optim not successful that est, CDR', len(est_cdr - so_cdr)),
                # ('# kmers optim not successful that est, subst', len(est_subst - so_subst)),
                # ('# kmers optim not successful and good est that gen, FR', len(gen_kmers_fr - so_fr)),
                # ('# kmers optim not successful and good est that gen, CDR', len(gen_kmers_cdr - so_cdr)),
                # ('# kmers optim not successful and good est that gen, subst', len(gen_kmers_all - so_subst)),
                # ('% kmers optim not successful and good est that gen, FR', len(gen_kmers_fr - (so_fr & est_fr)) / len(gen_kmers_fr)),
                # ('% kmers optim not successful and good est that gen, CDR', len(gen_kmers_cdr - (so_cdr & est_cdr)) / len(gen_kmers_cdr)),
                # ('% kmers optim not successful and good est that gen, subst', len(gen_kmers_all - (so_subst & est_subst)) / len(gen_kmers_all))
            ])


if __name__ == '__main__':
    input_config = read_config(parse_args().input)
    models = read_models(input_config.kmer_model_estimating.outdir)
    chains = [Chains.IGH, Chains.IGK, Chains.IGL]
    mut_str = [MutationStrategies.NoKNeighbours]
    model_config = input_config.kmer_model_estimating
    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    outfile = os.path.join(outdir, model_config.kmer_coverage_analysis)
    smart_makedirs(outdir)
    apply_to_each_model(
        models,
        verbose=False,
        chains=chains,
        strategies=mut_str).to_csv(outfile, sep=',')
