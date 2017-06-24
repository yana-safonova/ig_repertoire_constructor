import os

import numpy as np
from scipy.stats import mannwhitneyu

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from config.config import config, read_config
from config.parse_input_args import parse_args

from shm_kmer_model import cab_shm_model, yale_model, cab_shm_nonproductive_model
from mutation_strategies import mutation_strategies
from chains.chains import Chains

from shm_kmer_model.cab_shm_model import Region
from likelihood_calculator.likelihood_calculator import LikelihoodCalculator

from tree_test_utilities import flu_trees_statistics_calculator, tree_test_utilities

from special_utils.os_utils import smart_mkdir


def run_tree_test_chain_strategy(strategy, chain_type, log_dir):
    ymodel                  = yale_model.YaleSHM_Model()
    cabmodel                = cab_shm_model.CAB_SHM_Model(strategy, chain_type)
    cab_nonproductive_model = cab_shm_nonproductive_model.CAB_SHM_Nonproductive_Model(strategy, chain_type)

    flu_trees_paths = flu_trees_statistics_calculator.get_flu_trees_paths(chain_type=chain_type)
    tester = tree_test_utilities.TreeTester()
    yresults                  = tester.get_likelihood_statistics_trees(model=ymodel,
                                                                       tree_paths=flu_trees_paths)
    cabresults                = tester.get_likelihood_statistics_trees(model=cabmodel,
                                                                       tree_paths=flu_trees_paths)
    cab_nonproductive_results = tester.get_likelihood_statistics_trees(model=cab_nonproductive_model,
                                                                       tree_paths=flu_trees_paths)

    filename = '%s_%s.' % (strategy.name, chain_type.name)
    log_filename = os.path.join(log_dir, filename + 'log')

    smart_mkdir(log_dir)
    with open(log_filename, 'w') as f:
        f.write("Median for Yale: %lf\n" % np.nanmedian(yresults.accuracies))
        f.write("Median for CAB: %lf\n" % np.nanmedian(cabresults.accuracies))
        f.write("Median for CAB nonproductive: %lf\n" % np.nanmedian(cab_nonproductive_results.accuracies))

        f.write("Yale Accuracies: %s\n" % yresults.accuracies)
        f.write("CAB Accuracies: %s\n" % cabresults.accuracies)
        f.write("CAB nonproductive Accuracies: %s\n" % cab_nonproductive_results.accuracies)

        f.write("Full accuracy for Yale: %lf\n" % yresults.full_accuracy)
        f.write("Full accuracy for CAB: %lf\n" % cabresults.full_accuracy)
        f.write("Full accuracy for CAB nonproductive: %lf\n" % cab_nonproductive_results.full_accuracy)

        mw_ycab = mannwhitneyu(yresults.accuracies, cabresults.accuracies, use_continuity=True)
        f.write("Mann Whitney test for equality of means Yale vs CAB: %s\n" % str(mw_ycab))

        mw_ycabn = mannwhitneyu(yresults.accuracies, cab_nonproductive_results.accuracies, use_continuity=True)
        f.write("Mann Whitney test for equality of means Yale vs CAB nonprod: %s\n" % str(mw_ycabn))

        mw_cabcabn = mannwhitneyu(cabresults.accuracies, cab_nonproductive_results.accuracies, use_continuity=True)
        f.write("Mann Whitney test for equality of means CAB vs CAB nonprod: %s\n" % str(mw_cabcabn))

    def draw(title, yale, cab, cabn, fig_filename, bins=10):
        def draw_distplot(x, color):
            return sns.distplot(x[~np.isnan(x)], bins=bins, rug=True, kde=True)
        ypict = draw_distplot(yale.accuracies, color='blue')
        cabpict = draw_distplot(cab.accuracies, color='green')
        cabnpict = draw_distplot(cabn.accuracies, color='red')
        plt.legend(['Yale', 'CAB', 'CAB nonprod'])
        ax = plt.gca()
        leg = ax.get_legend()
        leg.legendHandles[0].set_color('blue')
        leg.legendHandles[1].set_color('green')
        leg.legendHandles[2].set_color('red')
        plt.title(title)
        plt.savefig(fig_filename, format='pdf')
        plt.clf()

    figures_dir = os.path.join(log_dir, 'figures')
    smart_mkdir(figures_dir)
    draw('Accuracy', yresults, cabresults, cab_nonproductive_results,
          os.path.join(figures_dir, filename) + "pdf")
    return {'Yale': yresults, 'CAB': cabresults, 'CAB_nonprod': cab_nonproductive_results}


def main():
    parsed_args = parse_args()
    input_config = read_config(parsed_args.input)
    test_config = input_config.kmer_model_tree_test
    for strategy in mutation_strategies.MutationStrategies:
        for chain_type in Chains:
            if chain_type.name == 'IG':
                continue
            print(strategy, chain_type.name)
            run_tree_test_chain_strategy(strategy, chain_type, test_config.outdir)



if __name__ == "__main__":
    main()
