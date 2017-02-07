import sys
sys.path.insert(0, "..")

import os
from glob import glob

import pandas as pd

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies

model_dir = "/Sid/abzikadze/shm_models/"


def get_shm_dict():
    x = dict.fromkeys(MutationStrategies)
    for strategy in x:
        x[strategy] = dict.fromkeys(Chains)
    return x


def read_models():
    models_path = glob(os.path.join(model_dir, "*.csv"))
    models = get_shm_dict()

    for model_path in models_path:
        model_filename = os.path.basename(model_path)
        model_filename = model_filename.split('.')[0]
        strategy, chain = model_filename.split('_')
        strategy, chain = MutationStrategies[strategy], Chains[chain]
        models[strategy][chain] = pd.read_csv(model_path)

    return models


def apply_to_each_model(models, func):
    res = get_shm_dict()
    for strategy in res:
        for chain in res[strategy]:
            res[strategy][chain] = func(models[strategy][chain])
    return res


if __name__ == '__main__':
    models = read_models()
    print(apply_to_each_model(models, lambda x: x))
