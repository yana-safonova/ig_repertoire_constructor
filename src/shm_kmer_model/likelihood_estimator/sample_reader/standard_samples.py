""" This module contains standard samples from AbVitro and age datasets.
Func concatenate_samples creates a Panel that is as full as possible. """

from sample_reader import SampleReader

import pandas as pd

samples_fv = SampleReader().read('AbVitro/flu_time_course/FV/', ['25'])
samples_gmc = SampleReader().read('AbVitro/flu_time_course/GMC/', ['8'])
samples_ido = SampleReader().read('AbVitro/flu_time_course/IDO/')
samples_age = SampleReader().read('age/')
samples_paired = SampleReader().read('AbVitro/paired/')


def concatenate_samples(chain_type='IGH', strategy='NoKNeighbours'):
    assert chain_type in ['IGH', 'IGL', 'IGK']
    result = pd.concat((samples_fv[strategy][chain_type],
                        samples_ido[strategy][chain_type],
                        samples_gmc[strategy][chain_type],
                        samples_paired[strategy][chain_type]))
    if chain_type != 'IGH':
        return result
    else:
        return pd.concat((result, samples_age[strategy][chain_type]))
