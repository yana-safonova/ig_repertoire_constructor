import pickle


def load_model(filename):
    with open(filename, 'rb') as f:
        model = pickle.load(f)
    return model


def full_igh_model():
    return load_model('shm_kmer_model_estimator/full_igh_model.dump')


# This one lacks samples_age because age contains only IGH
def full_all_chains_model():
    return load_model('shm_kmer_model_estimator/full_all_chains_model.dump')


def full_igh_model_extended_base():
    return load_model('shm_kmer_model_estimator/full_igh_model_extended_base.dump')


def dump_model(filename, tuple_datasets, chains=None, strategies=None):
    import shm_kmer_model_estimator
    estimator = shm_kmer_model_estimator.ShmKmerModelEstimator()
    model = estimator.estimate_models_of_one_type(tuple_datasets, chains=chains, strategies=strategies)
    with open(filename, 'wb') as f:
        pickle.dump(model, f)


def dump_full_igh_model():
    import sample_reader.standard_samples as standard_samples
    dump_model(filename='shm_kmer_model_estimator/full_igh_model.dump',
               tuple_datasets=(standard_samples.kmer_freq_matrices_fv,
                               standard_samples.kmer_freq_matrices_gmc,
                               standard_samples.kmer_freq_matrices_ido,
                               standard_samples.kmer_freq_matrices_age,
                               standard_samples.kmer_freq_matrices_paired),
               chains=['IGH'])


def dump_full_all_chains_model():
    import sample_reader.standard_samples as standard_samples
    dump_model(filename='shm_kmer_model_estimator/full_all_chains_model.dump',
               tuple_datasets=(standard_samples.kmer_freq_matrices_fv,
                               standard_samples.kmer_freq_matrices_gmc,
                               standard_samples.kmer_freq_matrices_ido,
                               standard_samples.kmer_freq_matrices_paired))


def dump_full_igh_model_extended_base():
    import sample_reader.standard_samples as standard_samples
    dump_model(filename='shm_kmer_model_estimator/full_igh_model_extended_base.dump',
               tuple_datasets=(standard_samples.kmer_freq_matrices_extended_fv,
                               standard_samples.kmer_freq_matrices_extended_gmc,
                               standard_samples.kmer_freq_matrices_extended_ido,
                               standard_samples.kmer_freq_matrices_extended_age,
                               standard_samples.kmer_freq_matrices_extended_paired),
               chains=['IGH'])


def dump_all_models():
    dump_full_all_chains_model()
    dump_full_igh_model()
    dump_full_igh_model_extended_base()
