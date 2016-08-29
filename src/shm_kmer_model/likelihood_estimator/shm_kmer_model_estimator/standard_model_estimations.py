import pickle

def full_igh_model():
    with open('shm_kmer_model_estimator/full_igh_model.dump', 'rb') as f:
        model = pickle.load(f)
    return model

# This one lacks samples_age because age contains only IGH
def full_all_chains_model():
    with open('shm_kmer_model_estimator/full_all_chains_model.dump', 'rb') as f:
        model = pickle.load(f)
    return model
