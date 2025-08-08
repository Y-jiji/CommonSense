import numpy as np

def subsample(data, nsamples, rank=1, seed=None):
    rng = np.random.default_rng(seed)

    n = data.shape[0]
    ngroups = n // rank
    assert ngroups >= nsamples

    selected_id = rng.choice(ngroups, nsamples, replace=False)
    row_ids = [id*rank + i for id in selected_id for i in range(rank)]

    data_selected = data[row_ids]
    return data_selected, selected_id

    
