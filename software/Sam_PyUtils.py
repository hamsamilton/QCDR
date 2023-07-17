


def prefix_dict(pre,dct):
    """        
     adds a prefix to the keys of a dictionary

     Inputs: 
            pre  = the prefix to add
            dct = the dictionary to add the names to the keys
    """
    newdict= {pre + str(key): val for key, val in dct.items()}
    
    return newdict
 

def get_axis_range(axis_lim):                            
     """
     Get the range size  of an axis matplotlib lim object
     Input:
         axis_lim: a matplotlib lim object
     Output: The size of the range
     """
     rng = axis_lim[1] - axis_lim[0]
     return rng

def repl_missing_values_indict(_indict,_repldict):
    for key, value in _indict.items():
    # If the value is NaN, replace it with the corresponding automatically generated keys
        if value != value:
            _indict.update({key: _repldict[key]})

def vec2CDF(_vec):
    """
    Convert a vector of numbers into a CDF (Must be valid all positive etc.)

    Input: A vector of #

    Return:A vector of probs signifying a CDF
    """
        
    # Normalize the vector to make it a probability distribution
    prob_dist = _vec / np.sum(_vec)

    # Compute the cumulative distribution function (CDF)
    cdf = np.cumsum(prob_dist)

    return(cdf)


def bootstrap_samples(data, n_samples, sample_size):
    """
    Create bootstrap samples from the original data.
    
    Args:
        data (list): A list of numerical values.
        n_samples (int): The number of bootstrap samples to generate.
        sample_size (int): The size of each bootstrap sample.
    
    Returns:
        bootstrap (list): A list of bootstrap samples.
    """
    bootstrap = [np.random.choice(data, size=sample_size, replace=True) for _ in range(n_samples)]
    
    return bootstrap

# writing a custom function to calculate zscores as the implementation in scipi leaves something to be desired
def calc_zscore(_newval,_comparevec):
    
    # get distribution information
    _mn_comparevec = np.mean(_comparevec)
    _std_comparevec= np.mean(_comparevec)

    # calculate zscore
    _zscr = (_newval - _mn_comparevec) / _std_comparevec

    return _zscr

def estimate_p_values(data, bootstrap_means):
    """
    Estimate p-values for the original values based on the bootstrap means.
    
    Args:
        data (list): A list of numerical values.
        bootstrap_means (list): A list of means from the bootstrap samples.
    
    Returns:
        p_values (list): A list of p-values corresponding to the input values.
    """
    p_values = [np.sum(np.array(bootstrap_means) <= value) / len(bootstrap_means)  for value in data]

    return p_values


# This is a simple function designed to take a list of variables and transform it into a dictionary with the variable names as keys and their associated values as values
def varlist_2dict(_list):
    _dict = {}
    for item in _list:
        _dict.update({str(item),item})
    
    return _dict
