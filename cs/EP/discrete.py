# Expectation Propagation

import numpy




def init_vars(vars, factors):
    """Initializes

    vars: the variables, as a dict keyed by name, and containing
        the size of each variable
    factors: the factors, as a dict keyed by name,
        each of which consists of a dict with keys:
            vars: names of the m variables
            f: the factor, as a tensor with m dimensions
    Returns a dict with:
        x: the initial beliefs (uniform random, in logits)
        from_factor: the messages from the factors, as a dict
            keyed by (var_name, factor_name)
    """
    x = {name: np.logsumexp(random(size))
        for name, size in vars.items
    }
    def 




    return {"x": x,
        "from_factor": {



