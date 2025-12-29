import networkx as nx
import numpy as np
import jax.numpy as jnp
from typing import Tuple, List

def get_clique_stats(G: nx.Graph, k: int) -> Tuple[bool, int]:
    """
    Computes ground truth for clique detection and clique parity.
    
    Args:
        G: NetworkX graph.
        k: Size of the clique to look for.
        
    Returns:
        exists: True if at least one clique of size k exists.
        parity: Parity of the number of cliques of size k (count % 2).
    """
    # enumerate_all_cliques finds cliques by size in increasing order.
    # It is efficient for reasonably sized graphs.
    count = 0
    exists = False
    
    for clique in nx.enumerate_all_cliques(G):
        if len(clique) == k:
            exists = True
            count += 1
        elif len(clique) > k:
            # Since enumerate_all_cliques yields in increasing order of size,
            # we can stop once we see cliques larger than k.
            break
            
    return exists, count % 2

def generate_random_graph_data(n: int, p: float, k: int):
    """Generates a single random graph and its labels."""
    G = nx.fast_gnp_random_graph(n, p)
    exists, parity = get_clique_stats(G, k)
    
    # Convert graph to adjacency matrix for the NN
    adj = nx.to_numpy_array(G)
    return adj, exists, parity

def create_dataset(num_samples: int, n: int, k: int, p_range: Tuple[float, float] = (0.2, 0.8)):
    """Creates a dataset of graphs with detection and parity labels."""
    adjs = []
    exists_labels = []
    parity_labels = []
    
    for _ in range(num_samples):
        p = np.random.uniform(*p_range)
        adj, exists, parity = generate_random_graph_data(n, p, k)
        adjs.append(adj)
        exists_labels.append(float(exists))
        parity_labels.append(float(parity))
        
    return (jnp.array(adjs), 
            jnp.array(exists_labels), 
            jnp.array(parity_labels))
