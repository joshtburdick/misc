import jax
import jax.numpy as jnp
from flax import linen as nn
from typing import Callable

class GNNLayer(nn.Module):
    """A simple Message Passing layer."""
    features: int
    activation: Callable = nn.relu

    @nn.compact
    def __call__(self, x, adj):
        # x: (nodes, features)
        # adj: (nodes, nodes)
        
        # Message passing: sum of neighbor features
        messages = jnp.matmul(adj, x)
        
        # Combine self and neighbor information
        h_neighbor = nn.Dense(self.features, name="neighbor_dense")(messages)
        h_self = nn.Dense(self.features, name="self_dense")(x)
        
        return self.activation(h_neighbor + h_self)

class CliqueModel(nn.Module):
    """
    A GNN for graph-level classification.
    Suitable for both detection and parity tasks.
    """
    num_layers: int = 3
    hidden_dim: int = 64
    out_dim: int = 1  # Single logit for binary classification

    @nn.compact
    def __call__(self, adj):
        # adj: (batch, nodes, nodes)
        batch_size, num_nodes, _ = adj.shape
        
        # Initial node features: just the degree of each node
        # (batch, nodes, 1)
        x = jnp.sum(adj, axis=-1, keepdims=True)
        
        # Message passing layers
        for i in range(self.num_layers):
            # Dense layers handle (batch, nodes, features) naturally.
            # jnp.matmul(adj, x) handles (batch, nodes, nodes) @ (batch, nodes, features).
            x = GNNLayer(features=self.hidden_dim, name=f"gnn_layer_{i}")(x, adj)
            
        # Global Pooling (Sum pooling)
        # (batch, hidden_dim)
        graph_repr = jnp.sum(x, axis=1)
        
        # Final MLP Head
        h = nn.Dense(self.hidden_dim)(graph_repr)
        h = nn.relu(h)
        h = nn.Dense(self.hidden_dim)(h)
        h = nn.relu(h)
        logits = nn.Dense(self.out_dim)(h)
        
        return jnp.squeeze(logits, axis=-1)
