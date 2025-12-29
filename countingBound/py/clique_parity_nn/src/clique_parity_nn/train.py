import jax
import jax.numpy as jnp
import optax
import argparse
import os
from flax.training import train_state
from clique_parity_nn.models import CliqueModel
from clique_parity_nn.data import create_dataset
from clique_parity_nn.plotting import plot_training_history, plot_k_scaling
import numpy as np

def compute_loss(params, model, batch_adj, batch_labels):
    logits = model.apply({'params': params}, batch_adj)
    # Binary cross entropy loss
    loss = optax.sigmoid_binary_cross_entropy(logits, batch_labels).mean()
    return loss

def compute_accuracy(params, model, batch_adj, batch_labels):
    logits = model.apply({'params': params}, batch_adj)
    preds = (logits > 0).astype(jnp.float32)
    return jnp.mean(preds == batch_labels)

def train_task(task_name, model, params, train_data, test_data, num_epochs=100, learning_rate=1e-3):
    """Trains the model on a specific task and returns accuracy history."""
    tx = optax.adam(learning_rate)
    state = train_state.TrainState.create(
        apply_fn=model.apply,
        params=params,
        tx=tx
    )
    
    train_adj, train_labels = train_data
    test_adj, test_labels = test_data
    
    @jax.jit
    def train_step(state, batch_adj, batch_labels):
        grad_fn = jax.value_and_grad(compute_loss)
        loss, grads = grad_fn(state.params, model, batch_adj, batch_labels)
        state = state.apply_gradients(grads=grads)
        return state, loss

    history = []
    for epoch in range(num_epochs):
        state, loss = train_step(state, train_adj, train_labels)
        test_acc = float(compute_accuracy(state.params, model, test_adj, test_labels))
        history.append(test_acc)
        
        if epoch % 20 == 0 or epoch == num_epochs - 1:
            print(f"  {task_name} Epoch {epoch:3d}: Loss {loss:.4f}, Test Acc {test_acc:.4f}")
            
    return state, history

def run_experiment(n, k_range, num_samples, num_epochs, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    results = []

    for k in k_range:
        print(f"\nEvaluating k={k} (n={n}, samples={num_samples})...")
        adjs, exists, parity = create_dataset(num_samples, n, k)
        
        # Split data
        split = int(0.8 * num_samples)
        train_adj, test_adj = adjs[:split], adjs[split:]
        train_exists, test_exists = exists[:split], exists[split:]
        train_parity, test_parity = parity[:split], parity[split:]
        
        # Initialize Model and Params once (fairness)
        # Using 3 layers as a baseline, or k
        num_layers = max(3, k)
        model = CliqueModel(num_layers=num_layers)
        key = jax.random.PRNGKey(0)
        init_params = model.init(key, train_adj[:1])['params']
        
        # Task 1: Detection
        _, history_exists = train_task(
            "Detection", model, init_params, 
            (train_adj, train_exists), (test_adj, test_exists),
            num_epochs=num_epochs
        )
        
        # Task 2: Parity
        _, history_parity = train_task(
            "Parity", model, init_params, 
            (train_adj, train_parity), (test_adj, test_parity),
            num_epochs=num_epochs
        )

        results.append({
            'k': k,
            'final_test_acc_exists': history_exists[-1],
            'final_test_acc_parity': history_parity[-1]
        })

        # Plot training history for this specific k
        plot_training_history(history_exists, history_parity, n, k, output_dir)

    # Plot performance vs k
    plot_k_scaling(results, n, k_range, output_dir)
    print(f"\nAll plots saved to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Clique Detection vs Clique Parity")
    parser.add_argument("--n", type=int, default=15, help="Number of nodes in the graph")
    parser.add_argument("--k_min", type=int, default=3, help="Min clique size to evaluate")
    parser.add_argument("--k_max", type=int, default=5, help="Max clique size to evaluate")
    parser.add_argument("--samples", type=int, default=2000, help="Number of samples per k")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    parser.add_argument("--out", type=str, default="results", help="Output directory for plots")
    
    args = parser.parse_args()
    k_range = list(range(args.k_min, args.k_max + 1))
    
    run_experiment(args.n, k_range, args.samples, args.epochs, args.out)
