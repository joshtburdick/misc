import matplotlib.pyplot as plt
import os

def plot_training_history(history_exists, history_parity, n, k, save_path):
    """Plots accuracy vs epochs for both tasks."""
    plt.figure(figsize=(10, 6))
    
    epochs = range(len(history_exists))
    plt.plot(epochs, history_exists, label="Detection (Acc)")
    plt.plot(epochs, history_parity, label="Parity (Acc)")
    
    plt.axhline(y=0.5, color='r', linestyle='--', label="Random Guess (0.5)")
    
    plt.title(f"Training Accuracy Comparison (n={n}, k={k})")
    plt.xlabel("Epochs")
    plt.ylabel("Test Accuracy")
    plt.legend()
    plt.grid(True)
    
    plt.savefig(os.path.join(save_path, f"training_history_n{n}_k{k}.png"))
    plt.close()

def plot_k_scaling(results, n, k_range, save_path):
    """Plots final accuracy vs k for both tasks."""
    plt.figure(figsize=(10, 6))
    
    ks = [r['k'] for r in results]
    acc_exists = [r['final_test_acc_exists'] for r in results]
    acc_parity = [r['final_test_acc_parity'] for r in results]
    
    plt.plot(ks, acc_exists, marker='o', label="Detection")
    plt.plot(ks, acc_parity, marker='s', label="Parity")
    
    plt.axhline(y=0.5, color='r', linestyle='--', label="Random Guess (0.5)")
    
    plt.title(f"Performance Scaling with Clique Size k (n={n})")
    plt.xlabel("Clique Size (k)")
    plt.ylabel("Final Test Accuracy")
    plt.xticks(ks)
    plt.legend()
    plt.grid(True)
    
    plt.savefig(os.path.join(save_path, f"k_scaling_n{n}.png"))
    plt.close()
