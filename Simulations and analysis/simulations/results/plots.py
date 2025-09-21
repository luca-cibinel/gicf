import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import os

# %% HEADER

simulation = "time"
key = "p"
results_file = f"simulation_{simulation}.csv"
figures_folder = os.path.join("figures", simulation)

line_styles = {
        "MLE": "--",
        "LASSO": "--",
        "RIDGE": "-",
        "LRIDGE": "-"
    }

xaxis_start = {
        "MLE": 1,
        "LASSO": 1,
        "RIDGE": 0,
        "LRIDGE": 0
    }
    
# %% UTILITY

def csv_to_matrix(results, metric, model, method):
    rows = []
    
    simulations = np.unique(results.loc[:, "simulation"])
    
    results_loc = results.loc[results["metric"] == metric, :]
    results_loc = results_loc.loc[results_loc["n_bands"] == model, :]
    results_loc = results_loc.loc[results_loc["method"] == method, :]
    
    for s in simulations:
        rows += [results_loc.loc[results_loc["simulation"] == s, "value"]]
        
    return np.array(rows)

# %% MAIN
if not os.path.exists(figures_folder):
    os.makedirs(figures_folder)

print("Reading results...")
results = pd.read_csv(results_file, sep = " ")
print("Results retrieved!")

metrics = np.unique(results.loc[:, "metric"])
models = np.unique(results.loc[:, "n_bands"])
methods = np.unique(results.loc[:, "method"])
n = np.unique(results.loc[:, key])

print(f"Metrics: {metrics}")
print(f"Models: {models}")
print(f"Methods: {methods}")
print(f"{key}: {n}")

colors = ["black", "orange", "blue"]
markers = ["o", "v", "^", "P", "*", "X"]
require_01_ylims = ["F1", "d", "ePPV", "eTNR", "eTPR", "lambdar0"]
titles = {
        "F1": "$F_1$ score",
        "EL": "Entropy loss",
        "condnum": "Condition number",
        "cv.iters": "Average n. of iterations (CV)",
        "cv.time.elapsed": "Elapsed time (CV)",
        "cv.time.sys.self": "System time (CV)",
        "cv.time.user.self": "User time (CV)",
        "cv.all.iters": "Total n. of iterations (CV)",
        "cv.all.time.elapsed": "Total elapsed time (CV)",
        "cv.all.time.sys.self": "Total system time (CV)",
        "cv.all.time.user.self": "Total user time (CV)",
        "time.elapsed": "Elapsed time",
        "time.sys.self": "System time",
        "time.user.self": "User time",
        "d": "Estimated density",
        "iters": "Average number of iterations",
        "kappa": "$\kappa$",
        "lambda": "$\lambda$",
        "lambdar0": "$\lambda / \lambda_{MAX}(0)$"
    }

for metric in metrics:
    if metric.endswith(".child"):
        continue
    
    fig, ax = plt.subplots()
    
    for i, model in enumerate(models):
        
        for method in methods:
            
            if metric == "kappa" and method not in ["LRIDGE", "RIDGE"]:
                continue
            
            X = csv_to_matrix(results, metric, model, method)
            
            ax.plot(
                n[xaxis_start[method]:], 
                X.mean(0)[xaxis_start[method]:], 
                line_styles[method], 
                color = colors[i],
                linewidth = 1,
                label = f"{model} bands" if line_styles[method] == "-" else None
            )
            
            ax.scatter(
                n[xaxis_start[method]:], 
                X.mean(0)[xaxis_start[method]:], 
                marker = markers[i],
                s = 35,
                color = colors[i],
                facecolor = colors[i] if line_styles[method] == "-" else (0,0,0,0),
                linewidth = 0.5
            )
            
            ci = X.std(0)[xaxis_start[method]:]
            ax.fill_between(
                n[xaxis_start[method]:], 
                X.mean(0)[xaxis_start[method]:] - ci,  
                X.mean(0)[xaxis_start[method]:] + ci,
                color = colors[i],
                alpha = 0.1
            )
    
    if metric in require_01_ylims:
        ax.set_ylim(-0.05, 1.05)
    
    ax.set_xscale("log")
    ax.set_xticks(n)
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    
    ax.set_ylabel(titles.get(metric, metric))
    ax.set_xlabel(key)
    
    ax.legend()
    
    fig.savefig(os.path.join(figures_folder, f"{metric}.png"), format = "png")
    fig.savefig(os.path.join(figures_folder, f"{metric}.pdf"), format = "pdf")
            
            
    
    