import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import os

# %% HEADER

simulation = "mle"
key = "p" if simulation == "time" else "n"
fig_suffix = "_p_50" if simulation != "time" else ""
results_file = f"simulation_{simulation}.csv"
figures_folder = os.path.join("figures", simulation)
legends_folder = os.path.join("figures", "legends")

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

method_name = {
        "MLE": "MLE",
        "LASSO": "covglasso",
        "RIDGE": "GICF",
        "LRIDGE": "GICF"
    }

colors = ["#000000", "#ffa500", "#0000ff"]
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
    
if not os.path.exists(legends_folder):
    os.makedirs(legends_folder)

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

legend_drawn = False

for metric in metrics:
    if metric.endswith(".child"):
        continue
    
    fig, ax = plt.subplots()
    
    for i, model in enumerate(models):
        
        for method in methods:
            
            if metric == "kappa" and method not in ["LRIDGE", "RIDGE"]:
                continue
            
            X = csv_to_matrix(results, metric, model, method)
            
            xstart = xaxis_start[method] if key == "n" else 0
            
            ax.plot(
                n[xstart:], 
                X.mean(0)[xstart:], 
                line_styles[method], 
                color = colors[i],
                linewidth = 1,
                marker = markers[i],
                mfc = colors[i] + ("" if line_styles[method] == "-" else "00"),
                mec = colors[i],
                #ms = 35,
                label = f"{model} band{'s' if model > 1 else ''}"
            )
            
            ci = X.std(0)[xstart:]
            ax.fill_between(
                n[xstart:], 
                X.mean(0)[xstart:] - ci,  
                X.mean(0)[xstart:] + ci,
                color = colors[i],
                alpha = 0.1
            )
            
            if metric == "condnum":
                ax.axhline(y=50, color='gray', linestyle='-')
    
    if metric in require_01_ylims:
        ax.set_ylim(-0.05, 1.05)
    
    ax.set_xscale("log")
    ax.set_xticks(n)
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    
    ax.set_ylabel(titles.get(metric, metric))
    ax.set_xlabel(key)

    if not legend_drawn:
        h, l = ax.get_legend_handles_labels()
        handles = [plt.plot([],marker="", ls="")[0]]*2 + h
        labels = [f"{method_name[m]}:" for m in methods] + l
        
        _, ax_temp = plt.subplots()
        ax_temp.axis("off")
        legend = ax_temp.legend(handles, labels, frameon = False, ncols = len(models) + 1, facecolor = "white", framealpha = 0)
        
        fig_l  = legend.figure
        fig_l.canvas.draw()
        bbox_l  = legend.get_window_extent().transformed(fig_l.dpi_scale_trans.inverted())
        
        fig_l.savefig(os.path.join(legends_folder, f"legend_{simulation}_horiz.png"), bbox_inches = bbox_l, format = "png")
        fig_l.savefig(os.path.join(legends_folder, f"legend_{simulation}_horiz.pdf"), bbox_inches = bbox_l, format = "pdf")
    
        legend.remove()
        
        legend_drawn = True
    
    fig.savefig(os.path.join(figures_folder, f"{simulation}_{metric}{fig_suffix}.png"), format = "png")
    fig.savefig(os.path.join(figures_folder, f"{simulation}_{metric}{fig_suffix}.pdf"), format = "pdf")
            
            
    
    