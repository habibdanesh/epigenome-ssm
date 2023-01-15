import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import codecs
import json

### snakemake parameters
model_json = snakemake.input.model_json
tracks = snakemake.params.tracks
n_features = snakemake.params.n_features
err_plot_file = snakemake.output.err_plot_file
time_plot_file = snakemake.output.time_plot_file
emission_plot_file = snakemake.output.emission_plot_file
transition_plot_file = snakemake.output.transition_plot_file

### other parameters
fig_dpi = 600

### read model parameters
obj_text = codecs.open(model_json, 'r', encoding="utf-8").read()
model_params = json.loads(obj_text)
model_type = model_params["model_type"]
emission_mat = np.asmatrix(model_params["theta_m"])
transition_mat = np.asmatrix(model_params["lambda_m"])
error_list = model_params["error_m"]
time_list = model_params["opt_time_m"]
epigenomes = model_params["epigenomes"]
assays = model_params["assays"]
model_params = None

### plot model error
err_plot = plt.figure().gca()
error_list = error_list[1:] # remove the initial error which might be too high
err_plot.plot(range(1, len(error_list)+1), error_list, label="K={}".format(n_features))
err_plot.xaxis.set_major_locator(MaxNLocator(integer=True))
err_plot.legend(loc="best")
plt.savefig(err_plot_file, dpi=fig_dpi)

### plot optimization time
time_plot = plt.figure().gca()
time_list = [t/60 for t in time_list] # convert to minutes
time_list = [t/60 for t in time_list] # convert to hours
time_plot.plot(time_list, error_list, label="K={}".format(n_features))
time_plot.xaxis.set_major_locator(MaxNLocator(integer=False))
time_plot.legend(loc="best")
plt.savefig(time_plot_file, dpi=fig_dpi)

f_ticks = []
for i in range(1, n_features+1):
    f_ticks.append("F{}".format(i))

### plot emission matrix
plt.clf()
if model_type == "stacked":
    y_tick_labels = tracks
if model_type == "concatenated":
    y_tick_labels = assays
ax = sns.heatmap(emission_mat.T, linewidth=0.5, cmap=sns.cm.rocket_r, 
                xticklabels=f_ticks,
                yticklabels=y_tick_labels)
plt.savefig(emission_plot_file, bbox_inches="tight", dpi=fig_dpi)

### plot transition matrix
plt.clf()
ax = sns.heatmap(transition_mat, linewidth=0.5, cmap=sns.cm.rocket_r, 
                xticklabels=f_ticks,
                yticklabels=f_ticks)
plt.savefig(transition_plot_file, bbox_inches="tight", dpi=fig_dpi)