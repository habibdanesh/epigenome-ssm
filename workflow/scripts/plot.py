import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import codecs
import json

### snakemake parameters
model_json = snakemake.input.model_json
model_err_file = snakemake.input.model_err_file
tracks = snakemake.params.tracks
n_features = snakemake.params.n_features
err_plot_file = snakemake.output.err_plot_file
emission_plot_file = snakemake.output.emission_plot_file
transition_plot_file = snakemake.output.transition_plot_file

### other parameters
fig_dpi = 600

### plot model error
err_k = np.fromfile(model_err_file)
err_plot = plt.figure().gca()
err_k = err_k[1:] # remove the initial error which might be too high
err_plot.plot(range(1, len(err_k)+1), err_k, label="K={}".format(n_features))
err_plot.xaxis.set_major_locator(MaxNLocator(integer=True))
err_plot.legend(loc="best")
plt.savefig(err_plot_file, dpi=fig_dpi)

### read model parameters
obj_text = codecs.open(model_json, 'r', encoding="utf-8").read()
model_param = json.loads(obj_text)
emission_mat = np.asmatrix(model_param["theta_m"])
transition_mat = np.asmatrix(model_param["lambda_m"])

f_ticks = []
for i in range(1, n_features+1):
    f_ticks.append("F{}".format(i))

### plot emission matrix
plt.clf()
ax = sns.heatmap(emission_mat.T, linewidth=0.5, cmap=sns.cm.rocket_r, 
                xticklabels=f_ticks,
                yticklabels=tracks)
plt.savefig(emission_plot_file, bbox_inches="tight", dpi=fig_dpi)

### plot transition matrix
plt.clf()
ax = sns.heatmap(transition_mat, linewidth=0.5, cmap=sns.cm.rocket_r, 
                xticklabels=f_ticks,
                yticklabels=f_ticks)
plt.savefig(transition_plot_file, bbox_inches="tight", dpi=fig_dpi)