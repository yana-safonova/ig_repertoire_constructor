import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from py.igquast_impl import save_plot, initialize_plot

min_tau = 10
max_tau = 50
step = 10

tau_missing = []
missing = []
tau_extra = []
extra = []

for tau in range(min_tau, max_tau + step, step):
    # tau_missing = []
    # missing = []
    tau_extra = []
    extra = []
    cluster_size = []
    total_clusters = 0
    zero = 0
    # for line in open("dist_%d/missing" % tau):
    #     tau_missing.append(tau)
    #     missing.append(float(line.strip()))
    for line in open("dist_%d/extra" % tau):
        total_clusters += 1
        pair = line.strip().split(" ")
        value = float(pair[0])
        size = int(pair[1])
        if value == 0:
            zero += 1
        else:
            tau_extra.append(tau)
            extra.append(value)
            cluster_size.append(size)
    # plot = sns.distplot(missing, bins=20, kde=False, ax=ax)
    # plot.set_xlabel("extra reads")
    # plot.set_ylabel("# collisions")
    # plt.savefig("missing_%d.png" % tau, format="png")
    initialize_plot()
    print max(extra)
    plot = sns.distplot(np.log(np.array(extra) + 1), bins=20, kde=False)
    plt.xlim(0, int(max(np.log(np.array(extra) + 1))) + 2)
    plot.set_xlabel("extra reads")
    # plot.set_ylabel("# collisions")
    plt.title("tau %d: # collisions" % tau)
    plt.savefig("extra_%d.png" % tau, format="png")
    print float(zero) / total_clusters

    initialize_plot()
    plot = sns.jointplot(np.log(np.array(extra) + 1), np.log(np.array(cluster_size)), stat_func=None)
    plot.set_axis_labels("constructed abundance / reference abundance", "constructed abundance")
    plt.savefig("extra_joint_%d.png" % tau, format="png")


# plot = sns.regplot(np.array(tau_missing), np.array(missing), fit_reg = False)
# plot.set_xlabel("clustering tau")
# plot.set_ylabel("missing reads")
# plt.savefig("missing.png", format="png")
#
# # sns.heatmap(zip(tau_extra, extra))
# plot = sns.regplot(np.array(tau_extra), np.log(np.array(extra) + 1), fit_reg = False)
# plot.set_xlabel("clustering tau")
# plot.set_ylabel("extra reads")
# plt.savefig("extra.png", format="png")
