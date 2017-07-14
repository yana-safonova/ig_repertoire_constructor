import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from py.igquast_impl import save_plot, initialize_plot

# ------------

initialize_plot()

# plt.title("umi tau choice: age3")

data = [45278, 40915, 40424, 40287]
plt.plot(range(len(data)), data, "b-", color="blue", label="label")

plt.xlabel("barcode tau")
plt.ylabel("# clusters")

save_plot("umi_tau_choice_age3", format="png")

# ------------

initialize_plot()

# plt.title("umi tau choice: simulation (2 errors per read)")

data = [26489, 19203, 19055, 18923]
plt.plot(range(len(data)), data, "b-", color="blue", label="label")

plt.xlabel("barcode tau")
plt.ylabel("# clusters")

save_plot("umi_tau_choice_sim", format="png")

# ------------

initialize_plot()

data = [39.13532894, 67.98549154, 68.53164346, 68.10431085]
plt.plot(range(len(data)), data, "b-", color="blue", label="label")
# plt.ylim(0, 1)

plt.xlabel("barcode tau")
plt.ylabel("% reads with corrected UMIs")

save_plot("umi_tau_choice_sim_reads", format="png")

# ------------

initialize_plot()

# data = [152, 1056, 3731, 9474, 14092, 15561, 14491, 12075, 9233, 6642]
data = [0.3476669716, 2.453075636, 9.200305773, 26.9960677, 45.92322232, 53.09109519, 47.78250404, 37.0330614, 26.21670737, 17.66160555]
tau = range(0, len(data) * 5, 5)
plt.plot(tau, data, "b-", color="blue", label="label")
# plt.ylim(0, 1)

plt.xlabel("chimera tau")
plt.ylabel("% detected chimeras")

save_plot("chimera_tau_choice", format="png")

