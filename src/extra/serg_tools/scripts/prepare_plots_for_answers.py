import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from py.igquast_impl import save_plot, initialize_plot

# ------------

initialize_plot()

plt.title("umi tau choice: age3")

data = [45278, 40915, 40424, 40287]
plt.plot(range(len(data)), data, "b-", color="blue", label="label")

plt.xlabel("barcode tau")
plt.ylabel("# clusters")

save_plot("umi_tau_choice_age3", format="png")

# ------------

initialize_plot()

plt.title("umi tau choice: simulation (2 errors per read)")

data = [26489, 19203, 19055, 18923]
plt.plot(range(len(data)), data, "b-", color="blue", label="label")

plt.xlabel("barcode tau")
plt.ylabel("# clusters")

save_plot("umi_tau_choice_sim", format="png")

