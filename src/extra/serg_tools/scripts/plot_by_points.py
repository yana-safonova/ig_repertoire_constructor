import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from py.igquast_impl import save_plot, initialize_plot

initialize_plot()
# plt.ioff()
# sns.set_style("darkgrid")
# sns.set(font_scale=font_scale)

data = [45278, 40915, 40424, 40287]
plt.plot(range(len(data)), data, "b-", color="blue", label="label")

# plt.xlim((0 - eps, 1 + eps))
# plt.ylim((0 - eps, 1 + eps))

plt.xlabel("barcode tau")
plt.ylabel("# clusters")

save_plot("plot", format="png")

