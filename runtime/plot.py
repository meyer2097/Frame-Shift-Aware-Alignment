import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

medianDF = pd.read_csv("median.csv", index_col=0, header=0)
meanDF = pd.read_csv("mean.csv", index_col=0, header=0)

df = medianDF

df = df.applymap(lambda x: x*1000)
df = df[::-1]

ax = sns.heatmap(df, cbar_kws={'label': 'Time in ms'})
ax.set(title="Median time to align two random sequences\nSingle Intel i5-8250U @ 3.400GHz core")
ax.set(ylabel="Number ob Nucleotides")
ax.set(xlabel="Number of Aminoacids")

plt.show()
