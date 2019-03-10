import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import traceback
import seaborn as sns; sns.set()
sns.set_style("ticks")
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import sys
import pandas as pd

def plot_boxes(df):
    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=False, sharex=True, figsize=(10,10))
    titles="1% 2% 3% 3%/Q".split()
    categories="01 02 03 03_q".split()
    i = 0
    flat = axes.flatten()
    for axi in flat:
        cat = categories[i]
        df2 = df.loc[df['Category'] == cat]
        for tool in ["MASH", "VAPOR", "BLAST"]:
            df3 = df2.loc[df2['Tool'] == tool] 
            print(cat, tool, np.mean(df3['Score']))
            for perc in [75, 95, 99]:
                print("\t",np.percentile(df3['Score'], perc))
        ax = sns.boxplot(x="Tool", y='Score', data=df2, ax=axi)
        ax.set_title(titles[i])
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.yaxis.set_major_locator(ticker.MultipleLocator(max(df["Score"])/10))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        i += 1

    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07)
    plt.xlabel("Tool",labelpad=7, size=14)
    plt.ylabel("$L_A$",labelpad=7, size=14)
    plt.savefig("tool_scores_simulation.pdf", format="pdf", dpi=300)

headers = ['Category', 'Tool', 'Score']
table = []
fname = sys.argv[1]
with open(fname) as f:
    lines = [l.strip() for l in f]
    c = 0
    for line in lines[1:]:
        spl = line.split()
        cls = spl[0]
        original_distance = float(spl[6])
        mash_scores = [float(j) for j in spl[9].split(",")]
        vapor_scores = [float(j) for j in spl[12].split(",")]
        blast_scores = [float(j) for j in spl[15].split(",")]
        table.append([cls, "MASH", mash_scores[0]-original_distance])
        table.append([cls, "VAPOR", vapor_scores[0]-original_distance])
        table.append([cls, "BLAST", blast_scores[0]-original_distance])
#        if blast_scores[0] < vapor_scores[0]:
#            print(spl)
    print("Warning:", c, "of", len(lines[1:]), "datapoints were invalid. Please review the above lines.")
df = pd.DataFrame(table, columns=headers)

plot_boxes(df)
