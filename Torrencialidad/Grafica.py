#%%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Database de cuencas torrenciales y no torrenciales
df = pd.read_csv('./torrencialidad.csv')
df_cuencas = pd.read_csv('./morfometria.csv', sep=';')

# Set the seaborn style and context
sns.set(context="poster", style="whitegrid")

# Create subplots with 1 row and 3 columns
fig, axes = plt.subplots(1, 3, figsize=(20, 7), sharey=True)

# Plot 1 - A vs Rh
sns.scatterplot(x="A (sq. km)", y="Rh", hue="FFR", style="FFR",
                markers={"nT": "o", "T": "X"},
                data=df[df['FFR'].isin(['nT', 'T'])],
                s=60, palette=["red", "blue"], alpha=0.8, ax=axes[0], legend=False)
sns.scatterplot(x="Area", y="Rh", data=df_cuencas,
                s=100, color='darkgreen', alpha=0.8, ax=axes[0], markers='.')
axes[0].grid(alpha=.4)

# Plot 2 - Su vs Rh
sns.scatterplot(x="Su", y="Rh", hue="FFR", style="FFR",
                markers={"nT": "o", "T": "X"},
                data=df[df['FFR'].isin(['nT', 'T'])],
                s=60, palette=["red", "blue"], alpha=0.8, ax=axes[1], legend=False)
sns.scatterplot(x="Su", y="Rh", data=df_cuencas,
                s=100, color='darkgreen', alpha=0.8, ax=axes[1], markers='.')
axes[1].grid(alpha=.4)

# Plot 3 - M vs Rh
sns.scatterplot(x="M", y="Rh", hue="FFR", style="FFR",
                markers={"nT": "o", "T": "X"},
                data=df[df['FFR'].isin(['nT', 'T'])],
                s=60, palette=["red", "blue"], alpha=0.8, ax=axes[2], legend=False)
sns.scatterplot(x="M", y="Rh", data=df_cuencas,
                s=100, color='darkgreen', alpha=0.8, ax=axes[2], markers='.')

# Create empty scatter plot artists for legend items
legend_markers = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=15, alpha=0.8),
    Line2D([0], [0], marker='X', color='w', markerfacecolor='red', markersize=15, alpha=0.8),
    Line2D([0], [0], marker='.', color='w', markerfacecolor='green', markersize=25, alpha=0.8)
]

# Add custom legend to the third subplot
axes[2].legend(legend_markers, ['No torrencial', 'Torrencial', 'Cuencas estudio'], fontsize=15)

axes[2].grid(alpha=.4)
plt.tight_layout()
plt.show()


#%%
# Gráfica Melton
plt.figure(figsize=(10,8))
sns.scatterplot(x="M", y="Lb",
                s=80, markers='.',
                data=df_cuencas, color='darkgreen', alpha=0.8)
# custom = [Line2D([], [], marker='^', color='Lime', linestyle='None', alpha=0.8)]
# plt.arrow(0.7, 20, 0.02, 17.3)
x1=(0.3,0.3);y1=(0,90);plt.plot(x1,y1,'k', zorder=0)
x2=(0.6,0.6);y2=(0,2.7);plt.plot(x2,y2,'k')
x3=(0.6,0.8);y3=(2.7,2.7);plt.plot(x3,y3,'k', zorder=0)
plt.axis([0, 0.8, 0, 90])
plt.annotate('Floods',(0.1,55))
plt.annotate('Debris floods',(0.4,40))
plt.annotate('Debris flows',(0.6,10))

plt.tight_layout()
# %%
