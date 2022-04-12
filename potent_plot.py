import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks")

PATH_FIG = ""


# Function to generates the figures associated to potentiation
def plot_tbpotent(df_potent, structures=("VAL", "CL", "M1"), tt="Saline", name=None):
    if tt == "Saline":
        comp = "Saline_postpre"
    else:
        comp = "Oxo_postpre"
    gen = ["WT", "HET"]
    cols = {"WT": "gray", "HET": "orange"}

    df_potent = df_potent[(df_potent["Comparaison"] == comp) & (df_potent["Structure"].isin(structures))]
    fig = plt.figure(figsize=(10, 5))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    sns.boxplot(data=df_potent, x="Structure", y="Potent", hue="Genotype", order=structures, hue_order=gen, ax=ax1,
                palette=cols)
    sns.swarmplot(data=df_potent, x="Structure", y="Potent", hue="Genotype", order=structures, hue_order=gen, ax=ax1,
                  dodge=True, palette=cols, edgecolor="black", linewidth=1)

    sns.barplot(data=df_potent, x="Structure", y="Potent", hue="Genotype", order=structures, hue_order=gen, ax=ax2,
                palette=cols)
    ax1.axhline(0, linestyle=":", color="black")
    ax2.axhline(0, linestyle=":", color="black")
    ax1.set_ylabel("DeltaFR induced by TB (Hz)")
    plt.tight_layout()
    if name is not None:
        plt.savefig(PATH_FIG + name + ".pdf", format="pdf")
    plt.show()