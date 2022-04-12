import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks")

PATH_FIG = ""
COLORS_SESSION = {"Naive": "C0", "Oxo": "red", "Exposed": 'C2'}


# Function to generate the figures associated to the description of the different pathological states
def endophenotype_plot(df, param="PSTH_deltastim", name=None):
    sessions = ["Naive", "Oxo", "Exposed"]
    df = df[df["Prepost"] != "post"]
    df.loc[(df["Session"] == 'BslOxo') & (df["Treatment"] == 'Saline'), 'Session'] = 'Naive'
    df.loc[(df["Session"] == 'BslOxo') & (df["Treatment"] == 'Oxo'), 'Session'] = 'Oxo'
    df.loc[(df["Session"] == 'Theta1') & (df["Treatment"] == 'Saline'), 'Session'] = 'Exposed'

    df_c = df[(df["Session"].isin(sessions)) & (df["Structure"].isin(["VAL", "CL", "M1"]))]
    fig = plt.figure(figsize=(12, 5))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132, sharey=ax1)
    ax3 = plt.subplot(133, sharey=ax1)
    sns.boxplot(data=df_c[df_c["Structure"] == "VAL"], x="Genotype", y=param, hue="Session", order=["WT", "HET"],
                palette=COLORS_SESSION, hue_order=sessions, ax=ax1)
    sns.swarmplot(data=df_c[df_c["Structure"] == "VAL"], x="Genotype", y=param, hue="Session", order=["WT", "HET"],
                  palette=COLORS_SESSION, hue_order=sessions, ax=ax1, dodge=True, edgecolor="black", linewidth=1)

    sns.boxplot(data=df_c[df_c["Structure"] == "CL"], x="Genotype", y=param, hue="Session", order=["WT", "HET"],
                palette=COLORS_SESSION, hue_order=sessions, ax=ax2)
    sns.swarmplot(data=df_c[df_c["Structure"] == "CL"], x="Genotype", y=param, hue="Session", order=["WT", "HET"],
                  palette=COLORS_SESSION, hue_order=sessions, ax=ax2, dodge=True, edgecolor="black", linewidth=1)

    sns.boxplot(data=df_c[df_c["Structure"] == "M1"], x="Genotype", y=param, hue="Session", order=["WT", "HET"],
                palette=COLORS_SESSION, hue_order=sessions, ax=ax3)
    sns.swarmplot(data=df_c[df_c["Structure"] == "M1"], x="Genotype", y=param, hue="Session", order=["WT", "HET"],
                  palette=COLORS_SESSION, hue_order=sessions, ax=ax3, dodge=True, edgecolor="black", linewidth=1)
    if param == "MeanFR":
        ax1.set_ylabel("Firing Rate (Hz)")
        ax2.set_ylabel(None)
        ax3.set_ylabel(None)
    ax1.set_title("VAL")
    ax2.set_title("CL")
    ax3.set_title("M1")
    ax1.axhline(0, linestyle=":", color="black")
    ax2.axhline(0, linestyle=":", color="black")
    ax3.axhline(0, linestyle=":", color="black")
    if name is not None:
        plt.savefig(PATH_FIG + name + ".pdf", format="pdf")
    plt.show()