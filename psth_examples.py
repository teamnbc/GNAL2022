import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

PATH_FIG = ""


# Function to display the PSTH and rasterplot of the cell "mcc"
def plot_psth(df, df_stims, mcc, session= "Theta1", tt="Saline", moment="before", name=None, ylim=None):
    df = df[(df["Mcc"] == mcc) & (df["Treatment"] == tt) & (df["Session"] == session) & (df["Prepost"] == moment)]
    x = list(df["PSTH_x"])[0]
    y = np.array(list(df["PSTH_raw"])[0])
    spk = list(list(df["Spiketrain"])[0])
    mouse = list(df["Mouse"])[0]
    date = list(df["Date"])[0]
    bl = list(df["Block"])[0]
    stims = np.array(
        df_stims[(df_stims["Mouse"] == mouse) & (df_stims["Date"] == int(date)) & (df_stims["Block"] == bl)][
            "Timestamp"])
    func = interpolate.interp1d(x, y, kind="previous")
    buffer = []
    for s in stims:
        mini = np.searchsorted(spk, s - 300)
        maxi = np.searchsorted(spk, s + 200)
        buffer.append(list(np.array(spk[mini:maxi]) - s))

    x_cropped = np.linspace(0, 100, 1000)
    y_cropped = func(x_cropped)

    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)

    ax1.eventplot(buffer, color="black")
    ax1.set_ylim([len(buffer) + 1, -1])
    ax2.step(x, y, color="black", where="post")
    ax2.fill_between(x_cropped, y_cropped, 0, color="blue", alpha=0.5, label="Cb stim")
    ax2.axvline(0, color="blue")
    ax2.axvline(100, color="blue")

    ax1.axvline(0, color="blue")
    ax1.axvline(100, color="blue")
    ax1.set_ylabel("Trial")
    ax2.set_ylabel("Firing Rate (Hz)")
    ax2.set_xlabel("Time (ms)")
    ax1.set_title(tt)
    ax1.set_xlim([-300, 150])
    if ylim is not None:
        ax2.set_ylim(ylim)
    plt.tight_layout()
    if name is not None:
        plt.savefig(PATH_FIG + "psth_{n}_{m}.pdf".format(n=name, m=mcc), format="pdf")
    plt.show()

# Function to display the overlay of PSTH before and after Theta-burst stimulation of the cell "mcc"
def plot_potentiation(df, mcc, tt="Saline", name=None):
    df_bef = df[(df["Mcc"] == mcc)&(df["Treatment"] == tt)&(df["Session"] == "Theta1")&(df["Prepost"] == "before")]
    df_af = df[(df["Mcc"]==mcc) & (df["Treatment"] == tt)&(df["Session"] == "Theta1")&(df["Prepost"] == "post")]
    x = list(df_bef["PSTH_x"])[0]
    y_bef = np.array(list(df_bef["PSTH_centered"])[0])
    y_af = np.array(list(df_af["PSTH_centered"])[0])
    int_bef = interpolate.interp1d(x, y_bef, kind="previous")
    int_af = interpolate.interp1d(x, y_af, kind="previous")

    x_cropped = np.linspace(0, 100, 1000)
    y_bef_cropped = int_bef(x_cropped)
    y_af_cropped = int_af(x_cropped)
    pot = np.where(y_af_cropped > y_bef_cropped, y_af_cropped - y_bef_cropped, 0)
    dep = np.where(y_af_cropped < y_bef_cropped, y_af_cropped - y_bef_cropped, 0)

    f, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})
    ax1.step(x, y_bef, color="black", label="pre-TB", where="post")
    ax1.step(x, y_af, color="green", label="post-TB", where="post")
    ax1.fill_between(x_cropped, y_bef_cropped, y_bef_cropped + pot, color="green", alpha=0.5, label="potentiation")
    ax1.fill_between(x_cropped, y_bef_cropped, y_bef_cropped + dep, color="black", alpha=0.5, label="depression")
    ax1.axhline(0, color="black", linestyle=":")
    ax1.axvline(0, color="blue")
    ax1.axvline(100, color="blue")

    ax2.step(x, y_bef, color="black", label="pre-TB", where="post")
    ax2.step(x, y_af, color="green", label="post-TB", where="post")
    ax2.fill_between(x_cropped, y_bef_cropped, y_bef_cropped + pot, color="green", alpha=0.5, label="potentiation")
    ax2.fill_between(x_cropped, y_bef_cropped, y_bef_cropped + dep, color="black", alpha=0.5, label="depression")
    ax2.axhline(0, color="black", linestyle=":")
    ax2.axvline(0, color="blue")
    ax2.axvline(4, color="black", linestyle=":")
    ax2.axvline(14, color="black", linestyle=":")

    ax1.set_ylabel("Delta FR from baseline (Hz)")
    ax1.set_xlabel("Time (ms)")
    ax1.set_title(tt)
    ax1.set_xlim([-300, 150])
    ax2.set_xlim([-50, 30])
    ax1.legend()
    if name is not None:
        plt.savefig(PATH_FIG + "{n}_{m}.pdf".format(n=name, m=mcc), format="pdf")
    plt.show()