import numpy as np
from utilities.peri_stimulus import PSTH

THRESH = 3
VOLLEY = 10

# Function to compute the PSTH over the 100ms of optogenetic stimulations
def compute_psth(df_spk_stims, df_stims, binsize=10, window=(-300, 200), stim_start=0, stim_end=100, reject_ranges=None):
    df_spk_stims_copy = df_spk_stims.copy()
    psth_list = []
    psth_raw = []
    psth_centered = []
    x_lags = []
    z_score = []
    dfr_stim = []
    for index, row in df_spk_stims.iterrows():
        mouse = str(row["Mouse"])
        date = int(row["Date"])
        bl = int(row["Block"])
        stims = np.array(df_stims[(df_stims["Mouse"]==mouse)&(df_stims["Date"]==date)&(df_stims["Block"]==bl)]["Timestamp"])
        psth = PSTH(row["Spiketrain"], stims, binsize=binsize, window=window, reject_ranges=reject_ranges)
        psth_list.append(psth)
        psth_raw.append(psth.get_rate())
        psth_centered.append(psth.get_centered())
        z_score.append(psth.get_zscore())
        x_lags.append(psth.get_timestamps())
        rate = psth.get_rate()
        x = psth.get_timestamps()
        mini = np.searchsorted(x, stim_start, side="left")
        maxi = np.searchsorted(x, stim_end, side="left")
        dfr_stim.append(np.nanmean(rate[mini:maxi])-np.nanmean(rate[:mini]))
    df_spk_stims_copy["PSTH"] = psth_list
    df_spk_stims_copy["PSTH_raw"] = psth_raw
    df_spk_stims_copy["PSTH_centered"] = psth_centered
    df_spk_stims_copy["PSTH_z"] = z_score
    df_spk_stims_copy["PSTH_x"] = x_lags
    df_spk_stims_copy["PSTH_deltastim"] = dfr_stim
    return df_spk_stims_copy


def compute_psth_afferent_volley(df_spk_stims, df_stims, binsize=1, window=(-300, 200), volley=VOLLEY, reject_ranges=None):
    df_spk_stims_copy = df_spk_stims.copy()
    psth_list = []
    psth_raw = []
    psth_centered = []
    x_lags = []
    z_score = []
    dfr_stim = []
    thal = 4
    cort_str = 7
    stim_starts = {"CL":thal, "VAL":thal, "M1":cort_str, "DLS":cort_str, "DN":0}
    for index, row in df_spk_stims.iterrows():
        mouse = str(row["Mouse"])
        date = int(row["Date"])
        bl = int(row["Block"])
        structure = str(row["Structure"])
        stim_start = stim_starts[structure]
        stims = np.array(df_stims[(df_stims["Mouse"]==mouse)&(df_stims["Date"]==date)&(df_stims["Block"]==bl)]["Timestamp"])
        stims =  stims + float(stim_start)
        psth = PSTH(row["Spiketrain"], stims, binsize=binsize, window=window, reject_ranges=reject_ranges)
        psth_list.append(psth)
        psth_raw.append(psth.get_rate())
        psth_centered.append(psth.get_centered())
        z_score.append(psth.get_zscore())
        x_lags.append(psth.get_timestamps())
        rate = psth.get_rate()
        x = psth.get_timestamps()
        mini = np.searchsorted(x, 0, side="left")
        maxi = np.searchsorted(x, volley, side="left")
        dfr_stim.append(np.nanmean(rate[mini:maxi])-np.nanmean(rate[:mini]))
    df_spk_stims_copy["PSTH"] = psth_list
    df_spk_stims_copy["PSTH_raw"] = psth_raw
    df_spk_stims_copy["PSTH_centered"] = psth_centered
    df_spk_stims_copy["PSTH_z"] = z_score
    df_spk_stims_copy["PSTH_x"] = x_lags
    df_spk_stims_copy["PSTH_deltastim"] = dfr_stim
    return df_spk_stims_copy


# Function to get the list of identifiers of responsive cells
def get_mcc_stimresp(df, thr=THRESH):
    mcc_buffer = []
    x = list(df["PSTH_x"])[0]
    mini = np.searchsorted(x, 0, side="left")
    maxi = np.searchsorted(x, 100, side="right")
    for index, row in df.iterrows():
        if np.nanmax(list(row["PSTH_z"])[mini:maxi]) >= thr:
            print(np.nanmax(list(row["PSTH_z"])[mini:maxi]))
            mcc_buffer.append(str(row["Mcc"]))
    return mcc_buffer
