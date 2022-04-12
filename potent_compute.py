import numpy as np
import scipy
from tqdm.notebook import tqdm
import pandas as pd

def compute_potentiation(df_spk_stims, window=(0, 100)):
    df_buffer = []
    df_header = ["Mouse", "Genotype", "Sex", "Date", "Session", "Mcc", "Channel", "Cell", "Grade", "Structure",
                 "Comparaison", "PSTH_x", "PSTH_diff", "Potent", "Potentiated", "Depressed"]
    df_spk_stims_copy = df_spk_stims.copy()
    df_spk_stims_copy = df_spk_stims_copy[df_spk_stims_copy["Session"]=="Theta1"]
    mcc_list = df_spk_stims_copy["Mcc"].unique()
    nb_cells = len(mcc_list)
    for m in tqdm(range(nb_cells), desc="Iterating over cells"):
        mcc = mcc_list[m]
        df_mcc = df_spk_stims[df_spk_stims["Mcc"] == mcc]
        mouse = str(list(df_mcc["Mouse"])[0])
        date = str(list(df_mcc["Date"])[0])
        gen = str(list(df_mcc["Genotype"])[0])
        sex = str(list(df_mcc["Sex"])[0])
        ch = int(list(df_mcc["Channel"])[0])
        c = int(list(df_mcc["Cell"])[0])
        grade = int(list(df_mcc["Grade"])[0])
        session = str(list(df_mcc["Session"])[0])
        struct = str(list(df_mcc["Structure"])[0])
        line_base = [mouse, gen, sex, date, session, mcc, ch, c, grade, struct]
        x = list(df_mcc["PSTH_x"])[0]
        comparaisons = ["preOxo_preSaline", "Saline_postpre", "Oxo_postpre", "pre_OxoSaline", "post_OxoSaline",
                        "preOxo_postSaline"]
        comp_dict = {"preOxo_preSaline": [2, 0], "Saline_postpre": [1, 0], "Oxo_postpre": [3, 2],
                     "pre_OxoSaline": [2, 0], "post_OxoSaline": [3, 1], "preOxo_postSaline": [2, 1]}

        psth_mcc_sal_before = \
        list(df_mcc[(df_mcc["Prepost"] == "before") & (df_mcc["Treatment"] == "Saline")]["PSTH_centered"])[0]
        psth_mcc_sal_post = \
        list(df_mcc[(df_mcc["Prepost"] == "post") & (df_mcc["Treatment"] == "Saline")]["PSTH_centered"])[0]
        psth_mcc_oxo_before = \
        list(df_mcc[(df_mcc["Prepost"] == "before") & (df_mcc["Treatment"] == "Oxo")]["PSTH_centered"])[0]
        psth_mcc_oxo_post = \
        list(df_mcc[(df_mcc["Prepost"] == "post") & (df_mcc["Treatment"] == "Oxo")]["PSTH_centered"])[0]
        list_psth = [psth_mcc_sal_before, psth_mcc_sal_post, psth_mcc_oxo_before, psth_mcc_oxo_post]
        for i in range(len(comparaisons)):
            comp = comp_dict[comparaisons[i]]
            comp_label = comparaisons[i]
            bins, psth_d, potent = get_potent(x, x, list_psth[comp[1]], list_psth[comp[0]], window=window)
            if potent > 0:
                dep = 0
                pot = 1
            else:
                dep = 1
                pot = 0
            df_buffer.append(line_base + [comp_label, bins, psth_d, potent, pot, dep])
    df_diff = pd.DataFrame(df_buffer, columns=df_header)
    return df_diff



def get_potent(x1, x2, psth1, psth2, window=(0, 100)):
    mini1 = np.searchsorted(x1, window[0], side="left")
    mini2 = np.searchsorted(x2, window[0], side="left")
    maxi1 = np.searchsorted(x1, window[1], side="right")
    maxi2 = np.searchsorted(x2, window[1], side="right")
    diff1 = maxi1-mini1
    diff2 = maxi2-mini2
    if diff1 != diff2:
        maxi = min([mini1+diff1, mini2+diff2])
    else:
        maxi = maxi1
    psth1_crop = np.array(psth1[mini1:maxi])
    psth2_crop = np.array(psth2[mini2:maxi])
    diff_psth = psth2_crop - psth1_crop
    potent = np.nanmean(diff_psth)
    return x1[mini1:maxi], diff_psth, potent