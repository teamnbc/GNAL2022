import numpy as np
import scipy
import scipy.stats
import itertools

class PSTH():
    """
        Computes PSTH for a given spiketrain and a list of events
        Arguments:
            - spiketrain: list of spikes [ms] (array-like)
            - events: list of events [ms] (array-like)
            - binsize: width of histogram bins [ms] (float)
            - window: start and end of window (start, end) [ms] (array-like)
            - scale_factor: scale related to time unit [1000 for ms] (float)
            - timestamp_select: bin edge selection (left, right or center) (string)
            - tokens: (dict)
        Attributes:
            - bins: histogram bins [ms] (array-like)
            - timestamps: time associated to each bin (array-like)
            - spikecount: spikecount timeseries (array-like, nbins)
        Functions:
            - get_timestamps(): returns timestamps
            - get_spikecount(): returns spikecount timeseries
            - get_rate(): returns rate timeseries
            - get_norm(method=np.average): returns normalized timeseries (method applied on baseline)
            - get_zscore(method="classical"): returns zscore timeseries (method: classical X-mean/SD, robust X-med/MAD)
    """
    def __init__(self, spiketrain, events, binsize, window=(-500, 500), scale_factor=1000, timestamp_select="left", reject_ranges=None, tokens={}):
        self.events = np.array(events)
        self.window = window
        self.scale_factor = scale_factor
        self.timestamp_select = timestamp_select
        self.tokens = tokens
        self.binsize = binsize
        self.reject_ranges = reject_ranges
        self.bins = np.arange(window[0], window[1]+binsize, binsize)
        self.timestamps = self.compute_timestamps(self.bins, self.timestamp_select)
        self.spikecount = self.compute_histogram(spiketrain, events, window, self.bins)


    def get_timestamps(self):
        return self.timestamps

    def get_spikecount(self):
        return self.spikecount

    def get_rate(self):
        return self.spikecount * self.scale_factor / (self.binsize * len(self.events))

    def get_norm(self, method=np.average):
        maxi = np.searchsorted(self.timestamps, 0, side="left")
        baseline = self.spikecount[:maxi]
        return self.spikecount / method(baseline)

    def get_centered(self, method=np.average, rate_norm=True):
        maxi = np.searchsorted(self.timestamps, 0, side="left")
        baseline = self.spikecount[:maxi]
        centered = self.spikecount - method(baseline)
        if rate_norm:
            centered =  centered * self.scale_factor / (self.binsize * len(self.events))
        return centered

    def get_zscore(self, method="classical"):
        maxi = np.searchsorted(self.timestamps, 0, side="left")
        new_mat = np.zeros_like(self.spikecount)
        baseline = self.spikecount[:maxi]
        if method == "classical":
            new_mat = (self.spikecount - np.nanmean(baseline)) / np.nanstd(baseline)
        elif method == "robust":
            new_mat = (self.spikecount - np.nanmedian(baseline)) / scipy.stats.median_abs_deviation(baseline, nan_policy="omit")
        return new_mat

    def compute_histogram(self, spiketrain, events, window, bins):
        crops = list(map(lambda x: self.crop_window(spiketrain, x, window), events))
        merged = list(itertools.chain.from_iterable(crops))
        hist, bin_edges = np.histogram(merged, bins)
        return hist

    def crop_window(self, spiketrain, event, window):
        mini = np.searchsorted(spiketrain, event + window[0], side="left")
        maxi = np.searchsorted(spiketrain, event + window[1], side="right")
        spk = np.array(spiketrain[mini:maxi]) - event
        if self.reject_ranges is not None:
            if len(np.array(self.reject_ranges).shape) == 1: # Only one range to reject
                ranges = [self.reject_ranges]
            else:
                ranges = self.reject_ranges
            for range in ranges:
                mini = np.searchsorted(spk, range[0], side="left")
                maxi = np.searchsorted(spk, range[1], side="right")
                spk2 = np.array(spk[mini:maxi])
                spk = np.setdiff1d(spk, spk2)
        return spk

    def compute_timestamps(self, bins, method):
        if method == "left":
            return bins[:-1]
        elif method == "right":
            return bins[1:]
        elif method == "center":
            diff = np.diff(bins)
            return np.array(bins) + diff