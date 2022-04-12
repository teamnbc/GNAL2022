import numpy as np

# Computes basic spiketrain statistics
class Spiketrain():
    def __init__(self, spikes, start=None, stop=None, tokens={}):
        if type(spikes[0]) == list or type(spikes[0]) == np.array:
            spikes = list(map(self.sort, spikes))
            spikes = list(map(lambda x: self.cleanSpk(x, start, stop), spikes))
        else:
            spikes.sort()
            spikes = self.cleanSpk(spikes, start, stop)
        self.spk = spikes
        self.start = start
        self.stop = stop
        self.tokens = tokens

    def meanFR(self, x=None, start=None, stop=None, factor=1000):
        if x is None:
            x = self.spk
        if start is None:
            start = self.start
        if stop is None:
            stop = self.stop
        return len(x) * factor / (stop - start)

    def isi(self, x=None):
        if x is None:
            x = self.spk
        return np.diff(x)

    def cv(self, x=None):
        isi = self.isi(x)
        std = np.std(isi)
        e = np.mean(isi)
        return std / e

    def cv2(self, x=None):
        isi = self.isi(x)
        summ_buffer = 0
        for i in range(0, len(isi) - 1):
            summ_buffer += 2 * (isi[i + 1] - isi[i]) / (isi[i + 1] + isi[i])
        return summ_buffer / len(isi)

    # Fano Factor
    # Warning: needs to have a trial like structure
    def ff(self, x=None, start=None, stop=None, factor=1000):
        if x is None:
            x = self.spk
        if start is None:
            start = self.start
        if stop is None:
            stop = self.stop
        x = np.array(list(map(len, x))) * factor / (stop - start)
        var = np.var(x)
        e = np.mean(x)
        return var / e

    def cleanSpk(self, spiketrain, start, stop):
        if start is not None:
            mini = np.searchsorted(spiketrain, start, side="left")
        else:
            mini = 0
        if stop is not None:
            maxi = np.searchsorted(spiketrain, stop, side="right")
        else:
            maxi = None

        if maxi is not None:
            return spiketrain[mini:maxi]
        else:
            return spiketrain[mini:]

    def sort(self, x):
        x.sort()
        return x