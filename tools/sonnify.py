import numpy as np

from scipy.io import wavfile

from obspy import Trace


def tr_to_wav(
    tr: Trace, 
    multiplier: int = 10, 
    sample_rate: float = 10000,
    normalizer: float = None,
    outfile: str = "tr_to_wav.wav"
):
    tr = tr.copy()
    # Speed-up
    tr.stats.delta /= multiplier
    tr.interpolate(sample_rate)
    
    normalizer = normalizer or np.std(tr.data)
    tr.data /= normalizer
    tr.data *= np.iinfo(np.int16).max
    tr.data = tr.data.astype(np.int16)  # Enforce 16-Bit audio
    wavfile.write(outfile, sample_rate, tr.data)


if __name__ == "__main__":
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime

    client = Client("GEONET")
    tr_raw = client.get_waveforms(
        network="NZ", station="GVZ", location="10", channel="HHZ",
        starttime=UTCDateTime(2016, 11, 13, 11, 2, 10),
        endtime=UTCDateTime(2016, 11, 13, 11, 3, 10)
    ).merge()[0]

    tr = tr_raw.copy()

    tr.detrend().taper(0.2).filter("highpass", freq=2)
    # tr.trim(
    #     UTCDateTime(2016, 11, 13, 11, 2, 10), 
    #     UTCDateTime(2016, 11, 13, 11, 3, 4, 900000)
    #     # UTCDateTime(2016, 11, 13, 11, 3, 4, 89)
    # )
    # # tr.data ** 2

    tr_to_wav(
        tr, 
        multiplier=5,
        sample_rate = 1500,
        normalizer = 400,
        outfile=f"{tr.id}_kaikoura_foreshock.wav")