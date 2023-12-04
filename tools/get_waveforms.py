"""
Get the waveforms and correct to displacement for the East Cape earthquake
at Wellington and at Mt Maunganui
"""

from typing import Tuple

from obspy.clients.fdsn import Client
from obspy import Trace, Stream
from obspy.core.event import Event
from obspy.geodetics import degrees2kilometers, gps2dist_azimuth

import matplotlib.pyplot as plt


# Currently hard-coded to use GeoNet for simplicity.
CLIENT = Client("GEONET")

PAZ_WA = {'sensitivity': 2080, 'zeros': [0j], 'gain': 1,
          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}


def get_and_plot_waveforms(
    eventid: str, 
    equal_scale: bool = False,
    duration_multiplier: float = 2.0,
) -> plt.Figure:
    try:
        ev = CLIENT.get_events(eventid=eventid)[0]
    except Exception as e:
        print(f"Could not download {eventid} due to {e}")
        raise e
    stations = sorted(list({p.waveform_id.station_code for p in ev.picks}))
    print(f"Available stations:\n{' '.join([s for s in stations])}")
    print("Type the code(s) of the stations you want to plot. "
          "Separate stations with spaces (e.g. FOZ WVZ). "
          "Press return when done")
    stations = input("Enter stations: ")
    stations = stations.split(" ")
    print(f"You selected: {' '.join([s for s in stations])}")
    
    # Work out "optimal" duration
    vpvs = 1.76  # Assumed vpvs ratio
    last_p = max([p.time for p in ev.picks 
                  if p.waveform_id.station_code in stations 
                  and p.phase_hint.startswith("P")])
    origin = _get_origin(ev)
    p_tt = last_p - origin.time
    s_tt = p_tt * vpvs
    duration = s_tt * duration_multiplier  # Hardcoded bodge
    
    st = Stream()
    for station in stations:
        # Work out what channels are available
        picks = [p for p in ev.picks if p.waveform_id.station_code == station]
        if len(picks) == 0:
            print(f"Skipping station {station} due to no matching picks")
            continue
        codes = {p.waveform_id.channel_code[0:2] for p in picks}
        # Preference for HHZ, the EHZ, then HNZ
        if "HH" in codes:
            channel = "HHZ"
        elif "EH" in codes:
            channel = "EHZ"
        elif "HN" in codes:
            channel = "HNZ"
        else:
            channel = codes.pop()
        tr, distance = get_and_correct(
            station=station, channel=channel, duration=duration,
            event=ev)
        tr.stats.distance = distance
        st += tr
    
    fig = _plot_waveforms(st=st, event=ev, equal_scale=equal_scale)
    return fig


def _plot_waveforms(st: Stream, event: Event, equal_scale: bool = False) -> plt.Figure():
    for tr in st:
        assert hasattr(tr.stats, "distance"), "Requires trace.stats to have distances"
    origin = _get_origin(event)
    
    fig, axes = plt.subplots(len(st), 1, sharex=True, sharey=equal_scale, squeeze=False)
    axes = axes[:,0]
    st.traces.sort(key=lambda tr: tr.stats.distance)
    for tr, ax in zip(st, axes):
        times = tr.times(type="relative", reftime=origin.time)
        ax.plot(times, tr.data, color="k")
        ax.text(0, 0.85, tr.id, transform=ax.transAxes)
        ax.text(0,0.1, f"Distance: {tr.stats.distance:.2f} km", transform=ax.transAxes)
        ax.grid("on")
        ax.set_ylabel("Wood-Anderson\nAmplitude (m)")
    
    ax.set_xlabel("Seconds since origin time")
        
    fig.subplots_adjust(hspace=0, wspace=0)
    return fig
    

    
def _get_origin(event: Event):
    return event.preferred_origin() or event.origins[-1]
                 
                 

def get_and_correct(
    station: str, 
    channel: str,
    duration: float = 360.,
    event: Event = None,
    eventid: str = None,
) -> Tuple[Trace, float]:
    
    assert eventid or event, "Requires either event or eventid"
    
    ev = event or CLIENT.get_events(eventid=EVENTID)[0]

    origin = _get_origin(event)
    st = CLIENT.get_waveforms(
        network="NZ", station=station, location="*", channel=channel,
        starttime=origin.time, endtime=origin.time + duration)
    tr = st.merge()[0]
    inv = CLIENT.get_stations(
        network="NZ", station=station, location="*", channel=channel,
        startbefore=origin.time, endafter=origin.time + duration,
        level="response")
    
    tr = tr.remove_response(inventory=inv, output="DISP", 
                            pre_filt=(0.005, 0.006, 30.0, 35.0))
    tr.simulate(paz_simulate=PAZ_WA, water_level=10)[0]

    
    distance, _, _ = gps2dist_azimuth(
        lat1=origin.latitude, lon1=origin.longitude,
        lat2=inv[0][0].latitude, lon2=inv[0][0].longitude)
    distance = distance / 1000 
    return tr, distance
    
    
def plot_responses():
    info = dict()
    for station, channel in STATIONS:
        tr, distance = get_and_correct(station, EVENTID, channel)
        info.update({station: {"trace": tr, "distance": distance}})
