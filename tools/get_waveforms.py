"""
Get the waveforms and correct to displacement for the East Cape earthquake
at Wellington and at Mt Maunganui
"""

from typing import Tuple

from obspy.clients.fdsn import Client
from obspy import Trace, Stream
from obspy.geodetics import degrees2kilometers, gps2dist_azimuth


EVENTID = "2021p169083"
STATIONS = [("WEL", "HHZ"), ("OPRZ", "HHZ")]  #TBCS is strong motion
CLIENT = Client("GEONET")

PAZ_WA = {'sensitivity': 2080, 'zeros': [0j], 'gain': 1,
          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}


def get_and_correct(
    station: str, 
    eventid: str, 
    channel: str,
    duration: float = 360.
) -> Tuple[Trace, float]:
    ev = CLIENT.get_events(eventid=EVENTID)[0]

    origin = ev.preferred_origin() or ev.origins[-1]
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
        lat1=origin.latitude, lon1=origin.latitude,
        lat2=inv[0][0].latitude, lon2=inv[0][0].longitude)
    distance = distance / 1000 
    return tr, distance
    
    
def plot_responses():
    info = dict()
    for station, channel in STATIONS:
        tr, distance = get_and_correct(station, EVENTID, channel)
        info.update({station: {"trace": tr, "distance": distance}})