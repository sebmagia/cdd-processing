import obspy
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL

# We'll first create all the various objects.
# These strongly follow the hierarchy of StationXML files.
inv = Inventory(networks=[], source="ObsPy-Tutorial")
net = Network(code="XX",
              stations=[],
              description="A test station.",
              start_date=obspy.UTCDateTime(2016, 1, 2))
sta = Station(code="ABC",
              latitude=1.0,
              longitude=2.0,
              elevation=345.0,
              creation_date=obspy.UTCDateTime(2016, 1, 2),
              site=Site(name="First station"))
cha = Channel(code="HHZ",
              location_code="",
              latitude=1.0,
              longitude=2.0,
              elevation=345.0,
              depth=10.0,
              azimuth=0.0,
              dip=-90.0,
              sample_rate=200)

# By default this accesses the NRL online.
nrl = NRL()
# We assume that the end point of data logger and sensor are known:
response = nrl.get_response(sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
                            datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])
cha.response = response   #tie all together
sta.channels.append(cha)
net.stations.append(sta)
inv.networks.append(net)
# write to stationXML file & force validation against the StationXML schema
inv.write("new_station.xml", format="stationxml", validate=True)