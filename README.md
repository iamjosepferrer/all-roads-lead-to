# All roads lead to…

An interactive map of the fastest driving routes to any city across the EU, the UK, EFTA and the Western Balkans. Pick a destination and every main road in the network routes to it, with line width showing how many people use each road on their way there.

This first part is the data pipeline: it turns OpenStreetMap roads and a population grid into one compact routing graph that a browser can load. The map itself comes next.

## How it works

**1. Extract** (`pipeline/extract.py`). DuckDB reads the public [Overture Maps](https://overturemaps.org) transportation data straight from S3 and keeps motorways, trunk roads, primary roads and ferry routes inside the bounding box: about 2.7 million segments.

**2. Build** (`pipeline/build.py`, about 2 minutes).
- Keeps the countries in scope; bridges such as Øresund and the Great Belt are kept even though their midpoints are over water.
- Sets a speed per road: the signed limit where there is one, otherwise 120, 100 and 80 km/h by road class and 60 on slip roads. Ferries run at 35 km/h with 20 minutes to board at each end; the Channel Tunnel shuttle takes 65 minutes.
- Bridges short gaps left by the missing minor roads with straight links at 40 km/h, links ferry terminals to the nearest road, drops cruise-like ferry loops, and keeps the largest connected network.
- Removes short dead ends and merges chains of road pieces into single links (661k junctions, 1.0M links).
- Gives every road the people who live closest to it, from the GHSL 2025 population grid (538 million people in total).
- Packs everything into `data/graph.bin.gz` (6.5 MB): node coordinates in the European LAEA projection, link times, lengths and people, and simplified geometry, stored as byte-shuffled columns so gzip compresses them well. `data/meta.json` holds the cities and country outlines.

**3. Check** (`pipeline/check.py`). Routes from Brussels with and without ferries:

```
from Brussels     ferries  roads only
London               4.1 h       4.1 h
Warsaw              11.0 h      11.0 h
Lisbon              17.4 h      17.4 h
Stockholm           14.1 h      14.2 h
Helsinki            23.5 h      33.3 h
Dublin              11.5 h     no road
people cut off without ferries: 3.6%
```

## Rebuild the data

```bash
pip install -r pipeline/requirements.txt
bash pipeline/download.sh      # Natural Earth + GHSL population into data_raw/
python pipeline/extract.py     # road network from Overture into data_raw/ (~490 MB)
python pipeline/build.py       # writes data/graph.bin.gz and data/meta.json
python pipeline/check.py       # drive times from Brussels as a sanity check
```

`build.py` peaks at about 6 GB of RAM. To use a newer Overture release, change `RELEASE` at the top of `extract.py`.

## Repository layout

```
data/
  graph.bin.gz        routing graph (built)
  meta.json           cities and country outlines (built)
pipeline/
  download.sh         static inputs
  extract.py          road network from Overture
  build.py            graph builder
  check.py            sanity check
  requirements.txt
```

## Data and credits

- Roads and ferries: © [OpenStreetMap](https://www.openstreetmap.org/copyright) contributors, via [Overture Maps](https://docs.overturemaps.org/attribution/). Open Database License (ODbL).
- Population: [GHSL GHS-POP R2023A](https://human-settlement.emergency.copernicus.eu/), epoch 2025, European Commission Joint Research Centre. CC BY 4.0.
- Country outlines and cities: [Natural Earth](https://www.naturalearthdata.com/). Public domain.
- Idea: [Todos los caminos llevan a…](https://ntizar.github.io/CuencasTrafico/) (Spain) and Arterial Mapping's *All Roads Lead to Richmond*.

Made by Josep Ferrer.
