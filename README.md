# pulsarsurveyscraper_server

## requirements:
* flask, flask_wtf, astropy, pandas, scipy, bs4
* [PyGEDM](https://github.com/telegraphic/pygedm): for DM models
* [pulsarsurveyscraper](https://github.com/dlakaplan/pulsarsurveyscraper)

## configuration:
Read from environment variables in `server_config.py`:
* `SECRET_KEY`: created using something like `python -c 'import os; print(os.urandom(16))'`
* `SURVEY_DATA`: location of survey cache files (default is current directory)

## API:
Queries can come via forms, or via an API that returns a JSON payload.  

* To query the survey scraper, the URL should be of the form:
```
https://pulsar.cgca-hub.org/api?type=search&ra=123.4&dec=56.7&radius=5.0
```
The usage is:
```
parameters:

    type: str ("compute" or "search")

    if type == "search":
    ra: float in degrees
    dec: float in degrees
    or
    l: float in degrees
    b: float in degrees

    radius: float in degrees
    dm: float (optional)
    dmtol: float (optional)
```
and the payload looks like:
```
{"searchra": {"display_name": "Search RA (deg)", "value": 123.4}, "searchdec": {"display_name": "Search Dec (deg)", "value": 56.7}, "searchcoord": {"display_name": "Search Coord", "value": "08:13:36 +56:42:00"}, "searchrad": {"display_name": "Search Radius (deg)", "value": 5.0}, "nmatches": 6, "J0750+57": {"ra": {"display_name": "RA (deg)", "value": 117.49999999999999}, "dec": {"display_name": "Dec (deg)", "value": 57.0}, "period": {"display_name": "Spin Period (ms)", "value": 1174.875}, "dm": {"display_name": "DM (pc/cc)", "value": 27.0}, "survey": {"display_name": "Survey", "value": "ATNF"}, "url": {"display_name": "Survey URL", "value": "https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.64&Name=Name&RaJ=RaJ&DecJ=DecJ&P0=P0&DM=DM&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Short+without+errors&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=51&table_bottom.y=23"}, "distance": {"display_name": "Distance (deg)", "value": 3.2392066983588004}, "date": {"display_name": "Retrieval Date", "value": "2021-04-01 00:00"}}, "J0749+57": {"ra": {"display_name": "RA (deg)", "value": 117.24999999999999}, "dec": {"display_name": "Dec (deg)", "value": 57.0}, "period": {"display_name": "Spin Period (ms)", "value": 1175.0}, "dm": {"display_name": "DM (pc/cc)", "value": 27.0}, "survey": {"display_name": "Survey", "value": "GBT820"}, "url": {"display_name": "Survey URL", "value": "http://astro.phys.wvu.edu/GBNCC/"}, "distance": {"display_name": "Distance (deg)", "value": 3.3752178169661717}, "date": {"display_name": "Retrieval Date", "value": "2021-04-29 00:00"}}, "J0746+55": {"ra": {"display_name": "RA (deg)", "value": 116.69583333333333}, "dec": {"display_name": "Dec (deg)", "value": 55.2}, "period": {"display_name": "Spin Period (ms)", "value": 2893.4700000000003}, "dm": {"display_name": "DM (pc/cc)", "value": 10.5}, "survey": {"display_name": "Survey", "value": "CHIME"}, "url": {"display_name": "Survey URL", "value": "http://catalog.chime-frb.ca/galactic"}, "distance": {"display_name": "Distance (deg)", "value": 4.040256384585515}, "date": {"display_name": "Retrieval Date", "value": "2021-04-29 00:00"}}, "J0827+53": {"ra": {"display_name": "RA (deg)", "value": 126.94999999999997}, "dec": {"display_name": "Dec (deg)", "value": 53.0}, "period": {"display_name": "Spin Period (ms)", "value": 13.5}, "dm": {"display_name": "DM (pc/cc)", "value": 23.1}, "survey": {"display_name": "Survey", "value": "ATNF"}, "url": {"display_name": "Survey URL", "value": "https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.64&Name=Name&RaJ=RaJ&DecJ=DecJ&P0=P0&DM=DM&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Short+without+errors&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=51&table_bottom.y=23"}, "distance": {"display_name": "Distance (deg)", "value": 4.225635766208534}, "date": {"display_name": "Retrieval Date", "value": "2021-04-01 00:00"}}}
```
with data:
```
searchra
searchdec
searchcoord
searchrad
searchdm (if specified)
searchdmtolerance (if specified
nmatches
```
and then for every match:
```
display_name
ra
dec
survey
distance
period
dm
url
```

* To query the DM tool, the URL should be of the form:
```
https://pulsar.cgca-hub.org/api?type=compute&ra=123.4&dec=56.7&dm=30.0&dmmodel=ymw16
```
The usage is:
```
parameters:

    type: str ("compute" or "search")

    if type == "compute":
    ra: float in degrees
    dec: float in degrees
    or
    l: float in degrees
    b: float in degrees

    dmmodel: str ("ne2001" or "ymw16")

    one of:
    d: float
    dm: float       
```
and the payload looks like:
```
{"searchra": {"display_name": "Search RA (deg)", "value": 123.4}, "searchdec": {"display_name": "Search Dec (deg)", "value": 56.7}, "searchcoord": {"display_name": "Search Coord", "value": "123.4 56.7"}, "dmmodel": {"display_name": "DM Model", "value": "YMW16"}, "searchdm": {"display_name": "Search DM", "value": 30.0}, "computed_d": {"display_name": "Computed Distance (pc)", "value": 2139.48095703125}, "max_dm": {"display_name": "Max DM", "value": 44.41266632080078}}
```
with data:
```
searchra
searchdec
searchcoord
dmmodel
searchdm (if specified)
or
searchd (if specified)
max_dm
computed_d
or
computed_dm
```
## deployment:
Based on [this tutorial](https://www.digitalocean.com/community/tutorials/how-to-serve-flask-applications-with-gunicorn-and-nginx-on-ubuntu-18-04).  Also used [this tutorial](https://github.com/jupyterhub/the-littlest-jupyterhub/issues/272) to keep TLJH running behind a proxy.

To restart the server: `sudo systemctl restart pulsarsurveyscraper`
