from flask import Flask, render_template, request, url_for, redirect, send_file
import flask
from forms import (
    SearchForm,
    DMForm,
    parse_equcoord_and_validate,
    parse_galcoord_and_validate,
    lb_label,
    radec_label,
)

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
from astropy import units as u
import pandas
import json
import os

from bs4 import BeautifulSoup, Tag

import pulsarsurveyscraper
import pulsarsurveyscraper.output

import pygedm

degree_symbol = "\N{DEGREE SIGN}"
arcmin_symbol = "\N{PRIME}"
arcsec_symbol = "\N{DOUBLE PRIME}"

"""
Among others, I used these tutorials:
https://codeloop.org/flask-tutorial-flask-forms-with-flask-wtf/
"""

# create the Flask object
app = Flask(__name__)

app.config.from_object("server_config")
# turn off the automatic sorting of the output
app.config["JSON_SORT_KEYS"] = False
# instantiate the pulsar survey table
# this needs the directory where the data are stored
# do it at the top level so that it's accessible within the functions below
pulsar_table = pulsarsurveyscraper.PulsarTable(
    directory=app.config["DATA_DIR"],
)
s, n = np.unique(pulsar_table.data["survey"], return_counts=True)
# this is a basic dictionary to be used
# for posting on the form (survey names/links)
survey_data = {}
for survey, number in zip(s, n):
    survey_data[survey] = {
        "url": pulsarsurveyscraper.Surveys[survey]["url"],
        "number": number,
        "date": pulsar_table.data[pulsar_table.data["survey"] == survey][0][
            "retrieval date"
        ],
    }

coordinate_type = "equatorial"


def format_lb(coord):
    return (
        "{}{}".format(
            coord.galactic.l.to_string(
                decimal=True,
                precision=3,
            ),
            degree_symbol,
        ),
        "{}{}".format(
            coord.galactic.b.to_string(
                decimal=True,
                alwayssign=True,
                precision=3,
            ),
            degree_symbol,
        ),
    )


def format_radec_decimal(coord):
    return (
        "{}{}".format(
            coord.icrs.ra.to_string(
                decimal=True,
                precision=3,
            ),
            degree_symbol,
        ),
        "{}{}".format(
            coord.icrs.dec.to_string(decimal=True, alwayssign=True, precision=3),
            degree_symbol,
        ),
    )


def format_radec(coord):
    sra = coord.icrs.ra.to_string(u.hour, decimal=False, sep="hms", precision=2)
    sdec = coord.icrs.dec.to_string(
        u.degree, decimal=False, sep="dms", precision=1, pad=True, alwayssign=True
    )
    sra = (
        sra.replace("s", "<sup>s</sup>")
        .replace("h", "<sup>h</sup>")
        .replace("m", "<sup>m</sup>")
    )
    sdec = (
        sdec.replace("d", degree_symbol)
        .replace("m", arcmin_symbol)
        .replace("s", arcsec_symbol)
    )
    return "{}, {}".format(sra, sdec)


def format_time(t, digits=2):
    """return a string with a time formatted properly (allow ms, us, ...)"""
    # seems like there has to be a better way
    exponent = np.log10(t.to_value(u.s))
    format_type = "f"
    if exponent >= 0:
        unit = u.s
        if exponent >= 3:
            # don't do ks, MS, ...
            format_type = "e"
    elif exponent >= -3:
        unit = u.ms
    elif exponent >= -6:
        unit = u.us
    elif exponent >= -9:
        unit = u.ns
    else:
        unit = u.s
        format_type = "e"
    format = f"{{:.{digits}{format_type}}} {unit}"
    return format.format(t.to_value(unit))


def format_frequency(f, digits=2):
    """return a string with a time formatted properly (allow kHz, MHz, ...)"""
    exponent = np.log10(f.to_value(u.Hz))
    format_type = "f"
    if exponent < 0:
        unit = u.Hz
        format_type = "e"
    elif exponent >= 9:
        unit = u.GHz
    elif exponent >= 6:
        unit = u.MHz
    elif exponent >= 3:
        unit = u.kHz
    else:
        unit = u.Hz
    format = f"{{:.{digits}{format_type}}} {unit}"
    return format.format(f.to_value(unit))


@app.route("/", methods=["GET", "POST"])
@app.route("/search", methods=["GET", "POST"])
def Search():
    """
    search form route

    """
    form = SearchForm()

    # this is a basic dictionary to be used
    # for posting on the form (survey names/links)
    surveys = {}
    for survey in pulsarsurveyscraper.Surveys:
        surveys[survey] = pulsarsurveyscraper.Surveys[survey]["url"]

    # if the button has been pressed and the input is valid:
    if form.validate_on_submit():
        # use the validation routine
        # to parse the coordinates
        if form.lb_or_radec.data:
            coord = parse_equcoord_and_validate(None, form.coordinates)
        else:
            coord = parse_galcoord_and_validate(None, form.coordinates)
        # figure out if we need to deal with DM as well
        if form.dm.data is not None:
            DM = float(form.dm.data)
        else:
            DM = form.dm.data
        if form.dmtol.data is not None:
            DMtol = float(form.dmtol.data)
        else:
            DMtol = form.dm.data

        if form.deduplicate.data:
            deduplicate = "hide"
        else:
            deduplicate = False

        # if we've pressed the "API" button
        # instead of showing the tabular output redirect
        # to the API result
        if form.api.data:
            if DM is None:
                if form.lb_or_radec.data:
                    return redirect(
                        url_for(
                            "API",
                            type="search",
                            ra=coord.icrs.ra.deg,
                            dec=coord.icrs.dec.deg,
                            radius=float(form.radius.data),
                        )
                    )
                else:
                    return redirect(
                        url_for(
                            "API",
                            type="search",
                            l=coord.galactic.l.deg,
                            b=coord.galactic.b.deg,
                            radius=float(form.radius.data),
                        )
                    )

            else:
                if form.lb_or_radec.data:
                    return redirect(
                        url_for(
                            "API",
                            type="search",
                            ra=coord.icrs.ra.deg,
                            dec=coord.icrs.dec.deg,
                            radius=float(form.radius.data),
                            dm=DM,
                            dmtol=DMtol,
                        )
                    )
                else:
                    return redirect(
                        url_for(
                            "API",
                            type="search",
                            l=coord.galactic.l.deg,
                            b=coord.galactic.b.deg,
                            radius=float(form.radius.data),
                            dm=DM,
                            dmtol=DMtol,
                        )
                    )

        # or, clear the form if desired
        if form.clear.data:
            return redirect(url_for("Search"))

        # if we've made it this far, we want the basic search
        # including HTML table output

        # first get the astropy Table
        if form.lb_or_radec.data:
            result = pulsar_table.search(
                coord,
                radius=float(form.radius.data) * u.deg,
                DM=DM,
                DM_tolerance=DMtol,
                deduplicate=deduplicate,
            )
        else:
            result = pulsar_table.search(
                coord,
                radius=float(form.radius.data) * u.deg,
                DM=DM,
                DM_tolerance=DMtol,
                return_native=True,
                deduplicate=deduplicate,
            )

        result["P"][result["P"] < 0] = np.nan
        # make a nice string for output
        if form.lb_or_radec.data:
            coord_string = "Searching <strong>{:.1f}{}</strong> around <strong>{} = {} = {}, {}</strong> ...".format(
                float(form.radius.data),
                degree_symbol,
                radec_label.replace(" ", ","),
                format_radec(coord),
                *format_radec_decimal(coord),
            )
        else:
            coord_string = "Searching <strong>{:.1f}{}</strong> around <strong>{} = {},{}</strong> ...".format(
                float(form.radius.data),
                degree_symbol,
                lb_label.replace(" ", ", "),
                *format_lb(coord),
            )

        if DM is not None:
            coord_string += (
                "<br>Also requiring DM = <strong>{:.1f}+/-{:.1f} pc/cc</strong>".format(
                    DM,
                    DMtol,
                )
            )

        if form.PNG.data or form.PDF.data:
            # we want a PDF or PNG output
            # Get the global URL for the query web page:
            # do it without actually querying
            if DM is None:
                if form.lb_or_radec.data:
                    query_url = url_for(
                        "API",
                        ra=coord.icrs.ra.deg,
                        dec=coord.icrs.dec.deg,
                        radius=float(form.radius.data),
                    )

                else:
                    query_url = url_for(
                        "API",
                        l=coord.galactic.l.deg,
                        b=coord.galactic.b.deg,
                        radius=float(form.radius.data),
                    )

            else:
                if form.lb_or_radec.data:
                    query_url = url_for(
                        "API",
                        ra=coord.icrs.ra.deg,
                        dec=coord.icrs.dec.deg,
                        radius=float(form.radius.data),
                        dm=DM,
                        dmtol=DMtol,
                    )

                else:
                    query_url = url_for(
                        "API",
                        l=coord.galactic.l.deg,
                        b=coord.galactic.b.deg,
                        radius=float(form.radius.data),
                        dm=DM,
                        dmtol=DMtol,
                    )

            format = "pdf" if form.PDF.data else "png"
            if form.lb_or_radec.data:
                search_query_txt = (
                    "Searching {:.1f}deg around RA,Dec = {} = {}d,{}d".format(
                        float(form.radius.data),
                        coord.to_string("hmsdms", sep=":"),
                        coord.ra.to_string(decimal=True),
                        coord.dec.to_string(decimal=True, alwayssign=True),
                    )
                )

            else:
                search_query_txt = "Searching {:.1f}deg around l,b = {}d,{}d".format(
                    float(form.radius.data),
                    coord.l.to_string(decimal=True),
                    coord.b.to_string(decimal=True, alwayssign=True),
                )

            output = pulsarsurveyscraper.output.make_output(
                result,
                format,
                search_query_txt,
                request.base_url + query_url,
                directory=app.config["OUTPUT_DIR"],
            )
            return send_file(output, as_attachment=True)

        # go from astropy Table -> pandas dataframe -> HTML table
        df = result.to_pandas()
        if len(df) > 0 and not isinstance(df["PSR"][0], str):
            # turn the "PSR" column from bytestring to string
            df["PSR"] = df["PSR"].str.decode("utf-8")
        html_table = df.to_html(
            formatters={
                "P": lambda x: "%.2f" % x,
                "Distance": lambda x: "%.2f" % x,
            },
            justify="left",
        )

        soup = BeautifulSoup(html_table, "html.parser")
        header_row = soup.find("tr")
        header_cols = header_row.find_all("th")
        # add units to the header row
        for col in header_cols:
            if col.text == "P":
                col.string = "P (ms)"
            elif col.text == "Distance":
                col.string = "Distance (deg)"
            elif col.text == "RA":
                col.string = "RA (deg)"
            elif col.text == "Dec":
                col.string = "Dec (deg)"
            elif col.text == "l":
                col.string = "l (deg)"
            elif col.text == "b":
                col.string = "b (deg)"
        # fix the alignment of various columns
        col_aligns = {3: "right", 4: "right", 5: "center", 7: "right"}
        rows = soup.find_all("tr")
        for row in rows[1:]:
            cols = row.find_all("td")
            # add links to survey column
            link_tag = soup.new_tag(
                "a",
                href=pulsarsurveyscraper.Surveys[cols[5].text]["url"],
            )
            link_tag.string = cols[5].text
            cols[5].string = ""
            cols[5].insert(0, link_tag)
            if cols[5].text == "ATNF":
                # link through to the catalog page
                link = f"https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.70&Name=Name&RaJ=RaJ&DecJ=DecJ&P0=P0&P1=P1&DM=DM&Survey=Survey&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names={cols[0].text}&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+last+digit+error&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=46&table_bottom.y=24"
                link_tag = soup.new_tag(
                    "a",
                    href=link,
                )
                link_tag.string = cols[0].text
                cols[0].string = ""
                cols[0].insert(0, link_tag)

            for i in col_aligns:
                cols[i]["align"] = col_aligns[i]
        html_table = soup

        return render_template(
            "search.html",
            form=form,
            surveys=surveys.items(),
            coord_string=coord_string,
            nfound=len(result),
            html_table=html_table,
        )

    # we can clear even if it's not valid
    if form.clear.data:
        return redirect(url_for("Search"))

    return render_template("search.html", form=form, surveys=surveys.items())


@app.route("/api", methods=["GET"])
def API():
    """
    API route

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

    returns JSON
    """

    type = "search"

    if "type" in request.args:
        type = request.args["type"].lower()

    # search for pulsars given a position
    if type == "search":
        # default values
        radius = 5
        dm = None
        dmtol = 10

        ra = None
        dec = None
        l = None
        b = None

        if "ra" in request.args:
            ra = float(request.args["ra"])
        if "dec" in request.args:
            dec = float(request.args["dec"])
        if "l" in request.args:
            l = float(request.args["l"])
        if "b" in request.args:
            b = float(request.args["b"])
        if "radius" in request.args:
            radius = float(request.args["radius"])

        coord = None

        if ra is not None and dec is not None:
            try:
                coord = SkyCoord(ra * u.deg, dec * u.deg)
            except ValueError as e:
                return "Unable to parse RA,Dec = '{},{}': {}".format(ra, dec, e)
        elif l is not None and b is not None:
            try:
                coord = SkyCoord(l * u.deg, b * u.deg, frame="galactic")
            except ValueError as e:
                return "Unable to parse l,b = '{},{}': {}".format(l, b, e)

        if coord is None:
            return "Error: must specify RA,Dec or l,b"

        if "dm" in request.args:
            dm = float(request.args["dm"])
        if "dmtol" in request.args:
            dmtol = float(request.args["dmtol"])

        result = pulsar_table.search(
            coord,
            radius=radius * u.deg,
            DM=dm,
            DM_tolerance=dmtol,
            return_json=True,
            return_native=True,
            deduplicate=True,
        )

    # compute DM at a distance
    elif type == "compute":
        dm = None
        d = None
        d = None
        dmmodel = "ne2001"

        ra = None
        dec = None
        l = None
        b = None

        if "ra" in request.args:
            ra = float(request.args["ra"])
        if "dec" in request.args:
            dec = float(request.args["dec"])
        if "l" in request.args:
            l = float(request.args["l"])
        if "b" in request.args:
            b = float(request.args["b"])

        coord = None

        if ra is not None and dec is not None:
            try:
                coord = SkyCoord(ra * u.deg, dec * u.deg)
            except ValueError as e:
                return "Unable to parse RA,Dec = '{},{}': {}".format(ra, dec, e)
            result = {
                "searchra": {
                    "display_name": "Search RA (deg)",
                    "value": coord.ra.value,
                },
                "searchdec": {
                    "display_name": "Search Dec (deg)",
                    "value": coord.dec.value,
                },
            }
        elif l is not None and b is not None:
            try:
                coord = SkyCoord(l * u.deg, b * u.deg, frame="galactic")
            except ValueError as e:
                return "Unable to parse l,b = '{},{}': {}".format(l, b, e)
            result = {
                "searchl": {"display_name": "Search l (deg)", "value": coord.l.value},
                "searchb": {"display_name": "Search b (deg)", "value": coord.b.value},
            }

        if coord is None:
            return "Error: must specify RA,Dec or l,b"

        result["searchcoord"] = {
            "display_name": "Search Coord",
            "value": coord.to_string(),
        }

        if "dm" in request.args:
            dm = float(request.args["dm"])
        if "d" in request.args:
            d = float(request.args["d"]) * u.pc
        if "dmmodel" in request.args:
            dmmodel = request.args["dmmodel"].lower()
            if not dmmodel in ["ymw16", "ne2001"]:
                return "Error: '{}' is not a valid DM model".format(dmmodel)
            result["dmmodel"] = dmmodel

        if dm is not None:
            model_d, _ = pygedm.dm_to_dist(
                coord.galactic.l, coord.galactic.b, dm, method=dmmodel.upper()
            )
            result["searchdm"] = {"display_name": "Search DM", "value": dm}
            result["computed_d"] = {
                "display_name": "Computed Distance (pc)",
                "value": model_d.to(u.pc).value,
            }
        elif d is not None:
            model_dm, _ = pygedm.dist_to_dm(
                coord.galactic.l,
                coord.galactic.b,
                d,
                method=dmmodel.upper(),
            )
            result["search_d"] = {
                "display_name": "Search D (pc)",
                "value": d.to(u.pc).value,
            }
            result["computed_dm"] = {
                "display_name": "Computed DM",
                "value": model_dm.value,
            }
        max_DM, _ = pygedm.dist_to_dm(
            coord.galactic.l,
            coord.galactic.b,
            100 * u.kpc,
            method=dmmodel.upper(),
        )
        result["max_dm"] = {"display_name": "Max DM", "value": max_DM.value}
        result["dmmodel"] = {"display_name": "DM Model", "value": dmmodel.upper()}
        # convert dict to JSON
        result = json.dumps(result)

    return result


@app.route("/compute", methods=["GET", "POST"])
def Compute():
    """
    compute form route

    """
    form = DMForm()

    # if the button has been pressed and the input is valid:
    if form.validate_on_submit():
        if form.model_selector.data:
            model = "NE2001"
        else:
            model = "YMW16"
        # use the validation routine
        # to parse the coordinates
        if form.lb_or_radec.data:
            coord = parse_equcoord_and_validate(None, form.coordinates)
        else:
            coord = parse_galcoord_and_validate(None, form.coordinates)
        if not form.d_or_dm_selector.data:
            DM = float(form.d_or_dm.data)
            distance, _ = pygedm.dm_to_dist(
                coord.galactic.l, coord.galactic.b, DM, method=model
            )
            Array = pygedm.ne2001_wrapper.dm_to_dist(
                coord.galactic.l.value,
                coord.galactic.b.value,
                DM,
                form.frequency.data * u.GHz,
                full_output=True,
            )
            nu = form.frequency.data * u.GHz
            V_perp = form.velocity.data * u.km / u.s
            sm = Array["smtau"]
            tauiss = Array["tau_sc"].value * u.s
            c1 = 1.16  # param for uniform, Kolmogorov medium

            # calculate scintillation bandwidth (scintbw, kHz)
            # Ref: Eq. 35 from Cordes & Lazio 1991
            scintbw = (c1 / (2 * np.pi * tauiss.to_value(u.ms))) * u.kHz
            # calculate scintillation time scale (scinttime, s)
            # Ref: Eq. 46 from Cordes & Lazio 1991 (2.3 change to 3.3)
            scinttime = (
                (3.3 * u.s)
                * nu.to_value(u.GHz) ** 1.2
                * sm ** (-0.6)
                * (100 / V_perp.to_value(u.km / u.s))
            )
            # and scale scattering time for frequency as well
            # this exponent is used in NE2001 TAUISS()
            tauiss *= nu.to_value(u.GHz) ** -4.4
        else:
            distance = float(form.d_or_dm.data) * u.pc
            DM, _ = pygedm.dist_to_dm(
                coord.galactic.l,
                coord.galactic.b,
                distance,
                method=model,
            )
            Array = pygedm.ne2001_wrapper.dist_to_dm(
                coord.galactic.l.value,
                coord.galactic.b.value,
                distance.to_value(u.kpc),
                form.frequency.data * u.GHz,
                full_output=True,
            )
            nu = form.frequency.data * u.GHz
            V_perp = form.velocity.data * u.km / u.s
            sm = Array["smtau"]
            tauiss = Array["tau_sc"].value * u.s
            c1 = 1.16  # param for uniform, Kolmogorov medium

            # calculate scintillation bandwidth (scintbw, kHz)
            # Ref: Eq. 35 from Cordes & Lazio 1991
            scintbw = (c1 / (2 * np.pi * tauiss.to_value(u.ms))) * u.kHz
            # calculate scintillation time scale (scinttime, s)
            # Ref: Eq. 46 from Cordes & Lazio 1991 (2.3 change to 3.3)
            scinttime = (
                (3.3 * u.s)
                * nu.to_value(u.GHz) ** 1.2
                * sm ** (-0.6)
                * (100 / V_perp.to_value(u.km / u.s))
            )
            # and scale scattering time for frequency as well
            # this exponent is used in NE2001 TAUISS()
            tauiss *= nu.to_value(u.GHz) ** -4.4

        max_DM, _ = pygedm.dist_to_dm(
            coord.galactic.l,
            coord.galactic.b,
            100 * u.kpc,
            method=model,
        )

        # if we've pressed the "API" button
        # instead of showing the tabular output redirect
        # to the API result
        if form.api.data:
            if form.d_or_dm_selector.data:
                if form.lb_or_radec.data:
                    return redirect(
                        url_for(
                            "API",
                            type="compute",
                            ra=coord.icrs.ra.deg,
                            dec=coord.icrs.dec.deg,
                            d=distance,
                            dmmodel=model.lower(),
                        )
                    )
                else:
                    return redirect(
                        url_for(
                            "API",
                            type="compute",
                            l=coord.galactic.l.deg,
                            b=coord.galactic.b.deg,
                            d=distance,
                            dmmodel=model.lower(),
                        )
                    )

            else:
                if form.lb_or_radec.data:
                    return redirect(
                        url_for(
                            "API",
                            type="compute",
                            ra=coord.icrs.ra.deg,
                            dec=coord.icrs.dec.deg,
                            dm=DM,
                            dmmodel=model.lower(),
                        )
                    )
                else:
                    return redirect(
                        url_for(
                            "API",
                            type="compute",
                            l=coord.galactic.l.deg,
                            b=coord.galactic.b.deg,
                            dm=DM,
                            dmmodel=model.lower(),
                        )
                    )

        # or, clear the form if desired
        if form.clear.data:
            return redirect(url_for("Compute"))

        # make a nice string for output
        if form.lb_or_radec.data:
            coord_string = "Computing for <strong>{} = {} = {}, {}</strong>".format(
                radec_label.replace(" ", ","),
                format_radec(coord),
                *format_radec_decimal(coord),
            )
            coord_string += "<br>= {} = {}, {} ...".format(
                lb_label.replace(" ", ","),
                *format_lb(coord),
            )
        else:
            coord_string = "Computing for <strong>{} = {}, {}</strong>".format(
                lb_label.replace(" ", ","), *format_lb(coord)
            )
            coord_string += "<br>= {} = {} = {}, {} ...".format(
                radec_label.replace(" ", ","),
                format_radec(coord),
                *format_radec_decimal(coord),
            )
        if not form.d_or_dm_selector.data:
            result_string = "For <strong>DM = {:.1f} pc/cc</strong>, find <strong>distance = {:.1f} pc</strong> with the {} model".format(
                DM, distance.to(u.pc).value, model
            )
        else:
            result_string = "For <strong>distance = {:.1f} pc</strong>, find <strong>DM = {:.1f} pc/cc</strong> with the {} model".format(
                distance.to(u.pc).value, DM.value, model
            )
        result_string += (
            "<br>Along this LOS, <strong>max DM = {:.1f} pc/cc</strong>".format(
                max_DM.value
            )
        )
        result_string += (
            "<br>Scintillation bandwidth (at {:.2f} GHz) = <strong>{}</strong>".format(
                nu.to_value(u.GHz), format_frequency(scintbw)
            )
        )
        result_string += (
            "<br>Scattering timescale (at {:.2f} GHz) = <strong>{}</strong>".format(
                nu.to_value(u.GHz), format_time(tauiss)
            )
        )
        result_string += "<br>Scintillation timescale (for v = {:.1f} km/s) = <strong>{}</strong>".format(
            V_perp.to_value(u.km / u.s), format_time(scinttime)
        )
        return render_template(
            "compute.html",
            form=form,
            coord_string=coord_string,
            result_string=result_string,
        )

    # we can clear even if it's not valid
    if form.clear.data:
        return redirect(url_for("Compute"))

    return render_template("compute.html", form=form)


@app.route("/get_coordtoggled_status")
def toggled_status():
    coordinate_status = flask.request.args.get("status")
    # coordinate_type = "equatorial" if coordinate_status == "true" else "galactic"
    # form.lb_or_radec_selector.data = coordinate_type
    return coordinate_type


@app.route("/surveys")
def Surveys():
    data = []
    for s in survey_data.keys():
        data.append([s, survey_data[s]["number"], survey_data[s]["date"]])
    tabledata = Table(rows=data, names=("Survey", "Number", "Retrieval Date"))
    # go from astropy Table -> pandas dataframe -> HTML table
    df = tabledata.to_pandas()
    html_table = df.to_html(
        justify="left",
    )
    soup = BeautifulSoup(html_table, "html.parser")
    header_row = soup.find("tr")
    header_cols = header_row.find_all("th")
    rows = soup.find_all("tr")
    # fix the alignment of various columns
    col_aligns = {1: "right"}
    for row in rows[1:]:
        cols = row.find_all("td")
        # add links to survey column
        link_tag = soup.new_tag(
            "a",
            href=pulsarsurveyscraper.Surveys[cols[0].text]["url"],
        )
        link_tag.string = cols[0].text
        cols[0].string = ""
        cols[0].insert(0, link_tag)
        for i in col_aligns:
            cols[i]["align"] = col_aligns[i]
    html_table = soup
    return render_template("surveys.html", html_table=html_table)


# run flask app
if __name__ == "__main__":
    app.run(debug=True)
