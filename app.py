from flask import Flask, render_template, request, url_for, redirect
import flask
from forms import (
    SearchForm,
    DMForm,
    parse_equcoord_and_validate,
    parse_galcoord_and_validate,
)

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas

from bs4 import BeautifulSoup, Tag

import pulsarsurveyscraper
import pygedm

"""
Among others, I used these tutorials:
https://codeloop.org/flask-tutorial-flask-forms-with-flask-wtf/
"""


# create the Flask object
app = Flask(__name__)

app.config.from_object("server_config")

# instantiate the pulsar survey table
# this needs the directory where the data are stored
# do it at the top level so that it's accessible within the functions below
pulsar_table = pulsarsurveyscraper.PulsarTable(directory=app.config["DATA_DIR"],)

coordinate_type = "equatorial"


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
        coord = parse_equcoord_and_validate(None, form.coordinates)
        # figure out if we need to deal with DM as well
        if form.dm.data is not None:
            DM = float(form.dm.data)
        else:
            DM = form.dm.data
        if form.dmtol.data is not None:
            DMtol = float(form.dmtol.data)
        else:
            DMtol = form.dm.data

        # if we've pressed the "API" button
        # instead of showing the tabular output redirect
        # to the API result
        if form.api.data:
            if DM is None:
                return redirect(
                    url_for(
                        "API",
                        ra=coord.ra.deg,
                        dec=coord.dec.deg,
                        radius=float(form.radius.data),
                    )
                )
            else:
                return redirect(
                    url_for(
                        "API",
                        ra=coord.ra.deg,
                        dec=coord.dec.deg,
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
        result = pulsar_table.search(
            coord, radius=float(form.radius.data) * u.deg, DM=DM, DM_tolerance=DMtol
        )
        result["P"][result["P"] < 0] = np.nan
        # make a nice string for output
        coord_string = "Searching {:.1f}deg around RA,Dec {} = {}d,{}d ...".format(
            float(form.radius.data),
            coord.to_string("hmsdms", sep=":"),
            coord.ra.to_string(decimal=True),
            coord.dec.to_string(decimal=True, alwayssign=True),
        )
        if DM is not None:
            coord_string += "<br>Also requiring DM with +/-{:.1f} of {:.1f} pc/cm**2".format(
                DMtol, DM
            )

        # go from astropy Table -> pandas dataframe -> HTML table
        df = result.to_pandas()
        # turn the "PSR" column from bytestring to string
        df["PSR"] = df["PSR"].str.decode("utf-8")
        html_table = df.to_html(
            formatters={"P": lambda x: "%.2f" % x, "Distance": lambda x: "%.2f" % x,},
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
        # fix the alignment of various columns
        col_aligns = {3: "right", 4: "right", 5: "center", 7: "right"}
        rows = soup.find_all("tr")
        for row in rows[1:]:
            cols = row.find_all("td")
            # add links to survey column
            print(cols[5].text, pulsarsurveyscraper.Surveys[cols[5].text]["url"])
            link_tag = soup.new_tag(
                "a", href=pulsarsurveyscraper.Surveys[cols[5].text]["url"],
            )
            link_tag.string = cols[5].text
            cols[5].string = ""
            cols[5].insert(0, link_tag)
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

    ra: float in degrees
    dec: float in degrees
    radius: float in degrees
    dm: float (optional)
    dmtol: float (optional)

    returns JSON
    """

    # default values
    radius = 5 * u.deg
    dm = None
    dmtol = 10

    if "ra" in request.args:
        ra = float(request.args["ra"])
    else:
        return "Error: no RA specified"
    if "dec" in request.args:
        dec = float(request.args["dec"])
    else:
        return "Error: no Dec specified"
    if "radius" in request.args:
        radius = float(request.args["radius"])
    if "dm" in request.args:
        dm = float(request.args["dm"])
    if "dmtol" in request.args:
        dmtol = float(request.args["dmtol"])

    try:
        coord = SkyCoord(ra * u.deg, dec * u.deg)
    except ValueError as e:
        return "Unable to parse RA,Dec = '{},{}': {}".format(ra, dec, e)

    result = pulsar_table.search(
        coord, radius=radius * u.deg, DM=dm, DM_tolerance=dmtol, return_json=True
    )

    return result


@app.route("/compute", methods=["GET", "POST"])
def Compute():
    """
    compute form route

    """
    form = DMForm()

    # if the button has been pressed and the input is valid:
    if form.validate_on_submit():
        print("coordinate = ", coordinate_type)
        print(form.lb_or_radec.data)
        # use the validation routine
        # to parse the coordinates
        if form.lb_or_radec.data:
            coord = parse_equcoord_and_validate(None, form.coordinates)
        else:
            coord = parse_galcoord_and_validate(None, form.coordinates)
        if form.d_or_dm_selector.data == "dm":
            DM = float(form.d_or_dm.data)
            distance, _ = pygedm.dm_to_dist(
                coord.galactic.l, coord.galactic.b, DM, method=form.model_selector.data
            )
        elif form.d_or_dm_selector.data == "distance":
            distance = float(form.d_or_dm.data) * u.pc
            DM, _ = pygedm.dist_to_dm(
                coord.galactic.l,
                coord.galactic.b,
                distance,
                method=form.model_selector.data,
            )

        # or, clear the form if desired
        if form.clear.data:
            return redirect(url_for("Compute"))

        # make a nice string for output
        if form.lb_or_radec.data:
            coord_string = "Computing for RA,Dec {} = {}d,{}d".format(
                coord.icrs.to_string("hmsdms", sep=":"),
                coord.icrs.ra.to_string(decimal=True),
                coord.icrs.dec.to_string(decimal=True, alwayssign=True),
            )
            coord_string += "<br>= l,b = {}d,{}d ...".format(
                coord.galactic.l.to_string(decimal=True),
                coord.galactic.b.to_string(decimal=True),
            )
        else:
            coord_string = "Computing for l,b = {}d,{}d ...".format(
                coord.galactic.l.to_string(decimal=True),
                coord.galactic.b.to_string(decimal=True),
            )
            coord_string += "<br>= RA,Dec {} = {}d,{}d".format(
                coord.icrs.to_string("hmsdms", sep=":"),
                coord.icrs.ra.to_string(decimal=True),
                coord.icrs.dec.to_string(decimal=True, alwayssign=True),
            )
        if form.d_or_dm_selector.data == "dm":
            result_string = "For DM = {:.1f} pc/cc, find distance = {:.1f} pc".format(
                DM, distance.to(u.pc).value
            )
        elif form.d_or_dm_selector.data == "distance":
            for model, model_label in form.model_selector.choices:
                if model == form.model_selector.data:
                    break
            result_string = "For distance = {:.1f} pc, find DM = {:.1f} pc/cc with the {} model".format(
                distance.to(u.pc).value, DM.value, model_label
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


# run flask app
if __name__ == "__main__":
    app.run(debug=True)
