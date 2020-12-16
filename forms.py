from flask_wtf import FlaskForm
from wtforms import StringField, DecimalField, SubmitField, RadioField
from wtforms.validators import InputRequired, ValidationError, Optional
import re
from astropy.coordinates import SkyCoord
from astropy import units as u

import pulsarsurveyscraper


def parse_equcoord_and_validate(form, data):
    """
    parse equatorial input form coordinates and validate them

    the input is a string (RA Dec pair in various forms that astropy can understand)
    if it can parse, returns astropy SkyCoord object
    if it cannot, raises ValidationError
    """

    # divide the string into equal RA, Dec pieces
    c = data.data.split()
    l = len(c)
    ra = " ".join(c[: (l // 2)])
    dec = " ".join(c[(l // 2) :])

    try:
        if (re.search(r"[^\d.+\-]", ra) is None) and (
            re.search(r"[^\d.+\-]", dec) is None
        ):
            # if it only has digits/sign/decimal, it's probably decimal degrees
            coord = SkyCoord(ra, dec, unit="deg")
        else:
            # see if astropy can figure out the units
            coord = SkyCoord(ra, dec)
    except (ValueError, u.core.UnitsError):
        try:
            # if that hasn't worked, try to assume hours/degrees
            coord = SkyCoord(ra, dec, unit=("hour", "deg"))
        except ValueError as e:
            raise ValidationError(
                "Unable to parse input coordinate = '{},{}': {}".format(ra, dec, e)
            )
    except u.core.UnitsError as e:
        raise ValidationError(
            "Unable to parse input coordinate = '{},{}': {}".format(ra, dec, e)
        )
    return coord


def parse_galcoord_and_validate(form, data):
    """
    parse galactic input form coordinates and validate them

    the input is a string (l b pair in various forms that astropy can understand)
    if it can parse, returns astropy SkyCoord object
    if it cannot, raises ValidationError
    """

    # divide the string into equal l, b pieces
    c = data.data.split()
    l = len(c)
    gal_l = " ".join(c[: (l // 2)])
    gal_b = " ".join(c[(l // 2) :])

    try:
        if (re.search(r"[^\d.+\-]", gal_l) is None) and (
            re.search(r"[^\d.+\-]", gal_b) is None
        ):
            # if it only has digits/sign/decimal, it's probably decimal degrees
            coord = SkyCoord(gal_l, gal_b, frame="galactic", unit="deg")
        else:
            # see if astropy can figure out the units
            coord = SkyCoord(gal_l, gal_b, frame="galactic")
    except ValueError as e:
        raise ValidationError(
            "Unable to parse input coordinate = '{},{}': {}".format(gal_l, gal_b, e)
        )
    except u.core.UnitsError as e:
        raise ValidationError(
            "Unable to parse input coordinate = '{},{}': {}".format(ra, dec, e)
        )
    return coord


class SearchForm(FlaskForm):
    """
    basic SearchForm

    contents:
        coordinates (string): gets parsed and validated using parse_coord_and_validate
        radius (decimal)
        dm (decimal): optional
        dmtol (decimal): optional
        search (button): for executing main search
        api (button): for submitting as API
        clear (button): for clearing the form
    """

    coordinates = StringField(
        "Search Coordinate (RA Dec)",
        validators=[InputRequired(), parse_equcoord_and_validate],
    )
    radius = DecimalField(
        "Search Radius (deg)", default=5, validators=[InputRequired()]
    )
    dm = DecimalField("Search DM (pc/cc)", validators=[Optional()])
    dmtol = DecimalField(
        "Search DM Tolerance (pc/cc)", default=10, validators=[Optional()]
    )
    search = SubmitField(label="Search")
    api = SubmitField(label="API")
    clear = SubmitField(label="Clear")


class DMForm(FlaskForm):
    """
    basic DMForm

    contents:
        coordinates (string): gets parsed and validated using parse_equcoord_and_validate
        lb_or_radec_selector (radio): is the input RA,Dec or l,b?
        d_or_dm (decimal): input value of distance (pc) or DM
        d_or_dm_selector (radio): is the input distance or DM?
        model_selector (radio): is the model NE2001 or YMW16?
        compute (button): for executing main search
        clear (button): for clearing the form
    """

    coordinates = StringField(
        "Search Coordinate (RA Dec or l b)",
        validators=[InputRequired(), parse_equcoord_and_validate],
    )
    d_or_dm = DecimalField("Distance (pc) or DM", validators=[InputRequired()])
    d_or_dm_selector = RadioField(
        "Input is Distance or DM",
        default="distance",
        choices=[("distance", "Distance"), ("dm", "DM")],
        validators=[InputRequired()],
    )
    lb_or_radec_selector = RadioField(
        "Equatorial (RA,Dec) or Galactic (l,b)",
        default="equatorial",
        choices=[("equatorial", "Equatorial (RA,Dec)"), ("galactic", "Galactic (l,b)")],
        validators=[InputRequired()],
    )
    model_selector = RadioField(
        "DM Model",
        default="ne2001",
        choices=[("ne2001", "NE2001"), ("ymw16", "YMW16")],
        validators=[InputRequired()],
    )
    compute = SubmitField(label="Compute")
    clear = SubmitField(label="Clear")
