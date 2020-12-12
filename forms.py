from flask_wtf import FlaskForm
from wtforms import StringField, DecimalField, SubmitField
from wtforms.validators import InputRequired, ValidationError, Optional
import re
from astropy.coordinates import SkyCoord
from astropy import units as u

import pulsarsurveyscraper


def parse_coord_and_validate(form, data):
    """
    parse input form coordinates and validate them

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
    except ValueError:
        try:
            # if that hasn't worked, try to assume hours/degrees
            coord = SkyCoord(ra, dec, unit=("hour", "deg"))
        except ValueError as e:
            raise ValidationError(
                "Unable to parse input RA,Dec = '{},{}': {}".format(ra, dec, e)
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
        validators=[InputRequired(), parse_coord_and_validate],
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
