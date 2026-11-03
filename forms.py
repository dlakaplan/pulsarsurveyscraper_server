from flask_wtf import FlaskForm
from wtforms import StringField, DecimalField, SubmitField, RadioField, BooleanField
from wtforms.validators import InputRequired, ValidationError, Optional
import re
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

import pulsarsurveyscraper

radec_label = "&#x3B1 &#x3B4"
lb_label = "<i>&#x2113</i> <i>b</i>"


def parse_equcoord_and_validate(form, data):
    """
    parse equatorial input form coordinates and validate them

    the input is a string (RA Dec pair in various forms that astropy can understand)
    if it can parse, returns astropy SkyCoord object
    if it cannot, raises ValidationError
    """

    # divide the string into equal RA, Dec pieces
    if "," in data.data:
        ra, dec = map(str.strip, data.data.split(","))
    else:
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
    if "," in data.data:
        gal_l, gal_b = map(str.strip, data.data.split(","))
    else:
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
            "Unable to parse input coordinate = '{},{}': {}".format(gal_l, gal_b, e)
        )
    return coord

def parse_sidelobe_arc_and_validate(form, data):
    """
    record whether to search along a CHIME-like sidelobe arc or in a cone

    the input is a boolean (True for arc, False for cone)
    sets form.parsed_cone_or_arc to "arc" or "cone" and returns the boolean value
    """

    form.parsed_cone_or_arc = "arc" if data.data else "cone"

    return data.data


def parse_arc_telescope(data):
    """
    parse arc telescope input

    the input is a string (either a site name or lon lat pair)
    if it can parse, returns astropy EarthLocation object
    if it cannot, raises ValidationError
    """
    
    def is_numeric(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    parts = data.split()

    if len(parts) == 2 and is_numeric(parts[0]) and is_numeric(parts[1]):
        try:
            telescope = EarthLocation(lon=parts[0], lat=parts[1])
        except Exception as _:
            raise ValidationError(f"Unable to parse telescope coordinates.")
    else:
        try:
            telescope = EarthLocation.of_site(parts[0].strip())
        except Exception as _:
            raise ValidationError(f"Unable to parse the input as a telescope site, or the telescope name is unknown.")
    
    return telescope


def parse_arc_telescope_and_validate(form, data):
    """
    parse arc telescope input and validate it

    the input is a string (either a site name or lon lat pair)
    if it can parse, returns astropy EarthLocation object
    if it cannot, raises ValidationError
    """

    # Check if the arc mode is on
    if form is not None and form.parsed_cone_or_arc != "arc":
        return

    # Check if the input is empty
    if not data.data:
        raise ValidationError("Telescope input is required for arc mode.")

    # Parse the telescope input using the helper function
    telescope = parse_arc_telescope(data.data)

    # Set the parsed telescope on the form for later use
    if form is not None:
        form.parsed_arc_telescope = telescope
    
    return telescope


def parse_arc_epoch(data):
    """
    parse arc epoch input

    the input is a string representing an arc epoch
    if it can parse, returns astropy Time object
    if it cannot, raises ValidationError
    """

    try:
        epoch = pulsarsurveyscraper.parse_epoch(data)
    except Exception as _:
        raise ValidationError(
            "Unable to parse arc epoch = '{}'".format(data)
        )

    return epoch

def parse_arc_epoch_and_validate(form, data):
    """
    parse arc epoch input and validate it

    the input is a string representing an arc epoch
    if it can parse, returns astropy Time object
    if it cannot, raises ValidationError
    """

    # Check if the arc mode is on
    if form is not None and form.parsed_cone_or_arc != "arc":
        return

    # Check if the input is empty
    if not data.data:
        raise ValidationError("Epoch input is required for arc mode.")

    # Parse the epoch input using the helper function
    epoch = parse_arc_epoch(data.data)

    # Set the parsed epoch on the form for later use
    if form is not None:
        form.parsed_arc_epoch = epoch

    return epoch


class SearchForm(FlaskForm):
    """
    basic SearchForm

    contents:
        coordinates (string): gets parsed and validated using parse_coord_and_validate
        lb_or_radec (bool): whether coordinates are equatorial or galactic
        radius (decimal)
        dm (decimal): optional
        dmtol (decimal): optional
        search (button): for executing main search
        api (button): for submitting as API
        clear (button): for clearing the form
        PDF (button): for downloading a PDF
        PNG (button): for downloading a PNG
    """

    coordinates = StringField(
        "Search Coordinates", validators=[InputRequired(), parse_equcoord_and_validate],
    )
    lb_or_radec = BooleanField(
        "Equatorial or Galactic", validators=[], id="coord-toggle", default=True,
    )
    sidelobe_arc = BooleanField(
        "Search along a CHIME-like sidelobe arc instead of a cone",
        validators=[parse_sidelobe_arc_and_validate],
        id="search-sidelobe-arc",
        default=False,
    )
    deduplicate = BooleanField("Deduplicate", validators=[], default=False)

    arc_telescope = StringField(
        "Telescope Name or Coordinates", validators=[parse_arc_telescope_and_validate], default="CHIME"
    )
    arc_epoch = StringField(
        "Observation Date", validators=[parse_arc_epoch_and_validate]
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
    PDF = SubmitField(label="PDF")
    PNG = SubmitField(label="PNG")


class DMForm(FlaskForm):
    """
    basic DMForm

    contents:
        coordinates (string): gets parsed and validated using parse_equcoord_and_validate
        lb_or_radec_selector (radio): is the input RA,Dec or l,b?
        d_or_dm (decimal): input value of distance (pc) or DM
        d_or_dm_selector (radio): is the input distance or DM?
        frequency (decimal): frequency for scintillation calculations
        velocity (decimal): velocity for scintillation calculations
        model_selector (radio): is the model NE2001 or YMW16?
        compute (button): for executing main search
        api (button): for submitting as API
        clear (button): for clearing the form
    """

    coordinates = StringField(
        radec_label, validators=[InputRequired(), parse_equcoord_and_validate],
    )
    d_or_dm = DecimalField("Distance (pc)", validators=[InputRequired()])
    d_or_dm_selector = BooleanField("Distance or DM", id="input-toggle", validators=[],)
    lb_or_radec = BooleanField(
        "Equatorial or Galactic", validators=[], default=True, id="coord-toggle",
    )
    frequency = DecimalField(
        "Frequency (GHz)", default=1, validators=[Optional()]
    )
    velocity = DecimalField(
        "Velocity (km/s)", default=100, validators=[Optional()]
    )
    model_selector = RadioField(
        "DM Model",
        default="ne2001",
        choices=[("ne2001", "NE2001"), ("ymw16", "YMW16")],
        validators=[InputRequired()],
    )
    model_selector = BooleanField("Model", validators=[], id="model-toggle")
    compute = SubmitField(label="Compute")
    api = SubmitField(label="API")
    clear = SubmitField(label="Clear")