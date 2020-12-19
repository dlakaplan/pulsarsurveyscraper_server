import flask
from flask import Flask, render_template, request, url_for, redirect
from flask_wtf import FlaskForm
from wtforms import StringField, DecimalField, SubmitField, RadioField
from wtforms.validators import InputRequired, ValidationError, Optional

# create the Flask object
app = Flask(__name__)

app.config.from_object("server_config")

coordinate_type = 'equatorial'

class DMForm(FlaskForm):
    coordinates = StringField(
        "Search Coordinate (RA Dec or l b)",
        validators=[InputRequired()],
    )

@app.route("/", methods=["GET", "POST"])
def Something():
    form=DMForm()
    if form.validate_on_submit():
        print('Coordinate type: {}'.format(coordinate_type))

    return render_template("toggle.html", form=form)

@app.route('/get_coordtoggled_status') 
def toggled_status():
  current_status = flask.request.args.get('status')
  print(current_status)
  coordinate_type = 'equatorial' if True else 'galactic'
  return coordinate_type


# run flask app
if __name__ == "__main__":
    app.run(debug=True)
