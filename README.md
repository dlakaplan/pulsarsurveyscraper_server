# pulsarsurveyscraper_server

## requirements:
* flask, flask_wtf, astropy, pandas, scipy
* [PyGEDM](https://github.com/telegraphic/pygedm): for DM models
* [pulsarsurveyscraper](https://github.com/dlakaplan/pulsarsurveyscraper)

## configuration:
Read from environment variables in `server_config.py`:
* `SECRET_KEY`: created using something like `python -c 'import os; print(os.urandom(16))'`
* `SURVEY_DATA`: location of survey cache files (default is current directory)

## deployment:
Based on (this tutorial)[https://www.digitalocean.com/community/tutorials/how-to-serve-flask-applications-with-gunicorn-and-nginx-on-ubuntu-18-04]
