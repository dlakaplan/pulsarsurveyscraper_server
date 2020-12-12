# pulsarsurveyscraper_server

## requirements:
* flask, flask_wtf, astropy, pandas, 
* [pulsarsurveyscraper](https://github.com/dlakaplan/pulsarsurveyscraper)

## configuration:
Read from environment variables in `server_config.py`:
* `SECRET_KEY`: created using something like `python -c 'import os; print(os.urandom(16))'`
* `SURVEY_DATA`: location of survey cache files (default is current directory)
