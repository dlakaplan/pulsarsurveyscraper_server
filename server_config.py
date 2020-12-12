import os

DEBUG = True

# python -c 'import os; print(os.urandom(16))'
SECRET_KEY = os.environ.get("SECRET_KEY")
if not SECRET_KEY:
    raise ValueError("No SECRET_KEY set for Flask application")

# export SURVEY_DATA=/path/to/survey/data
DATA_DIR = os.environ.get("SURVEY_DATA")
if not DATA_DIR:
    DATA_DIR = os.path.curdir
