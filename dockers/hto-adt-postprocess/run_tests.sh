## LOCAL
# using pyenv
# python 3.12 (!)

# install python 12
pyenv install 3.12.5
pyenv virtualenv 3.12.5 sharp_test

# install requirements
pyenv activate sharp_test
pip install -r requirements.txt

# run tests
pytest -vvv --local

# run ruff for linting
pip install ruff
pip install black

ruff check --fix
black *

# uninstall environment
pyenv uninstall sharp_test
