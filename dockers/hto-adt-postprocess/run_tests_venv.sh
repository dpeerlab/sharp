## LOCAL
# using pyenv
# python 3.12 (!)

# install python 12
# pyenv install 3.12.5
virtualenv sharp_test

# install requirements
# pyenv activate sharp_test
source sharp_test/bin/activate
pip install -r requirements.txt

# run tests
pytest -vvv --local

# run ruff for linting
pip install ruff
pip install black

ruff check --fix
black *

# uninstall environment
# pyenv uninstall sharp_test
rm -rf sharp_test