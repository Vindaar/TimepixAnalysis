from setuptools import setup

setup(
    name='InGrid',    # This is the name of your PyPI-package.
    version='0.1',                          # Update the version number for new releases
    packages=['ingrid'],
    scripts=['ingrid/ingrid_helper_functions.py']     # The name of your scipt, and also the command you'll be using for calling it
)
