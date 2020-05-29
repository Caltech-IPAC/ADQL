import os

from setuptools import setup

setup(
    name = "ADQL",
    version = '1.0',
    author = 'John Good',
    author_email = 'jcg@ipac.caltech.edu',
    description = 'Utility for converting ADQL to local SQL syntax.',
    long_description = open('README.txt').read(),
    license = 'LICENSE',
    keywords = 'astronomy database ADQL SQL',
    url = 'https://github.com/Caltech-IPAC/pyADQL',
    install_requires=['spatial_index']
)
