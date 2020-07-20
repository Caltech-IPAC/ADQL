from setuptools import setup

setup(
    name='ADQL',
    version='1.0.6',
    author='John Good',
    author_email='jcg@ipac.caltech.edu',
    description='Utility for converting ADQL to local SQL syntax.',
    long_description=open('README.md').read(),
    license='LICENSE',
    keywords='astronomy database ADQL SQL',
    url='https://github.com/Caltech-IPAC/ADQL',
    install_requires=['spatial_index'],
    packages=['ADQL']
)
