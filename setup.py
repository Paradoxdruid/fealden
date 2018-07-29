from setuptools import setup

setup(
    name='fealden',
    url='https://github.com/Paradoxdruid/fealden',
    author='Andrew J. Bonham',
    author_email='bonham@gmail.com',
    packages=['fealden'],
    dependency_links=['https://rna.urmc.rochester.edu/RNAstructure.html'],
    version='0.1',
    license='LGPL v3.0',
    scripts=['bin/fealden'],
    include_package_data=True,
    description='Tool for generating optimized nucleic acid biosensor sequences',
    long_description=open('README.md').read(),
)
