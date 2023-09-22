from setuptools import setup

setup(
    name='ontosunburst',
    packages=['ontosunburst'],
    author="Pauline Hamon-Giraud",
    author_email="pauline.hamon-giraud@irisa.fr",
    description="Sunbursts for ontologies",
    install_requires=['dash', 'numpy', 'padmet', 'plotly', 'scipy', 'SPARQLWrapper'],
    url='https://github.com/PaulineGHG/Ontology_sunburst.git',
)