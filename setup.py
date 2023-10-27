'''
setup script for mcflow
'''

from setuptools import setup, find_packages

# Read the requirements.txt file
with open('requirements.txt', 'r',encoding='utf-8') as f:
    requirements = f.read().splitlines()

setup(
    name="mcflow",
    version="0.1.0",
    description="radiological modelling of radioisotopes contained in activated water",
    author="Marco De Pietri",
    author_email="mdepietri@ind.uned.es",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[requirements],
    package_data={
        'mcflow': ['*.total', ],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mcflow=mcflow.mcflow:main',
        ],
    }
)
