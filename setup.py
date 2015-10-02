#!/usr/bin/env python2

from setuptools import setup, find_packages
import simrrls

requirements = [
    'numpy>1.7',
    ]

## loads __version__
#exec(open('simrrls/version.py').read())

setup(
    name="simrrls",
    version=simrrls.__version__,
    url="https://github.com/dereneaton/simrrls",

    author=simrrls.__author__,
    author_email=simrrls.__contact__,

    description="simulating RADseq data sets",
    long_description=open('README.rst').read(),

    packages=find_packages(),
    
    install_requires=[requirements],

    entry_points={
            'console_scripts': [
                'simrrls = simrrls.__main__:main',
            ],
    },

    license='GPL',

    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
