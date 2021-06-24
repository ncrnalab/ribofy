"""
A setuptools based setup module.
See:

https://packaging.python.org/en/latest/distributing.html

https://github.com/pypa/sampleproject

"""


# Always prefer setuptools over distutils

from setuptools import setup, find_packages
from os import path
import sys
from ribofy import __version__

if sys.version_info[0] == 2:
    sys.stderr.write("Ribofy only supports Python3")
    sys.exit(1)

here = path.abspath(path.dirname(__file__))


setup(

    name='ribofy',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html

    

    description='Ribofy: ORF detection using RiboSeq data',

    long_description='Ribofy is a fast, simple and conservative pipeline for ORFs detection using RiboSeq data.',

    # The project's main homepage.
    url='https://github.com/ncrnalab/ribofy',

    # Author details
    author='Thomas Hansen',
    author_email='tbh@mbg.au.dk',

    # Version
    version = __version__,

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are

        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable

        'Development Status :: 4 - Beta',
        #'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',

    ],
    project_urls ={
        'Documentation': 'https://github.com/ncrnalab/ribofy/blob/master/README.rst',
        'Source': 'https://github.com/ncrnalab/ribofy',
        'ChangeLog': 'https://github.com/ncrnalab/ribofy/blob/master/ChangeLog.rst',
        'Issues': 'https://github.com/ncrnalab/ribofy/issues',
    },

    # What does your project relate to?    
    keywords='ribo-seq ribosome-profiling ORF',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().

    packages=['ribofy'], #find_packages(exclude=['contrib', 'docs', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:

    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's

    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html

    install_requires=['pysam>0.8.4','matplotlib','numpy','scipy', 'networkx', 'pandas', 'tqdm', 'argparse'],

    #
    package_data = {},

    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa

    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow

    # pip to create the appropriate form of executable for the target platform.

    entry_points={

        'console_scripts': [
            #name determined the name of cmd line direct call
            'ribofy=ribofy.ribofy:main'
        ],

    },

)