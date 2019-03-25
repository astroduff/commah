#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
from codecs import open

# Package imports
from setuptools import find_packages, setup

# %% SETUP DEFINITION
# Read-in the README.rst file
with open('README.rst', 'r') as f:
    readme = f.read()

# Read-in the HISTORY.rst file
with open('HISTORY.rst', 'r') as f:
    history = f.read().replace(".. :changelog:", "")

# Combine readme and history to the long description
long_description = "\n\n".join([readme, history])

# Get the requirements list
with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()

# Get the test requirements list
with open('requirements_text.txt', 'r') as f:
    test_requirements = f.read().splitlines()

# Get the current package version
version = None
with open('commah/__version__.py', 'r') as f:
    exec(f.read())

# Setup function declaration
setup(
    name='commah',
    version=version,
    description='Get conc - mass or accretion history for any cosmology',
    long_description=long_description,
    author='Alan Duffy',
    author_email='mail@alanrduffy.com',
    url='https://github.com/astroduff/commah',
    download_url='https://github.com/astroduff/commah/tarball/%s' % (version),
    packages=find_packages(),
    package_dir={'commah': 'commah'},
    include_package_data=True,
    install_requires=requirements,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4',
    license="BSD",
    zip_safe=False,
    keywords=('commah', 'cosmology', 'astronomy', 'NFW', 'concentration',
              'accretion'),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
