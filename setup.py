#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = ['numpy>=1.9.1', 'CosmoloPy>=0.1.105', 'scipy>=0.15.1']

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='commah',
    version='0.1.0',
    description='Find concentration - mass relation MAH for any cosmology',
    long_description=readme + '\n\n' + history,
    author='Alan Duffy',
    author_email='mail@alanrduffy.com',
    url='https://github.com/astroduff/commah',
    packages=[
        'commah'
    ],
    package_dir={'commah':
                 'commah'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='commah',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements
)