#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages


requirements = [ ]

test_requirements = [ ]

setup(
    author="Wang Jiaxuan",
    author_email='poormouse@126.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="the wonderful cell type auto-annotation software",
    entry_points={
        'console_scripts': [
            'cellname=cellname.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    include_package_data=True,
    keywords='cellname',
    name='cellname',
    packages=find_packages(include=['cellname', 'cellname.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/wangjiaxuan666/cellname',
    version='0.1.0',
    zip_safe=False,
    package_data={'cellname': [' data/s_scanpy.h5ad ']},
)
