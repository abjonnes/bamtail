# -*- coding: utf-8 -*-

from setuptools import setup

from bamtail import __version__

setup(
    name='bamtail',
    version=__version__,
    description='tail for BAM files',
    author='Andrew Bjonnes',
    author_email='andrew@abjonnes.com',
    url='https://github.com/abjonnes/bamtail',
    py_modules=['bamtail'],
    entry_points={
        'console_scripts': ['bamtail = bamtail:main']
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
