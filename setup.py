#!/usr/bin/env python

from setuptools import setup, find_packages
import glob
import shutil
import os

script = 'peak/peak'
shutil.copy(f'{script}.py', script)
os.chmod(script, 0o755)

scripts = glob.glob('peak/*.pl')
scripts.append(script)

with open('README.md') as f:
    long_description = f.read()

setup(
    name='eclip-peak',
    version='1.0.6',
    description='Pipeline for using IDR to identify a set of reproducible peaks given eClIP dataset with '
                'two or three replicates.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/VanNostrandLab/peak',
    author='FEI YUAN',
    author_email='fei.yuan@bcm.edu',
    license='MIT',
    license_files='LICENSE',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX'
    ],
    keywords='eCLIP-seq, peaks, bioinformatics',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'cmder>=1.0',
        'pandas>=1.2.3',
        'inflect>=5.3.0',
        'seqflow>=0.0.3',
    ],
    scripts=scripts
)

