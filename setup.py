from setuptools import setup

setup(
name="NAnPack-Learner's-Edition",
version='1.0.0.dev1',
author='Vishal Sharma',
author_email='sharma_vishal14@hotmail.com',
#package=['NAnPack', 'NAnPack.test'],
scripts=['bin/script1','bin/script2'],
url='http://pypi.python.org/pypi/PackageName/',
license='LICENSE.md',
description='A package for simulating physical processes.',
long_description=open('README.md').read(),
long_description_content_type='text/markdown',
packages=setuptools.find_packages(),
python_requires='>=3.7',
classifiers=[
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 3.7',
    'License :: OSI Approved :: MIT License',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
],
keywords='numerical methods ' 'computational engineering '\
'computational fluid dynamics ' 'numerical simulation',
)
