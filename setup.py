from setuptools import setup
from setuptools import find_packages

setup(
name='NAnPack',
version='1.0.0-alpha4',
author='Dr. Vishal Sharma',
author_email='sharma_vishal14@hotmail.com',
url='https://github.com/vxsharma-14/project-NAnPack',
license='LICENSE.md',
description='A package of scientific computing tools for learning and teaching.',
long_description=open('README.md').read(),
long_description_content_type='text/markdown',
include_package_data=True,
packages=['nanpack'],
#packages=find_packages(),
install_requires=['numpy', 'matplotlib'],
python_requires='>=3.7',
classifiers=[
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 3.7',
    'License :: OSI Approved :: MIT License',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
],
keywords='numerical methods ' 'computational engineering '\
'computational fluid dynamics ' 'scientific computing',
)
