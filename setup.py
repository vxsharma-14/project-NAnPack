from setuptools import setup

setup(
name='NAnPack-Learners-Edition',
version='1.0.0.dev1',
author='Vishal Sharma',
author_email='sharma_vishal14@hotmail.com',
url='https://github.com/vxsharma-14/NAnPack',
license='LICENSE.md',
description='A package for simulating physical processes.',
long_description=open('README.md').read(),
long_description_content_type='text/markdown',
package_dir = {'':'Project_NAnPack'},
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
