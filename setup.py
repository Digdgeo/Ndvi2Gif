import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

classifiers = [
	'Intended Audience :: GIS & Remote Sensing Technicians',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Operating System :: OS Independent'
]


setup(
	name='ndvi2gif',
	version='0.0.3',
	description='Python package to create ndvi seasonal composites, and download them as gif and geotiff',
	long_description=read('README.rst'),
	url='https://github.com/Digdgeo/Ndvi2Gif',
	python_requires='>=3.5',
	author='Diego Garcia Diaz',
	author_email='digd.geografo@gmail.com',
	license='MIT',
	install_requires=['geemap >= 0.6.10'],
	packages=find_packages(include=['ndvi2gif', 'ndvi2gif.*']),
	zip_safe=False
)