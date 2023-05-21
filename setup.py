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
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Operating System :: OS Independent'
]


setup(
	name='ndvi2gif',
	version='0.0.4',
	description='Python package to create seasonal composites of NDVI or several more vegetation, water, snow, etc. indexes, and download them as pretty gifs or geotiff rasters',
	long_description=read('README.rst'),
	url='https://github.com/Digdgeo/Ndvi2Gif',
	python_requires='>=3.8',
	author='Diego Garcia Diaz',
	author_email='digd.geografo@gmail.com',
	license='MIT',
	install_requires=['geemap >= 0.19.0', 'deims >= 3.1'],
	packages=find_packages(include=['ndvi2gif', 'ndvi2gif.*']),
	zip_safe=False
)
