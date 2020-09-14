#!/usr/bin/python3

from setuptools import setup, find_packages

setup(name='aurauas_flightdata',
      version='1.2',
      description='Flight data management libs',
      author='Curtis L. Olson',
      author_email='curtolson@flightgear.org',
      url='https://github.com/AuraUAS',
      #package_dir = {'': 'aurauas'},
      #packages = ['aurauas_flightdata']
      packages = find_packages()
     )
