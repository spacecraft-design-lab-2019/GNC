# -*- coding: utf-8 -*-
"""
Script for testing pyIGRF
Created on Tue Oct  1 22:03:04 2019

@author: Paul DeTrempe
"""

import pyIGRF

lat = 0.0       # [-90,90] degrees
lon = 0.0       # east longitude, [0,360] degrees east of Greenwich 
alt = 400.0     # height, km
date = 2019.0   # Year


(D,I,H,Bx,By,Bz,F) = pyIGRF.igrf_variation(lat, lon, alt, date)
# Declination, Inclination, Horizontal intensity, North, East, Down, total intensity (nT)