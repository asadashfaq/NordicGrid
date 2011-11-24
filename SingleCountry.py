#
#  SingleCountry.py
#  Maximum local integration wind and solar power generation for a region/country.
#
#  Created by Gorm Bruun Andresen on 24/11/2011.
#  Copyright (c) 2011 Department of Engineering, University of Aarhus. All rights reserved.
#

#Standard modules
from pylab import *
from scipy import *
import os

#Custom functions
from Database_v1 import get_data_countries, get_data_regions

def generate_data_files():
    """Works from anywhere  if you setup an ssh tunnel first: ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk"""
    # t, L, GW, GS, datetime_offset, datalabels = get_data_countries(localhost=True);
    # t, l, Gw, Gs, datetime_offset, datalabels = get_data_countries(schema='norm_agg_avg_1hour_pdata_caps_eu2020',localhost=True);