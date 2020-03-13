# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 12:47:52 2020

@author: samur
"""


import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import scipy.io
import os.path
from numpy import linalg as LA
import pandas as pd

from os import listdir
from os.path import isfile, join
import re
import csv
from copy import deepcopy

yeast8=pd.read_excel('Yeast8.xlsx')
#list(yeast8.columns)
highflux=pd.read_table('flux_profile_yeast8_highglucose.txt',index_col=None,header=None)
lowflux=pd.read_table('flux_profile_yeast8_lowglucose.txt',index_col=None,header=None)

yeast8['lowglucoseflux']=lowflux
yeast8['highglucoseflux']=highflux

yeast8.to_csv('yeast8_glucose_fluxes.csv')

