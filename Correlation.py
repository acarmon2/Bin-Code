#!/usr/bin/env python
# -*- coding: utf-8 -*-

#########################################################################################################
# functions to calculate the cross correlation in time domain, in Spectral.csv has the power response
# of correlational parameters. 
# pandas is the library ro read .csv data (comma separate dataset, not point and comma)
# matplotlib contain the basis to plot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

Spectral0 = pd.read_csv("Spectral.csv")
# convert dataframe to numpy arrays
Specter = Spectral0.values
