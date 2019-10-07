#!/usr/bin/env python
# -*- coding: utf-8 -*-

###################################################################################################
# main function to calculate the generators
import numpy as np
import time
from Bin import *

start = time.time()

X0 = BinSearch(29, 8, 1000)
print np.shape(X0)

end = time.time()
print (end - start)