#
#  NordicGrid.py
#  Power flow within the Nordic region with Denmark East and West at the center.
#
#  Created by Gorm Bruun Andresen on 08/12/2011.
#  Copyright (c) 2011 Department of Engineering, University of Aarhus. All rights reserved.
#

#Standard modules
from pylab import *
from scipy import *
import os, sys

#Custom modules
sys.path.append( './zdcpf/' ) #This can be done in a more fancy way using __init__.py or some such.
from zdcpf import *


N = Nodes()