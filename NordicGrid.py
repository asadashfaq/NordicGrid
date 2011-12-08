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

#Special functions
from copy import deepcopy

#Custom modules
sys.path.append( './zdcpf/' ) #This can be done in a more fancy way using __init__.py or some such.
from zdcpf import *


#This function will later be replaced by some fancy save/load thing.
def get_nodes_and_flows_vs_year(year=linspace(1985,2053,21),incidence='incidence.txt',constraints='constraints.txt',setupfile='setupnodes.txt',coop=0,copper=0,path='./settings/',lapse=None):

    #Initialize
    N = Nodes()
    Gamma = get_basepath_gamma(year)

    data_nodes, data_flows = [], []
    for i in arange(len(year)):
        print 'Year: {0:.0f}'.format(year[i])
        sys.stdout.flush()

        #Apply year to nodes (Dummy gamma's)
        N.set_gammas(Gamma.transpose()[i])
        N.set_alphas([1.,1.,1.,1.,.9])
        
        print 'Gamma: ' + str(Gamma.transpose()[i])
        sys.stdout.flush()
        
        #Calculate flows
        N,F = zdcpf(N,incidence=incidence,constraints=constraints,setupfile=setupfile,path=path,coop=coop,copper=copper,lapse=lapse)
        
        #append data
        data_nodes.append(deepcopy(N))
        data_flows.append(deepcopy(F))

    return year, data_nodes, data_flows

def plot_generation_summary_vs_year(year=linspace(1985,2053,21),node_id=3,lapse=50*24,data=None):

    return 1
    
### Local utilities
def get_basepath_gamma(year,filename='basepath_gamma.npy'):
    """ File should be moved to subdir.
    """
    print "Loading: {0}. Warning columns not pre-labeled!!".format(filename)
    data = np.load(filename)

    Gamma = zeros((len(data)-1,len(year)))
    for i in arange(len(Gamma)):
        Gamma[i] = interp(year,data[0],data[i+1])

    return Gamma