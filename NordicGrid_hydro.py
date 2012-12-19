#
#  NordicGrid_hydro.py
#  Power flow within the Nordic region with Denmark East and West at the center. The influence of storage lakes in Norway and Sweeden has been added.
#
#  Created by Gorm Bruun Andresen on 28/08/2012.
#  Copyright (c) 2012 Department of Engineering, University of Aarhus. All rights reserved.
#

## Standard modules
from pylab import *
from scipy import *
import os, sys

## Custom modules and functions
sys.path.append( './zdcpf/' ) #This can be done in a more fancy way using __init__.py or some such.
sys.path.append( './EuropeanGridR/' )
from storage_classes import one_way_storage
from water import get_inflow, get_median_storage_level_Norway, get_median_storage_power_Norway
from shortcuts import *
from SingleCountry import get_ISET_country_data

#Specific functions
from scipy.optimize import brentq


## Standard figure size:
figure_size = [6.5,4.3]


def plot_inflow(t=arange(365*24)):

    inflow, inflow_basic = get_inflow(t,returnall=True)

    ### Plot inflow
    close(1);figure(1);clf()

    plot(inflow)
    plot(inflow_basic)
    
    tight_layout()
    save_file_name = 'plot_inflow.pdf'
    save_figure(save_file_name)
    
def plot_one_way_storage(t=None,volume=80e3,P_out=30,N=None):

    storage_lake = one_way_storage(volume=volume,P_out=P_out,inflow=get_inflow(t))
    if N==None:
        N = len(storage_lake.inflow)
    
    
    ### Plot storage time series for a one-way storage based on the Norwegian storage lakes.
    close(1);figure(1);clf()

    subplot(211)
    plot(storage_lake.inflow,label='Inflow')
    plot(storage_lake.default_output,label='Default output')
    plot(storage_lake.forced_overflow,'r',label='Forced overflow')
    
    
    xlabel('Time [h]')
    ylabel('Power [GWh/h]')
    
    legend()
    
    axis(xmin=0,xmax=N,ymin=0, ymax=75)
    
    subplot(212)
    fill_between(arange(len(storage_lake.virtual_two_way_storage.prefilling)),storage_lake.virtual_two_way_storage.prefilling/volume,edgecolor=None,lw=0,facecolor='r',alpha=0.5,label='Prefilling')
    
    xlabel('Time [h]')
    ylabel('Storage filling fraction')
    
    axis(xmin=0,xmax=N,ymin=0,ymax=.25)
    
    tight_layout()
    save_file_name = 'plot_one_way_storage.png'
    save_figure(save_file_name)
    
    return storage_lake

def plot_isolated_Norway(ISO='NO',volume=80e3,P_out=30,N_days=4*365,gamma=0,alpha_w=1):
    
    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    #L = L*1.1
    
    inflow, storage_level, t, datetime_offset = get_inflow(t,datetime_offset, returnall=True)
    
    if N_days==None:
        N_days = amax(t)
    
    #median_level = volume*get_median_storage_level_Norway()[0]
    median_level = volume*get_median_storage_power_Norway(returnall=True)[2][0]
    median_level = kron(ones(ceil(len(t)/365./24.)),median_level)
    median_level = copy(median_level[0:len(t)])
    
    storage_lake = one_way_storage(volume=volume,P_out=P_out,inflow=inflow,median_level=median_level)
    Gh = storage_lake.default_output  # Hydro default output
    
    mismatch = gamma*mean(L)*(alpha_w*Gw + (1.-alpha_w)*Gs) + Gh - L
    mismatch_ = mismatch
    #print 'Average hourly deficit: ' + str(sum(mismatch)/len(mismatch))
    deficit = sum(mismatch)
    print 'Mismatch sum (before balancing): ' + str(sum(mismatch))
    #balancing = lambda x: get_positive(-(mismatch_+x))
    
    #x_bal = brentq(lambda x: sum(balancing(x))+deficit, amin(mismatch), amax(mismatch))
    #print x_bal# - 1e-8 
    
    #mismatch = mismatch + balancing(x_bal)
    balancing = deficit/len(t)*ones_like(t)
    mismatch = mismatch - balancing
    print 'Mismatch sum (after balancing): ' + str(sum(mismatch))
    
    mismatch_r = zeros_like(mismatch)
    
    i, ii = 0, 0
    storage_lake.virtual_two_way_storage.level[-1] = storage_level[0]
    print storage_lake.virtual_two_way_storage.level[-1] - median_level[0]
    storage_lake.virtual_two_way_storage.power_balance[-1] = storage_lake.virtual_two_way_storage.level[-1] - median_level[0]
    while storage_lake.iscyclic()!=True:
    
        if mismatch[i]<0:
            charge = amax([-storage_lake.get_power_out(),mismatch[i]])
        else:
            charge = amin([storage_lake.get_power_in(),mismatch[i]])
    
        storage_lake.charge(charge)
        mismatch_r[i] = mismatch[i] - charge
        
        if mod(i+1,len(mismatch))==0:
            #storage_lake.virtual_two_way_storage.level[-1] = storage_level[0] # Testing reset of initial level.
            ii = ii +1
            print ii, i
            if ii==1:
                break
        i = mod(i+1,len(mismatch))
    print ii, i
    print storage_lake.virtual_two_way_storage.level[:3], storage_lake.virtual_two_way_storage.level[-3:]
    
    close(1);figure(1);clf()
    
    subplot(211)
    plot(t,-storage_lake.virtual_two_way_storage.ekstra_discharge,'y-',label='Ekstra discharge')
    plot(t,mismatch,label='Mismatch')
    plot(t,mismatch_r,label='Mismatch_r')
    fill_between(t,balancing,edgecolor=None,lw=0,facecolor='r')
    
    plot(t,storage_lake.virtual_two_way_storage.power_cap_in,'k-',lw=.5,alpha=.5)
    plot(t,-storage_lake.virtual_two_way_storage.power_cap_out,'k-',lw=.5,alpha=.5)
    
    plot(t,storage_lake.virtual_two_way_storage.default_charge,'k-',lw=.5,label='Default charge')
    plot(t,storage_lake.virtual_two_way_storage.default_charge-volume*0.5/365./24,'k--',lw=.5,label='Default charge (offset)')
    plot(t,storage_lake.virtual_two_way_storage.power_balance/volume*100,'r-',lw=.5,label='Power balance')
    
    
    xlabel('Time [h]')
    ylabel('Power [GWh/h]')
    
    legend()
    
    axis(xmin=0,xmax=N_days,ymin=-30,ymax=30)
    
    subplot(212)
    
    fill_between(t,storage_lake.get_level()/volume,edgecolor=None,lw=0,facecolor='r',alpha=0.5,label='Prefilling')
    
    ## Plot real storage levels
    #storage_level_real_=kron(np.genfromtxt('./raw_data/meannor'),ones(7*24))
    #storage_level_real = zeros_like(t)
    #for i in arange(len(storage_level_real)):
    #    storage_level_real[i] = storage_level_real_[mod(i,len(storage_level_real_))]
    plot(t,storage_level/volume,'k-',label='Historical level')
    plot(t,median_level/volume,'k--',label='Median level')
    
    min_level = zeros_like(median_level)
    for i in arange(len(median_level)):
        min_level[i] = -amin(median_level - median_level[i])
    plot(t,min_level/volume,'r-',label='Min level')
            
    xlabel('Time [h]')
    ylabel('Storage filling fraction')
    
    legend()
    
    axis(xmin=0,xmax=N_days,ymin=0,ymax=1.)
    
    tight_layout()
    save_file_name = 'plot_isolated_Norway.png'
    save_figure(save_file_name)    
    
    