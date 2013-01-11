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

#Colors. To be placed somewhere else
color_solar = (1.,0.8,0,1)

## Standard figure size:
figure_size = [6.5,4.3]


def plot_inflow_Norway(t=arange(365*24)):

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('NO')
    datetime_offset = datetime_offset
    inflow, storage_level, t, datetime_offset = get_inflow(t,datetime_offset, returnall=True)

    year = array([h.year for h in num2date(t+datetime_offset-1/24.)])
    
    inflow_year = zeros_like(inflow)
    for y in unique(year):
        inflow_year[year==y] = mean(inflow[year==y])

        print 'Inflow ' + str(y)
        print str(mean(inflow[year==y])/mean(inflow)) + '% '

    ### Plot inflow
    close(1);figure(1);clf()

    fill_between(t,inflow_year,edgecolor=None,lw=0,facecolor=(.8,.8,.8),label='Annual average')
    pp_inflow_year = Rectangle((0, 0), 1, 1, facecolor=(.8,.8,.8),lw=0)
    
    plot(t,inflow,'k-',label='Historical values')
    pp_inflow = Line2D([0,1], [0,0], color='k',ls='-')

    axhline(mean(inflow),ls='--',color='k',label='Average')
    pp_inflow_mean = Line2D([0,1], [0,0], color='k',ls='--')
    
    
    ylabel('Hourly inflow to reservoirs [GW]')
    
    t_min = amin(t)
    t_max = amax(t)
    
    axis(xmin=t_min,xmax=t_max,ymin=0)
    
    xticks(t[1+find(diff(year))],year[1+find(diff(year))],rotation=-45,ha='left')
    
    pp = [pp_inflow,pp_inflow_year,pp_inflow_mean]
    pp_text = ['Historical values','Annual average','Average']
    legend(pp,pp_text)
    
    tight_layout()
    save_file_name = 'plot_inflow_Norway.pdf'
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
    print 'Load incresed by: ' + str(100*(balancing[0] + 0.022*mean(L))/mean(L))
    mismatch = mismatch - balancing - 0.022*mean(L) 
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
    
    year = array([h.year for h in num2date(t+datetime_offset-1/24.)])
    
    
    close(1);figure(1);clf()
    
    fill_between(t,-get_positive(-(mismatch - mismatch_r)),edgecolor=None,lw=0,facecolor='r',alpha=0.5)
    pp_discharge = Rectangle((0, 0), 1, 1, edgecolor=None,lw=0,facecolor='r',alpha=0.5)
    
    fill_between(t,get_positive(mismatch - mismatch_r),edgecolor=None,lw=0,facecolor='g',alpha=0.7)
    pp_charge = Rectangle((0, 0), 1, 1, edgecolor=None,lw=0,facecolor='g',alpha=0.7)
    
    fill_between(t,-get_positive(-mismatch),-get_positive(-(mismatch - mismatch_r)),edgecolor=None,lw=0,facecolor='r')
    pp_balancing = Rectangle((0, 0), 1, 1, edgecolor=None,lw=0,facecolor='r')
    
    #plot(t,mismatch,'k-',lw=.1,label='Mismatch')
    
    plot(t,storage_lake.virtual_two_way_storage.default_charge,color='k',lw=1,label='Default charge')
    pp_default = Line2D([0,1], [0,0], color='k',ls='-')
    
    plot(t,storage_lake.virtual_two_way_storage.default_charge-storage_lake.virtual_two_way_storage.ekstra_discharge,'k-',alpha=0.5,label='Ekstra discharge')
    pp_allowed = Line2D([0,1], [0,0], color='k',ls='-',alpha=0.5)
    
    #plot(t,mismatch_r,label='Mismatch_r')
    #fill_between(t,balancing,edgecolor=None,lw=0,facecolor='r')
    
    #plot(t,storage_lake.virtual_two_way_storage.power_cap_in,'k-',lw=.5,alpha=.5)
    #plot(t,-storage_lake.virtual_two_way_storage.power_cap_out,'k-',lw=.5,alpha=.5)
    
    #plot(t,storage_lake.virtual_two_way_storage.default_charge,'k-',lw=.5,label='Default charge')
    #plot(t,storage_lake.virtual_two_way_storage.default_charge-volume*0.5/365./24,'k--',lw=.5,label='Default charge (offset)')
    #plot(t,storage_lake.virtual_two_way_storage.power_balance/volume*100,'r-',lw=.5,label='Power balance')
    
    
    ylabel('Hourly mismatch [GW]')
    
    pp = [pp_charge,pp_discharge,pp_balancing,pp_default,pp_allowed]
    pp_text = ['Charge','Discharge','Non-hydro balancing','Default charge/discharge','Allowed discharge']
    legend(pp,pp_text,ncol=3,title='Norway (NO)')
    
    axis(xmin=0,xmax=N_days,ymin=-33,ymax=39)
    
    xticks(t[1+find(diff(year))],year[1+find(diff(year))],rotation=-45,ha='left')
    
    tight_layout()
    save_file_name = 'plot_isolated_Norway_mismatch.png'
    save_figure(save_file_name)
    
    close(1);figure(1);clf()
    
    fill_between(t,storage_lake.get_level()/volume,edgecolor=None,lw=0,facecolor='r',alpha=0.5,label='Prefilling')
    pp_level = Rectangle((0, 0), 1, 1, edgecolor=None,lw=0,facecolor='r',alpha=0.5)
    
    ## Plot real storage levels
    #storage_level_real_=kron(np.genfromtxt('./raw_data/meannor'),ones(7*24))
    #storage_level_real = zeros_like(t)
    #for i in arange(len(storage_level_real)):
    #    storage_level_real[i] = storage_level_real_[mod(i,len(storage_level_real_))]
    plot(t,storage_level/volume,'k-',label='Historical level')
    pp_historical = Line2D([0,1], [0,0], color='k',ls='-')
    
    plot(t,storage_lake.virtual_two_way_storage.median_level/volume,'k--',label='Median level')
    pp_median = Line2D([0,1], [0,0], color='k',ls='--')
    
    plot(t,storage_lake.virtual_two_way_storage.min_level/volume,'r-',label='Minimum level')
    pp_min = Line2D([0,1], [0,0], color='r',ls='-')
    
    axhline(1,color=(.5,.5,.5))
                    
    ylabel('Storage filling fraction')
    
    pp = [pp_historical,pp_level,pp_median,pp_min]
    pp_text = ['Historical level','Model level','Median level','Minimum level']
    legend(pp,pp_text,ncol=2,title='Norway (NO)')
    
    xticks(t[1+find(diff(year))],year[1+find(diff(year))],rotation=-45,ha='left')
    
    axis(xmin=0,xmax=N_days,ymin=0,ymax=1.199)
    
    tight_layout()
    save_file_name = 'plot_isolated_Norway_level.png'
    save_figure(save_file_name)    
    
    
    
    
    close(1);figure(1);clf()
    
    subplot(211)
    plot(t,mismatch-storage_lake.virtual_two_way_storage.default_charge)
    
    subplot(212)
    hist(mismatch-storage_lake.virtual_two_way_storage.default_charge,25)
    
    tight_layout()
    save_file_name = 'plot_isolated_Norway_2.png'
    save_figure(save_file_name)
    