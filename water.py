#! /usr/bin/env python
from pylab import *
from scipy import *
import numpy as np

## Custom modules and functions
from shortcuts import get_high_low_fft, get_positive 
from SingleCountry import get_ISET_country_data
from scipy.stats.mstats import mquantiles

def get_inflow_one_year(path_data='./raw_data/',returnall=False):
    #ion()
    prod=np.genfromtxt(path_data+"prodnor",dtype=int)
    water=np.genfromtxt(path_data+"meannor")
    soco=np.zeros(8761)
    soc=np.zeros(8761)
    inflow=np.ones(8760)
    soco[0]=51951000 #2M more than next week, extrapolated
    for i in range(8760):
        soco[i+1]=soco[i]-prod[i]
    # plot(soco)

    correction=0
    for i in range(51):
        obj=1000*water[i]
        mes=np.mean(soco[168*i:168*(i+1)])
        x=(obj-mes-correction)/168
        if x<0:x=0
        inflow[i*168:168*(i+1)]*=x
        correction=obj-mes

    soc[0]=51951000
    for i in range(8760):
        soc[i+1]=soc[i]-prod[i]+inflow[i]

    #Fill something in first week
    inflow[0:7*24] = inflow[8*24]
    storage_level = soc

    if returnall:
        return inflow/1e3, storage_level/1e3
    else:
        return inflow/1e3
    
from datetime import datetime, timedelta, date
def tofirstdayinisoweek(year, week):
    ret = datetime.strptime('%04d-%02d-1' % (year, week), '%Y-%W-%w')
    if date(int(year), 1, 4).isoweekday() > 4:
        ret -= timedelta(days=7)
    return ret    
    
def get_storage_level_Norway(t=None,datetime_offset=None,path_data='./raw_data/',returnall=False):
    """ Returns the fractional storage filling level."""
    
    storage_level_data = np.genfromtxt(path_data+'ReservoirLevelsWeekly_1990_to_2012_NO.csv',skip_header=1,delimiter=',').transpose()
    
    year = storage_level_data[0]
    week = storage_level_data[1]
    level = storage_level_data[2]
    
    # The date should be changed to the Wednesday 1 pm CET in the week. See http://www.ssb.no/vannmag_en/
    date = date2num(array([tofirstdayinisoweek(year, week) for year, week in zip(year,week) ]))
    
    if datetime_offset==None:
        datetime_offset = date[0]
    if t==None:
        t = arange(24*(date[-1] - date[0]))/24.
    
    level = interp(t+datetime_offset,date,level,left=nan,right=nan)
    
    if returnall:
        return level, t, datetime_offset
    else:
        return level    
        
def get_median_storage_level_Norway(quantile=0.5, returnall=False):

    level, t, datetime_offset = get_storage_level_Norway(returnall=True)
    years = np.unique(array([date.year for date in num2date(t + datetime_offset)])) 
    
    t_year = arange(0,24*366,1)/24.
    level_year = zeros((len(years),len(t_year)))
    for i in arange(len(years)):
        start = date2num(datetime(year=years[i],month=1,day=1,hour=0))
        level_year[i] = get_storage_level_Norway(t=t_year,datetime_offset=start)
    
    level_year_ = copy(level_year)
    level_year_[isnan(level_year)] = -.1
    level_median = mquantiles(level_year_,quantile,limit=(0,1),axis=0)
    
    N_hour = 24*7
    for i in arange(len(level_median)):
        level_median[i] = get_high_low_fft(level_median[i],N_hour)[0]
        
    if returnall:
        return level_median, level_year
    else:
        return level_median

def get_median_storage_power_Norway(quantile=0.5, smooth=True, cyclic_level=True, returnall=False):

    level, t, datetime_offset = get_storage_level_Norway(returnall=True)
    years = np.unique(array([date.year for date in num2date(t + datetime_offset)])) 
    
    t_year = arange(0,24*366,1)/24.
    level_year, power_year = zeros((len(years),len(t_year))), zeros((len(years),len(t_year)))
    for i in arange(len(years)):
        start = date2num(datetime(year=years[i],month=1,day=1,hour=0))
        level_year[i] = get_storage_level_Norway(t=t_year,datetime_offset=start)
        power_year[i] = concatenate([diff(level_year[i]),[level_year[i][-1]-level_year[i][-2]]])
    
    power_year_ = copy(power_year)
    power_year_[isnan(power_year)] = 2
    
    power_median = mquantiles(power_year_,quantile,limit=(-1,1),axis=0)
    
    
    N_hour = 24*7
    level_median = get_median_storage_level_Norway(quantile)
    power_median_level = zeros_like(power_median)
    for i in arange(len(power_median)):
        if smooth:
            power_median[i] = get_high_low_fft(power_median[i],N_hour)[0]
        
        if cyclic_level:
            power_median[i] = power_median[i] - sum(power_median[i])/float(len(power_median[i]))
        
        power_median_level[i] = cumsum(power_median[i])+level_median[i][0]
                
    if returnall:
        return power_median, power_year, power_median_level
    else:
        return power_median
    
def plot_stuff(quantile=0.5):

    level_median, level_year = get_median_storage_level_Norway(returnall=True,quantile=quantile)
    power_median, power_year, power_median_level = get_median_storage_power_Norway(returnall=True,smooth=True,quantile=quantile)
    
    close(1); figure(1)
    
    subplot(311)
    for level_year_ in level_year:
        plot(level_year_)
        
    plot(level_median[0],'k-',lw=3)
    #plot(cumsum(power_median[0])+level_median[0][0],'r--',lw=2)
    plot(power_median_level[0],'r--',lw=2)

    
    subplot(312)
    
    for i in arange(len(power_year)):
        #plot(cumsum(power_year[i])+level_year[i][0])
        plot(power_year[i])
        
    #plot(cumsum(power_median[0])+level_median[0],'k-',lw=3)
    plot(power_median[0],'k-',lw=3)
    
    plot(diff(level_median[0]),'r--',lw=2)
    
    subplot(313)
    
    plot(power_median[0][:-1]-diff(level_median[0]),'k-',lw=2)
    #plot(cumsum(power_median[0][:-1])-cumsum(diff(level_median[0])),'r-',lw=2)
    #plot(cumsum(power_median[0][:-1])+level_median[0][0]-level_median[0][:-1],'g-',lw=.5)
    
    savefig('TestFigure.png')

def get_inflow(t=None,datetime_offset=None, storage_capacity=80e3, returnall=False, correction=1.04):
    """Estimated hourly inflow of water to the Norwigian storage lakes."""
    
    #Load data
    t_, L, Gw, Gs, datetime_offset_, datalabel = get_ISET_country_data('NO')
    
    if t == None:
        t = t_
    if datetime_offset == None:
        datetime_offset = datetime_offset_

    storage_level = storage_capacity*get_storage_level_Norway(t,datetime_offset)
    
    L_week = get_high_low_fft(L, tau=7*24)[0]
    L_week = interp(t,t_,L_week)
    
    inflow = get_positive(diff(storage_level) + L_week[:-1])
    inflow = concatenate([inflow,[inflow[-1]]])
    
    ## Correction for Norway's export of hydropower
    inflow = correction*inflow
    
    if returnall:
        return inflow, storage_level, t, datetime_offset
    else:
        return inflow
    
def get_inflow_old(t,datetime_offset=None, returnall=False):
    """Estimated hourly inflow of water to the Norwigian storage lakes."""


    inflow_one_year, storage_level_one_year = get_inflow_one_year(returnall=True)
    inflow = zeros_like(t)
    inflow_basic = zeros_like(t)
    storage_level_basic = zeros_like(t)

    # datetime_offset is actually not used at the moment. To use it would mean replacing mult_year by a year specific multiplier. (http://www.ssb.no/elektrisitetaar_en/)
    if datetime_offset == None:
        datetime_offset = datestr2num('1-1-2000')
    
    for i in arange(len(t)):
        if mod(i,7*24)==0:
            i_offset = randn(1)*7.*24.  # Mix up the time series somewhat.
            mult_week = 1 + randn(1)*.1  # Adds a random weekly variation.
        if mod(i,365*24)==0:   
            mult_year = 1 + randn(1)*.1  # Adds a random yearly variation.
            
        inflow[i] = mult_year*mult_week*inflow_one_year[int(mod(i+i_offset,len(inflow_one_year)))]
        inflow_basic[i] = inflow_one_year[mod(i,len(inflow_one_year))]
        storage_level_basic[i] = storage_level_one_year[mod(i,len(storage_level_one_year))]
        
    if returnall:
        return inflow, inflow_basic, storage_level_basic
    else:
        return inflow
        
    