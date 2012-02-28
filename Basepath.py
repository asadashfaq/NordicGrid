#Standard modules
from pylab import *
from scipy import *
import os
import sys

#Specific functions
import scipy.optimize as optimize

#Custom functions
from shortcuts import *
from SingleCountry import *

def main():

    year = array([2000,2010,2020,2030,2050])
    gamma = array([0.1,.25,0.5,.75,1.0])
    alpha_w = array([1.,1.,1.5/2.,2.2/3.,2.9/4.])
    #alpha_w = array([1.,1.,1.8/2.,2.5/3.,3.2/4.])
    gamma_wind = get_logistic_fit(year,alpha_w*gamma,year0=2000,plot_on=True,txtlabel='Wind',txttitle='Wind')
    gamma_solar = get_logistic_fit(year,(1-alpha_w)*gamma,year0=2000,plot_on=True,txtlabel='Solar',txttitle='Solar')

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1);figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.5])
    
    #plot(arange(amin(year),amax(year),1),gamma_wind+gamma_solar)
    plot(3.9*gamma_wind,3.9*gamma_solar)
    axis('scaled')
    axis(xmin=0,xmax=3.9*1.1,ymin=0,ymax=3.9*1.1)
    xlabel(r'Av. wind power [GWh/h]')
    ylabel(r'Av. solar PV power [GWh/h]')
    
    tight_layout()
    save_figure('TestFigure.png')


#
# get_optimal_path_and_logistic_fit([1980,2000,2010,2020,2050],[0.01,0.1,.25,0.5,1.0])
#
def get_optimal_path_and_logistic_fit(p_year, p_gamma, year=None, ISO='DK', CS=None, dgamma=.05):

    p_year, p_gamma = array(p_year,ndmin=1), array(p_gamma,ndmin=1)

    if year==None:
        year = arange(amin(p_year),amax(p_year),1)

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('DK')

    ## Find gamma vs. year
    gamma = get_logistic_fit(p_year,p_gamma,year=year,plot_on=True,txtlabel=ISO+'_test',txttitle=ISO+'_test')

    ## Find optimal mix vs gamma
    #gamma = arange(amin(p_gamma),amax(p_gamma),dgamma)
    alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L,Gw,Gs,gamma,returnall=True)
    p_alpha_w_opt = interp(p_gamma,gamma,alpha_w_opt)
    
    
    
    
    close(1);figure(1)
    plot(year,alpha_w_opt*gamma,'b')
    plot(year,(1-alpha_w_opt)*gamma,'y')
    
    plot(p_year,p_alpha_w_opt*p_gamma,'o')
    
    axis(ymin=0,ymax=1)
    
    tight_layout()
    save_figure('TestFigure.png')
    
    return p_alpha_w_opt

#
# plot_wind_solar_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0], p_alpha_w=[1.,1,.9,.75,.75])
#
def plot_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0=1980., year=None):

    year, gamma, alpha_w = get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0, year)

    close(1); figure(1)
    
    plot(year,gamma*alpha_w)
    plot(year,gamma*(1-alpha_w))
    
    plot(year,gamma,lw=2)
    
    plot(p_year,p_gamma,'o')
    
    axis(ymin=0,ymax=1)
    
    tight_layout()
    save_figure('TestFigure.png')

#
# get_wind_solar_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0], p_alpha_w=[1.,1,.9,.75,.75])
#
def get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0=1980., year=None):

    p_year, p_gamma, p_alpha_w = array(p_year,ndmin=1), array(p_gamma,ndmin=1), array(p_alpha_w,ndmin=1)

    if year==None:
        year = arange(amin(p_year),amax(p_year),1)

    f_logistic = lambda p, x, x0: p[0]*abs(p[2])*exp(p[1]*(x-x0))/(abs(p[2])+p[0]*(exp(p[1]*(x-x0))-1.)) # p is a lenth 3 array

    f_gamma_w = lambda pp, year: f_logistic(pp[0:3],year,year0)
    f_gamma_s = lambda pp, year: f_logistic(pp[3:6],year,year0)
    f_gamma = lambda pp, year: f_gamma_w(pp,year) + f_gamma_s(pp,year)
    f_alpha_w = lambda pp, year: f_gamma_w(pp,year)/f_gamma(pp,year)
    
    errfunc = lambda pp, year, gamma, alpha_w: concatenate([f_gamma(pp,year)-gamma,f_alpha_w(pp,year)-alpha_w]) # pp is a length 6 array.

    pp_0 = [amin(p_gamma), .01, amax(p_gamma), amin(p_gamma), .01, amax(p_gamma)] # Initial guess for the parameters
    pp_fit, success = optimize.leastsq(errfunc, pp_0[:], args=(p_year, p_gamma, p_alpha_w))

    return year, f_gamma(pp_fit,year), f_alpha_w(pp_fit,year)


def get_logistic_fit(p_year,p_gamma,year0=1980,year=None,plot_on=False,p_historical=None,txtlabel=None,txttitle=None):
    
    if year==None:
        year = arange(amin(p_year),amax(p_year),1)
    
    fitfunc = lambda p, x: p[0]*abs(p[3])*exp(p[1]*(x-year0))/(abs(p[3])+p[0]*(exp(p[1]*(x-year0))-1.))
    errfunc = lambda p, x, y, weight: (fitfunc(p, x) - y)/weight # Distance to the target function

    p_0 = [amin(p_gamma), .01, year0, amax(p_gamma)] # Initial guess for the parameters
    p_weight = ones(p_year.shape)
    p_fit, success = optimize.leastsq(errfunc, p_0[:], args=(p_year, p_gamma, p_weight))

    if plot_on==True:
        plot_logistic_fit(year,fitfunc(p_fit,year),p_year,p_gamma,p_historical,txtlabel,txttitle)

    return fitfunc(p_fit,year)
    
def plot_logistic_fit(year,gamma_fit,p_year,p_gamma,p_historical=None,txtlabel=None,txttitle=None):

    if p_historical==None:
        p_historical = array(zeros(p_year.shape),dtype=bool)

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    pp_hist = plot(array(p_year[find(p_historical)]),array(p_gamma[find(p_historical)]),'go')
    pp_target = plot(p_year[find(~p_historical)],p_gamma[find(~p_historical)],'ro')
    pp_fit = plot(year,gamma_fit,'k--',lw=1.5)
    
    axis(xmin=amin(year),xmax=2053,ymin=0,ymax=1.3)
    
    xlabel('Reference year')
    ylabel(r'Gross share of electricity demand ($\gamma_{'+txtlabel+'}$)')

    pp = concatenate([pp_hist,pp_target,pp_fit])
    pp_text = ['Historical values','Target values','Logistic fit']
    leg = legend(pp,pp_text,loc='upper left',title=txttitle)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    tight_layout(pad=.2)
    save_figure('plot_logistic_fit_' + txtlabel + '.png')