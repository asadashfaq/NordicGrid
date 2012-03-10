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

##Colors, should be moved to a color module
color_wind = (0.5,0.7,1.)
color_solar = (1.,.8,0.)

#
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2050], p_gamma=[0.023,0.121,.219,0.5,1.0])
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2050], p_gamma=[0.023,0.121,.219,0.5,1.0],CS=7.6,label='storage')
#
def plot_optimal_path_logistic_fit(p_year, p_gamma, ISO='DK', year0=1980., year=None, CS=None, rel_tol=1e-2,label=''):
    
    p_gamma = array(p_gamma,ndmin=1)
    
    year, gamma, alpha_w = get_optimal_path_logistic_fit(p_year, p_gamma, ISO, year0, year, CS, rel_tol)

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    wind_sum, solar_sum, surplus_sum, storage_sum = zeros_like(year), zeros_like(year), zeros_like(year), zeros_like(year)
    for i in arange(len(year)):
    
        wind = gamma[i]*alpha_w[i]*Gw*mean(L)
        solar = gamma[i]*(1-alpha_w[i])*Gs*mean(L)
        mismatch = (wind+solar) - L
    
        if CS!=None or CS==0:
            mismatch_r = get_policy_2_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = CS)[0]
        else:
            mismatch_r = mismatch
    
        curtailment = get_positive(mismatch_r)
        filling = get_positive(mismatch - mismatch_r)
        extraction = get_positive(-(mismatch - mismatch_r))
    
        wind_local = wind - (curtailment+filling)*wind/(wind+solar+1e-10)
        solar_local = solar - (curtailment+filling)*solar/(wind+solar+1e-10)
    
        wind_sum[i] = mean(wind_local)*365*24
        solar_sum[i] = mean(solar_local)*365*24
        surplus_sum[i] = mean(curtailment)*365*24
        storage_sum[i] = mean(extraction)*365*24
    

    #Set plot options	
    matplotlib.rcdefaults()
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1);
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.5])

    bar(year,wind_sum/1e3,color=color_wind,align='center',label='Wind')
    bar(year,solar_sum/1e3,bottom=wind_sum/1e3,color=color_solar,align='center',label='Solar')    
    if CS!=None:
        bar(year,storage_sum/1e3,bottom=(wind_sum+solar_sum)/1e3,color='green',align='center',label='Storage output')    
        bar(year,surplus_sum/1e3,bottom=(wind_sum+solar_sum+storage_sum)/1e3,color='red',align='center',label='Surplus - Storage output')    
    else:
        bar(year,surplus_sum/1e3,bottom=(wind_sum+solar_sum)/1e3,color='red',align='center',label='Surplus')

    plot(p_year,p_gamma*mean(L)*365*24/1e3,'ko',label='Target')

    for y, g in zip(p_year,p_gamma):
        text(y-1.,g*mean(L)*365*24/1e3,'{0:0.0f}%'.format(g*100),ha='right',va='bottom',weight='bold',fontsize=12)

    axis(xmin=amin(year)-1,xmax=amax(year)+1,ymin=0,ymax=1.05*amax(p_gamma*mean(L)*365*24/1e3))

    xlabel('Reference year')
    ylabel('Energy [TWh]')

    leg = legend(loc='upper left',title='Optimal build-up')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    tight_layout()
    save_file_name = 'plot_optimal_path_logistic_fit_'+ISO+'_'+label+'.pdf'
    save_figure(save_file_name)

#
# get_optimal_path_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0]);
#
def get_optimal_path_logistic_fit(p_year, p_gamma, ISO='DK', year0=1980., year=None, CS=None, rel_tol=1e-2):
    """Combines targets for gamma with the optimal path to get an optimal build-up of wind and solar. Both follow a logistic growth."""

    p_year, p_gamma = array(p_year,ndmin=1), array(p_gamma,ndmin=1)
    
    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)
        
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('DK')
    
    ## Initial estimate of gamma vs year
    p_alpha_w_opt = get_optimal_path_balancing(L,Gw,Gs,p_gamma,CS=CS)
    year, gamma, alpha_w = get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w_opt, year0, year)
    
    ## Iterate to get better match with the optimal mix. The gamma vs year estimate is used to provide a better and continous match to the optimal path. Otherwise few targets can result in large deviations from the optimal path.
    rel_err = 1000; rel_goodness_old = 1000; i=0
    while rel_err>rel_tol or i>10:
        i+=1
        gamma_w_old = gamma*alpha_w
        
        alpha_w_opt = get_optimal_path_balancing(L,Gw,Gs,gamma,CS=CS)
        year, gamma, alpha_w = get_wind_solar_logistic_fit(year, gamma, alpha_w_opt, year0, year)

        rel_goodness = amax(abs(gamma_w_old - gamma*alpha_w))
        rel_err = abs(rel_goodness - rel_goodness_old)
        rel_goodness_old = rel_goodness    
    
    return year, gamma, alpha_w

#
# get_wind_solar_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0], p_alpha_w=[1.,1,.9,.75,.75])
#
def get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0=1980., year=None):
    """Combine targets for both gamma and alpha_w in one fit. Both wind and solar are required to follow a logistic growth."""

    p_year, p_gamma, p_alpha_w = array(p_year,ndmin=1), array(p_gamma,ndmin=1), array(p_alpha_w,ndmin=1)

    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)

    f_logistic = lambda p, x, x0: abs(p[0])*abs(p[2])*exp(abs(p[1])*(x-x0))/(abs(p[2])+abs(p[0])*(exp(abs(p[1])*(x-x0))-1.)) # p is a lenth 3 array

    f_gamma_w = lambda pp, year: f_logistic(pp[0:3],year,year0)
    f_gamma_s = lambda pp, year: f_logistic(pp[3:6],year,year0)
    f_gamma = lambda pp, year: f_gamma_w(pp,year) + f_gamma_s(pp,year)
    f_alpha_w = lambda pp, year: f_gamma_w(pp,year)/f_gamma(pp,year)
    
    #errfunc = lambda pp, year, gamma, alpha_w: concatenate([(f_gamma(pp,year)-gamma),(f_alpha_w(pp,year)-alpha_w)]) # pp is a length 6 array.
    
    #This error function treats wind and solar in the same way. The in the one above emphasis can be placed on either mix or on gamma.
    errfunc = lambda pp, year, gamma, alpha_w: concatenate([(f_alpha_w(pp,year)*f_gamma(pp,year)-alpha_w*gamma),((1-f_alpha_w(pp,year))*f_gamma(pp,year)-(1-alpha_w)*gamma)]) # pp is a length 6 array.

    pp_0 = [amin(p_alpha_w*p_gamma), .01, amax(p_alpha_w*p_gamma), amin((1-p_alpha_w)*p_gamma), .01, amax((1-p_alpha_w)*p_gamma)] # Initial guess for the parameters
    pp_fit, success = optimize.leastsq(errfunc, pp_0[:], args=(p_year, p_gamma, p_alpha_w))

    return year, f_gamma(pp_fit,year), f_alpha_w(pp_fit,year)






# May be obsolete. I have to make sure though.
def get_logistic_fit(p_year,p_gamma,year0=1980,year=None,plot_on=False,p_historical=None,txtlabel=None,txttitle=None):
    
    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)
    
    fitfunc = lambda p, x: p[0]*abs(p[3])*exp(p[1]*(x-year0))/(abs(p[3])+p[0]*(exp(p[1]*(x-year0))-1.))
    errfunc = lambda p, x, y, weight: (fitfunc(p, x) - y)/weight # Distance to the target function

    p_0 = [amin(p_gamma), .01, year0, amax(p_gamma)] # Initial guess for the parameters
    p_weight = ones(p_year.shape)
    p_fit, success = optimize.leastsq(errfunc, p_0[:], args=(p_year, p_gamma, p_weight))

    if plot_on==True:
        plot_logistic_fit(year,fitfunc(p_fit,year),p_year,p_gamma,p_historical,txtlabel,txttitle)

    return fitfunc(p_fit,year)
    
# May be obsolete. I have to make sure though.    
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