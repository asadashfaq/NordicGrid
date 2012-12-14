#Standard modules
from pylab import *
from scipy import *
import os
import sys

#Specific functions
import scipy.optimize as optimize
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.special import erfc

#Custom functions
from shortcuts import *
from SingleCountry import *

##Colors, should be moved to a color module
color_wind = (0.5,0.7,1.)
color_solar = (1.,.8,0.)

#
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2050], p_gamma=[0.023,0.121,.219,0.5,1.0])
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2050], p_gamma=[0.023,0.121,.219,0.5,1.0],CS=7.6,label='storage',use_CS_opt=True)
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2050], p_gamma=[0.023,0.121,.219,0.5,1.0],CS=12,label='storage_12h_BalOptMix')
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2050], p_gamma=[0.023,0.121,.219,0.5,1.0],CS=12,label='storage_12h',use_CS_opt=True)
#
# plot_optimal_path_logistic_fit(p_year=[1990,2000,2010,2020,2030,2050], p_gamma=[0.023,0.121,.219,0.5,0.75,1.0],p_alpha_w=ones(6),label='WindOnly')
def plot_optimal_path_logistic_fit(p_year, p_gamma, ISO='DK', year0=1980., year=None, CS=None, rel_tol=1e-2, p_alpha_w=None,label='', use_CS_opt=False):
    
    p_gamma = array(p_gamma,ndmin=1)
    
    if p_alpha_w==None:
        if use_CS_opt:
            year, gamma, alpha_w = get_optimal_path_logistic_fit(p_year, p_gamma, ISO, year0, year, CS, rel_tol)
        else:
            year, gamma, alpha_w = get_optimal_path_logistic_fit(p_year, p_gamma, ISO, year0, year, CS=None, rel_tol=rel_tol)
    else:
        year, gamma, alpha_w = get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0, year)
    
    

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
    

    ## Energy produced
    close(1); figure(1);
    # gcf().set_size_inches([6.5,4.5])

    bar(year,wind_sum/1e3,color=color_wind,align='center',label='Wind')
    bar(year,solar_sum/1e3,bottom=wind_sum/1e3,color=color_solar,align='center',label='Solar PV')    
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


def get_wind_solar_use(year,gamma,alpha_w,ISO='DK',CS=None):

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
        
    return wind_sum, solar_sum, surplus_sum, storage_sum 

def plot_capacity_factor_wind(p_year=[1900,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1990,1991,1992,1993,1997,2002,2009,2015,2020,2030,2050,2100],p_CF_wind=[0.05,0.05,0.13,0.10,0.19,0.22,0.15,0.08,0.27,0.15,0.16,0.27,0.19,0.40,0.32,0.30,0.32,0.35,0.46,0.46,0.48,0.50,0.51,0.51],textlabel='Denmark'):

    
    

    CF_wind = get_capacity_factor_wind(p_year,p_CF_wind)
    
    year = linspace(amin(p_year),amax(p_year),100)
    
    close(1); figure(1)
    
    plot(year,CF_wind(year),'k--',label='Interpolation')
    plot(p_year,p_CF_wind,'kx',label='Data/predictions')
    
    xlabel('Reference year')
    ylabel('Wind capacity factor')
    
    axis(ymin=0,ymax=.52,xmin=1980,xmax=2051)
    
    legend(title=textlabel,loc='upper left')

    tight_layout()
    save_file_name = 'plot_capacity_factor_wind.pdf'
    save_figure(save_file_name)

def plot_capacity_factor_solar(p_year=[1960,1990,2010,2015,2030,2050,2100],p_CF_solar=[0.08,0.097,0.097,0.103,0.148,0.17,0.17],textlabel='Denmark'):

    CF_solar = get_capacity_factor_solar(p_year,p_CF_solar)
    
    year = linspace(amin(p_year),amax(p_year),100)
    
    close(1); figure(1)
    
    plot(year,CF_solar(year),'k--',label='Interpolation')
    plot(p_year,p_CF_solar,'kx',label='Data/predictions')
    
    xlabel('Reference year')
    ylabel('Solar capacity factor')
    
    axis(ymin=0,ymax=.27,xmin=1980,xmax=2051)
    
    legend(title=textlabel)

    tight_layout()
    save_file_name = 'plot_capacity_factor_solar.pdf'
    save_figure(save_file_name)

def get_capacity_factor_wind(p_year=[1900,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1990,1991,1992,1993,1997,2002,2009,2015,2020,2030,2050,2100],p_CF_wind=[0.05,0.05,0.13,0.10,0.19,0.22,0.15,0.08,0.27,0.15,0.16,0.27,0.19,0.40,0.32,0.30,0.32,0.35,0.46,0.46,0.48,0.50,0.51,0.51]):

    #CF_wind = UnivariateSpline(p_year,p_CF_wind,s=0.1)
    year = linspace(amin(p_year),amax(p_year),100)

    CF_wind = lambda x: interp(x,year,get_high_low_fft(interp(year,p_year,p_CF_wind),5)[0])

    return CF_wind

def get_capacity_factor_solar(p_year=[1960,1990,2010,2015,2030,2050,2100],p_CF_solar=[0.08,0.097,0.097,0.103,0.148,0.17,0.17]):

    #CF_wind = UnivariateSpline(p_year,p_CF_wind,s=0.1)
    year = linspace(amin(p_year),amax(p_year),100)

    CF_solar = lambda x: interp(x,year,get_high_low_fft(interp(year,p_year,p_CF_solar),10)[0])

    return CF_solar

# plot_installed_capacity_optimal_path(p_year=[1990,2000,2010,2020,2030,2050], p_gamma=[0.023,0.121,.219,0.5,0.75,1.0],savelabel='OptMix')
#
# plot_installed_capacity_optimal_path(p_year=[1990,2000,2010,2020,2030,2050], p_gamma=[0.023,0.121,.219,0.5,0.75,1.0],p_alpha_w=ones(6),savelabel='WindOnly')

def plot_installed_capacity_optimal_path(p_year, p_gamma, ISO='DK', year0=1980., year=None, CS=None, rel_tol=1e-2, p_alpha_w=None, savelabel=''):

    p_gamma = array(p_gamma,ndmin=1)
    
    # Installed capacity in DK
    p_installed = array([326,2390,3752,nan,nan,nan])
    
    # Move this around
    CF_wind = get_capacity_factor_wind() #UnivariateSpline([1980,1985,1988,1990,2000,2010,2011,2020,2050,2100],[0.0845,0.107,0.159,0.22,0.25,0.30,0.40,.40,.40,.45])
    CF_solar = get_capacity_factor_solar() #UnivariateSpline([1990,2010,2015,2030,2050,2100],[0.097,0.097,0.103,0.148,0.17,0.17])
    
    lifetime_wind = 20. #year.
    lifetime_solar = 30. #year.
    
    ## Load/get data 
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    if p_alpha_w==None:
        year, gamma, alpha_w = get_optimal_path_logistic_fit(p_year, p_gamma, ISO, year0, year, CS, rel_tol)
    else:
        year, gamma, alpha_w = get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0, year)
    
    annual_TWh = mean(L)*365*24/1e3
    
    energy_wind = gamma*alpha_w
    energy_solar = gamma*(1-alpha_w)
    
    installed_solar, installed_wind = zeros_like(year), zeros_like(year)
    energy_installed_solar, energy_installed_wind = zeros_like(year), zeros_like(year) #The energy produced by already installed capacity.
    d_capacity_solar, d_capacity_wind = zeros_like(year), zeros_like(year) # New installed capacity per year. Capacity that is dismanteled is not included.
    d_capacity_new_solar, d_capacity_new_wind = zeros_like(year), zeros_like(year)
    
    ## Initial conditions:
    f = lambda year, year_0, lifetime: (year>=year_0)*erfc((year-(year_0+lifetime))/sqrt(lifetime))/2.
    
    # energy_installed_solar, energy_installed_wind = energy_solar[0]*f(year,year[0],lifetime_solar), energy_wind[0]*f(year,year[0],lifetime_wind)
    energy_installed_solar, energy_installed_wind = energy_solar[0]*exp(-(year-year[0])**2/(sqrt(2)*lifetime_solar)**2), energy_wind[0]*exp(-(year-year[0])**2/(sqrt(2)*lifetime_wind)**2)
    installed_solar, installed_wind = energy_installed_solar*mean(L)*1e3/CF_solar(year[0]), energy_installed_wind*mean(L)*1e3/.2 #CF_wind(year[0])
    
    for i in arange(1,len(year)):
        
        ## Solar power
        d_energy_solar = energy_solar[i] - energy_installed_solar[i] # Additional energy to be produced in year[i].
        d_energy_new_solar = energy_solar[i] - energy_solar[i-1]
        d_capacity_solar[i] = d_energy_solar*mean(L)*1e3/CF_solar(year[i]) # Additional capacity requred to produce d_energy in year[i].
        d_capacity_new_solar[i] = d_energy_new_solar*mean(L)*1e3/CF_solar(year[i]) 
        energy_installed_solar = energy_installed_solar + d_energy_solar*f(year,year[i],lifetime_solar)
        installed_solar = installed_solar + d_capacity_solar[i]*f(year,year[i],lifetime_solar)  
        
        ## Wind power
        d_energy_wind = energy_wind[i] - energy_installed_wind[i] # Additional energy to be produced in year[i].
        d_energy_new_wind = energy_wind[i] - energy_wind[i-1]
        d_capacity_wind[i] = d_energy_wind*mean(L)*1e3/CF_wind(year[i]) # Additional capacity requred to produce d_energy in year[i].
        d_capacity_new_wind[i] = d_energy_new_wind*mean(L)*1e3/CF_wind(year[i]) 
        energy_installed_wind = energy_installed_wind + d_energy_wind*f(year,year[i],lifetime_wind)
        installed_wind = installed_wind + d_capacity_wind[i]*f(year,year[i],lifetime_wind)

    d_capacity_wind[0], d_capacity_new_wind[0] = 2*d_capacity_wind[1]-d_capacity_wind[2], 2*d_capacity_new_wind[1]-d_capacity_new_wind[2]
    d_capacity_solar[0], d_capacity_new_solar[0] = 2*d_capacity_solar[1]-d_capacity_solar[2], 2*d_capacity_new_solar[1]-d_capacity_new_solar[2]
    
    #installed_wind = gamma*alpha_w*mean(L)*1e3/CF_wind(year)
    #installed_solar = gamma*(1-alpha_w)*mean(L)*1e3/CF_solar(year)
    
    wind_sum, solar_sum, surplus_sum, storage_sum = get_wind_solar_use(year,gamma,alpha_w,ISO,CS)
    wind_sum, solar_sum, surplus_sum, storage_sum = wind_sum/1e3, solar_sum/1e3, surplus_sum/1e3, storage_sum/1e3
    
    close(1); figure(1);
    
    gcf().set_size_inches([6.7,7.5]) # Make room for three subplots.

    ## Annual energy production.
    subplot(311)
    
    surplus_color = (1,.2,.2)
    
    bar(year,wind_sum,color=color_wind,align='center',label='Wind (used)',lw=0.9)
    bar(year,solar_sum,bottom=wind_sum,align='center',color=color_solar,label='Solar PV (used)',lw=0.9)
    
    
    bar(year,energy_wind*annual_TWh-wind_sum,bottom=solar_sum+wind_sum,hatch='\\',color=surplus_color,align='center',label='Wind (surplus)',edgecolor='k',lw=0.9)
    
    
    bar(year,energy_solar*annual_TWh-solar_sum,bottom=solar_sum+wind_sum+(energy_wind*annual_TWh-wind_sum),hatch='//',color=surplus_color,align='center',label='Solar PV (surplus)',lw=0.9)
    
    plot(p_year,p_gamma*mean(L)*365*24/1e3,'ko',label='Target')
    for y, g in zip(p_year,p_gamma):
        text(y-1.,g*mean(L)*365*24/1e3,'{0:0.0f}%'.format(g*100),ha='right',va='bottom',weight='bold',fontsize=8)

    
    xlabel('Reference year')
    ylabel('Energy per year [TWh]')
    
    axis(xmin=1990,xmax=2051,ymin=0,ymax=annual_TWh*1.1)
    
    legend(loc='upper left')
    
    ## Accumulated installed capacity.
    subplot(312)
    
    bar(year,installed_wind/1e3,align='edge',color=color_wind,label='Wind',width=0.4,lw=0.9)
    bar(year+0.4,installed_solar/1e3,align='edge',color=color_solar,label='Solar PV',width=0.4,lw=0.9)
    
    # plot(p_year,p_installed/1e3,'ko',label='Target')

    
    xlabel('Reference year')
    ylabel('Installed capacity [GW]')
    
    axis(xmin=1990,xmax=2051)
    
    legend(loc='upper left')
    
    ## New capacity per year
    subplot(313)
    
    ## Wind
    bar(year,d_capacity_new_wind,color=color_wind,align='edge',label='Wind ({0:.0f} yr lifetime)'.format(lifetime_wind),width=0.4,lw=1.)
    bar(year,d_capacity_wind - d_capacity_new_wind,bottom=d_capacity_new_wind,color=color_wind,alpha=0.75,align='edge',label='Wind (replacement)',width=0.4,lw=0.5)
    
    ## Solar
    bar(year+0.4,d_capacity_new_solar,color=color_solar,align='edge',label='Solar PV ({0:.0f} yr lifetime)'.format(lifetime_solar),width=0.4,lw=1.)
    bar(year+0.4,d_capacity_solar-d_capacity_new_solar,bottom=d_capacity_new_solar,color=color_solar,alpha=0.75,align='edge',label='Solar PV (replacement)',width=0.4,lw=0.5)

    
    xlabel('Reference year')
    ylabel('New capacity [MW/yr]')
    
    axis(xmin=1990,xmax=2051,ymin=0,ymax=450)
    
    legend(loc='upper left')
    
    tight_layout()
    save_file_name = 'plot_installed_capacity_optimal_path_'+savelabel+'.pdf'
    save_figure(save_file_name)

#
# get_optimal_path_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0]);
#
def get_optimal_path_logistic_fit(p_year, p_gamma, ISO='DK', year0=1980., year=None, CS=None, rel_tol=1e-2):
    """Combines targets for gamma with the optimal path to get an optimal build-up of wind and solar. Both follow a logistic growth."""

    p_year, p_gamma = array(p_year,ndmin=1), array(p_gamma,ndmin=1)
    
    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)
        
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
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
# get_wind_solar_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0], p_alpha_w=[1,1,1,1,.99])
#
def get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0=1980., year=None):
    """Combine targets for both gamma and alpha_w in one fit. Both wind and solar are required to follow a logistic growth."""

    p_year, p_gamma, p_alpha_w = array(p_year,ndmin=1), array(p_gamma,ndmin=1), array(p_alpha_w,ndmin=1)

    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)

    f_logistic = lambda p, x, x0: abs(p[0])*abs(p[2])*exp(abs(p[1])*(x-x0))/(abs(p[2])+abs(p[0])*(exp(abs(p[1])*(x-x0))-1.)) # p is a lenth 3 array

    if all(p_alpha_w==0): # Pure solar scenario.
        f_gamma_w = lambda pp, year: zeros_like(year)
    else:
        f_gamma_w = lambda pp, year: f_logistic(pp[0:3],year,year0)
    
    if all(p_alpha_w==1): # Pure wind scenario.
        f_gamma_s = lambda pp, year: zeros_like(year)
    else:
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