#
#  SingleCountry.py
#  Maximum local integration wind and solar power generation for a region/country.
#
#  Created by Gorm Bruun Andresen on 24/11/2011.
#  Copyright (c) 2011 Department of Engineering, University of Aarhus. All rights reserved.
#

#Standard modules
from pylab import *
from scipy import *
import os
import sys

#Custom functions
from shortcuts import *
from Database_v1 import get_data_countries, get_data_regions
from MortenStorage import get_policy_2_storage, get_storage

#Specific functions, not sure if all are actually used
from scipy.optimize import brentq
from scipy.optimize import fmin
from scipy.optimize import leastsq
from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles
from mpl_toolkits.axes_grid1 import make_axes_locatable
import calendar #Used to place day labels on plots

## Colors, should be moved to a color module
color_wind = (0.5,0.7,1.)
color_solar = (1.,.8,0.)
bg_color = (.75,.0,.0)
color_edge = (.4,.4,.4)
color_opt = (0.53,0.73,0.37)
golden_ratio = 1.61803399


## Standard figure size:
figure_size = [6.5,4.3]



################################################################################################
################################################################################################



def plot_duration_curves(ISO='DK', gamma=[0,.25,.5,.75,1.], alpha_w=.8, CS=None, N_bins=211):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    duration_residual_load, power_residual_load = zeros((len(gamma),N_bins)), zeros((len(gamma),N_bins))
    duration_upregulation, power_upregulation = zeros((len(gamma),N_bins)), zeros((len(gamma),N_bins))
    duration_downregulation, power_downregulation = zeros((len(gamma),N_bins)), zeros((len(gamma),N_bins))
    for i in arange(len(gamma)):
        duration_residual_load[i], power_residual_load[i] = get_duration_curve_residual_load(ISO, gamma[i], alpha_w, CS, N_bins)
        duration_upregulation[i], power_upregulation[i] = get_duration_curve_upregulation(ISO, gamma[i], alpha_w, CS, N_bins)
        duration_downregulation[i], power_downregulation[i] = get_duration_curve_downregulation(ISO, gamma[i], alpha_w, CS, N_bins)
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    ### Residual load
    close(1);figure(1);clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    for i in arange(len(gamma)):
        plot(duration_residual_load[i],power_residual_load[i]*mean(L),label=str(gamma[i]))
    
    axis(xmin=1,xmax=1.05e4,ymin=0,ymax=1.05*max(L))
    
    xscale('log')
    
    xlabel('Hours')
    ylabel('Power')
    
    legend(title='Duration curve',loc='lower left')
    
    tight_layout()
    save_figure('plot_duration_curves_residual_load_' + ISO + '.pdf')
    
    ### Upregulation
    close(1);figure(1);clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    for i in arange(len(gamma)):
        plot(duration_upregulation[i],power_upregulation[i]*mean(L),label=str(gamma[i]))
    
    axis(xmin=1,xmax=1.05e4,ymin=0,ymax=1.05*max(L)/2.)
    
    xscale('log')
    
    xlabel('Hours')
    ylabel('Power')
    
    legend(title='upregulation',loc='upper right')
    
    tight_layout()
    save_figure('plot_duration_curves_upregulation_' + ISO + '.pdf')
    
    ### Downregulation
    close(1);figure(1);clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    for i in arange(len(gamma)):
        plot(duration_downregulation[i],power_downregulation[i]*mean(L),label=str(gamma[i]))
    
    axis(xmin=1,xmax=1.05e4,ymin=0,ymax=1.05*max(L)/2.)
    
    xscale('log')
    
    xlabel('Hours')
    ylabel('Power')
    
    legend(title='downregulation',loc='upper right')
    
    tight_layout()
    save_figure('plot_duration_curves_downregulation_' + ISO + '.pdf')



################################################################################################
################################################################################################


    
def get_duration_curve_residual_load(ISO='DK', gamma=0, alpha_w=.8, CS=None, N_bins=211):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    bin_edges = N_bins #linspace(0,2,N_bins+1)

    residual_load = get_positive(-get_mismatch(L, Gw, Gs, gamma, alpha_w,CS))

    bin_counts, bin_edges = np.histogram(residual_load,bin_edges)
    bin_center = bin_edges[:-1] + diff(bin_edges)/2.
    
    duration = 365.*24. - cumsum(bin_counts)*365.*24./len(L)
    power = bin_center

    return duration, power


################################################################################################
################################################################################################



def get_duration_curve_upregulation(ISO='DK', gamma=0, alpha_w=.8, CS=None, N_bins=211):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    bin_edges = N_bins #linspace(0,2,N_bins+1)

    residual_load = get_positive(-get_mismatch(L, Gw, Gs, gamma, alpha_w,CS))

    upregulation = get_positive(diff(residual_load))

    bin_counts, bin_edges = np.histogram(upregulation,bin_edges)
    bin_center = bin_edges[:-1] + diff(bin_edges)/2.
    
    duration = 365.*24. - cumsum(bin_counts)*365.*24./len(L)
    power = bin_center

    return duration, power



################################################################################################
################################################################################################



def get_duration_curve_downregulation(ISO='DK', gamma=0, alpha_w=.8, CS=None, N_bins=211):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    bin_edges = N_bins #linspace(0,2,N_bins+1)

    residual_load = get_positive(-get_mismatch(L, Gw, Gs, gamma, alpha_w,CS))

    downregulation = get_positive(-diff(residual_load))

    bin_counts, bin_edges = np.histogram(downregulation,bin_edges)
    bin_center = bin_edges[:-1] + diff(bin_edges)/2.
    
    duration = 365.*24. - cumsum(bin_counts)*365.*24./len(L)
    power = bin_center

    return duration, power


################################################################################################
################################################################################################



def plot_phase_space_residual_load(ISO='DK', gamma=0, alpha_w=.8, CS=None):
    
    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    residual_load = get_positive(-get_mismatch(L, Gw, Gs, gamma, alpha_w,CS))
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    ### Residual load
    close(1);figure(1);clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    plot(residual_load[:-1],diff(residual_load),'.')
    
    axis(xmin=0,xmax=2,ymin=-0.5,ymax=.5)
    
    xlabel('Residual load')
    ylabel('Ramp 1h')
    
    tight_layout()
    save_figure('plot_phase_space_residual_load_' + 'g{0:0.2f}_'.format(gamma) + ISO + '.png')



################################################################################################
################################################################################################



def plot_storage_balancing_synergy_old(ISO='DK', gamma=[.75,1.,1.25], CS=linspace(0,24,5)):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    balancing_mean = zeros((len(gamma),len(CS)))
    print 'Starting {0:.0f} calculations:'.format(len(gamma)*len(CS)),
    for j in arange(len(gamma)):
        for i in arange(len(CS)):
            print i,
            sys.stdout.flush()

            #Optimal mix
            alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = get_optimal_mix_balancing(L, Gw, Gs, gamma[j], CS=CS[i], returnall=True)

            balancing_mean[j][i] = get_positive(-mismatch_opt).mean()

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    figure(1);clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    for j in arange(len(gamma)):
        plot(CS,balancing_mean[j] - gamma[j]*(gamma[j]<=1.))
    
    axis(xmin=0,xmax=amax(CS),ymin=0)
    
    xlabel('Storage size [av.l.h.]')
    ylabel('Balancing power [av.l.h.]')
    
    tight_layout()
    save_figure('plot_storage_balancing_synergy_' + ISO + '.png')



################################################################################################
################################################################################################



def plot_country_balancing_at_optimal_mix_vs_gamma(ISO='DK', gamma=linspace(0,2.05,51), quantiles=[.50,.90,.99,1.00], linestyle=[':','-.','--','-'], CS=None, ROI=[.2,.5]):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    balancing_power_mean, excess_power_mean = zeros(len(gamma)), zeros(len(gamma))
    balancing_power_mean_wind, excess_power_mean_wind = zeros(len(gamma)), zeros(len(gamma))
    balancing_power_quantiles, excess_power_quantiles = zeros((len(gamma),len(quantiles))), zeros((len(gamma),len(quantiles)))
    balancing_power_quantiles_wind, excess_power_quantiles_wind = zeros((len(gamma),len(quantiles))), zeros((len(gamma),len(quantiles)))
    for i in arange(len(gamma)):
        #Optimal mix
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = get_optimal_mix_balancing(L, Gw, Gs, gamma[i], CS=CS, returnall=True)
        
        balancing_power_mean[i], excess_power_mean[i] = mean(get_positive(-mismatch_opt)), mean(get_positive(mismatch_opt))
        balancing_power_quantiles[i], excess_power_quantiles[i] = mquantiles(get_positive(-mismatch_opt),quantiles), mquantiles(get_positive(mismatch_opt),quantiles)
    
        #Wind only
        mismatch_wind = get_mismatch(L, Gw, Gs, gamma[i], alpha_w=1.,CS=CS)
        balancing_power_mean_wind[i], excess_power_mean_wind[i] = mean(get_positive(-mismatch_wind)), mean(get_positive(mismatch_wind))
        balancing_power_quantiles_wind[i], excess_power_quantiles_wind[i] = mquantiles(get_positive(-mismatch_wind),quantiles), mquantiles(get_positive(mismatch_wind),quantiles)
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10    
                
    figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,6])
    
    subplot(211)
    #Highlight ROI
    if not ROI==None:
        fill_betweenx([0,10],ROI[1]*ones(2),ROI[0],lw=0,color='k',alpha=.1)

    plot(gamma,balancing_power_mean_wind,'-',color=color_wind,lw=3, label='Wind')
    for i in arange(len(quantiles)):
        plot(gamma,balancing_power_quantiles_wind.transpose()[i],ls=linestyle[i],lw=2,color=color_wind,label=r'{0:.0f}% quantile'.format(quantiles[i]*100))
    
    
    plot(gamma,balancing_power_mean,'-',color=(.5,.5,.5),lw=2.5,label='Mean')
    for i in arange(len(quantiles)):
        plot(gamma,balancing_power_quantiles.transpose()[i],ls=linestyle[i],lw=1.5,color='k',label=r'{0:.0f}% quantile'.format(quantiles[i]*100))



    axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=1.65)
    
    ylabel('Balancing power [av.l.h.]')
    xlabel(r'Gross share of electricity demand $\gamma_{'+ISO+'}$')
    
    #leg = legend(loc='upper right',title='Balancing ('+ISO+'):');
    #ltext  = leg.get_texts();
    #setp(ltext, fontsize='small')    # the legend text fontsize

    add_duplicate_yaxis(gcf(),unit_multiplier=mean(L),label='[GW]')
    
    subplot(212)
    
    #Highlight ROI
    if not ROI==None:
        fill_betweenx([0,10],ROI[1]*ones(2),ROI[0],lw=0,color='k',alpha=.1)
    
    plot(gamma,excess_power_mean_wind,'-',color=color_wind,lw=3, label='Wind')
    for i in arange(len(quantiles)):
        plot(gamma,excess_power_quantiles_wind.transpose()[i],ls=linestyle[i],lw=2,color=color_wind,label=r'{0:.0f}% quantile'.format(quantiles[i]*100))
    
    plot(gamma,excess_power_mean,'-',color=(.5,.5,.5),lw=2.5,label='Mean')
    for i in arange(len(quantiles)):
        plot(gamma,excess_power_quantiles.transpose()[i],ls=linestyle[i],lw=1.5,color='k',label=r'{0:.0f}% quantile'.format(quantiles[i]*100))
        
    axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=2.05)

    ylabel('Excess power [av.l.h.]')
    xlabel(r'Gross share of electricity demand $\gamma_{'+ISO+'}$')

    #leg = legend(loc='upper left',title='Excess ('+ISO+'):');
    #ltext  = leg.get_texts();
    #setp(ltext, fontsize='small')    # the legend text fontsize

    add_duplicate_yaxis(gcf(),unit_multiplier=mean(L),label='[GW]')

    tight_layout(pad=0.5,h_pad=2)
    save_file_name = 'plot_country_balancing_at_optimal_mix_vs_gamma_'+'CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)


################################################################################################
################################################################################################



def add_duplicate_yaxis(figure_handle,unit_multiplier=10.,label='New label',tickFormatStr='%.1f',invert_axis=False):

    majorFormatter = FormatStrFormatter(tickFormatStr)

    figure(figure_handle.number)
    
    ax1 = gca()
    yticks_ = unit_multiplier*ax1.get_yticks()[find((ax1.get_yticks()<=ax1.get_ylim()[1]) * (ax1.get_yticks()>=ax1.get_ylim()[0]))]
    
    divider = make_axes_locatable(ax1)
    ax2 = divider.append_axes("right", "0%", pad="0%")

    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    
    axis(ymin=unit_multiplier*ax1.get_ylim()[0],ymax=unit_multiplier*ax1.get_ylim()[1])
    ylabel(label)
    xticks([0],[''])
    
    yticks(yticks_,yticks_)
    if invert_axis:
        ax2.invert_yaxis()
    
    ax2.yaxis.set_major_formatter(majorFormatter)
 

################################################################################################
################################################################################################


       
             
###
# plot_country_optimal_mix_vs_gamma('DK', gamma=linspace(0,2.05,31))
# 6h storage: plot_country_optimal_mix_vs_gamma('DK', gamma=linspace(0,2.05,31),CS=6)
#
# plot_country_optimal_mix_vs_gamma('DK', gamma=linspace(0,1.05,31),p_interval=[0.01,0.05],linespec=['-','-'])
#
#
#(1.,.53,.20)
def plot_country_optimal_mix_vs_gamma(ISO='DK', gamma=linspace(0,2.05,11), p_interval=[0.01,0.05,0.25], CS=None, ROI=[.2,.5], linespec=['--','-.',':'], color=[(0.53,0.73,0.37),(1.,.82,.20),(.90,.27,.20)]):
    """Legacy function from solarDK.py. Should be updated at some point."""

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, res_load_sum_1p = zeros(len(gamma)), zeros((len(gamma),2)), zeros(len(gamma)), zeros(len(gamma))
    lower_bound, upper_bound, res_load_sum_p = zeros((len(p_interval),len(gamma))), zeros((len(p_interval),len(gamma))), zeros((len(p_interval),len(res_load_sum_1p)))
    for j in arange(len(p_interval)):
        for i in arange(len(gamma)):
            alpha_w_opt[i], alpha_w_opt_1p_interval[i], res_load_sum_opt[i], mismatch_opt, res_load_sum_1p[i] = get_optimal_mix_balancing(L, Gw, Gs, gamma[i], p_interval=p_interval[j], CS=CS, returnall=True, normalized=True)

        lower_bound[j] = alpha_w_opt_1p_interval.transpose()[0]
        upper_bound[j] = alpha_w_opt_1p_interval.transpose()[1]
        res_load_sum_p[j] = res_load_sum_1p/len(L)

    #Filter optimal path
    mask = arange(len(gamma) - amax([argmin(lower_bound[0][::-1]),argmax(upper_bound[0][::-1])]),len(gamma))
    mask = mask[find(abs(diff(alpha_w_opt[mask],2))/diff(gamma[mask][:-1])<1)] #Filter noise in optimal path.

    res_load_sum_opt = res_load_sum_opt/len(L)

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    #gcf().set_size_inches([5.25,6])
    #gcf().set_size_inches([5.25,3.5])
    gcf().set_size_inches([6.5,4.3])
    
    #Upper panel
    #ax1 = axes([.11,.565,.885,.42])
    #subplot(211)
    plot(gamma[mask],alpha_w_opt[mask],'w-',lw=2)

    fill_between(gamma,ones(len(gamma)),color=bg_color,edgecolor=(0,0,0,0))

    pp = list(zeros(1+len(p_interval)))
    for j in arange(len(p_interval))[::-1]:
        fill_between(gamma,lower_bound[j],upper_bound[j],color=color[j],lw=1,edgecolor=(0,0,0,0))
        pp[j] = Rectangle((0, 0), 1, 1, color=color[j],lw=1)
        
        plot(gamma,lower_bound[j],linespec[j],color='w',lw=1)
        plot(gamma,upper_bound[j],linespec[j],color='w',lw=1)

    pp[-1] = Rectangle((0, 0), 1, 1, color=bg_color,lw=1)
    pp = list(pp)

    ## Present state
    gamma_present,alpha_w_present, label_present = .22, 1.0, '2010'
    plot(gamma_present,alpha_w_present,'o',ms=10,mfc='k',mec=(.8,.8,.8),mew=1)
    text(gamma_present+.025,alpha_w_present-.025,label_present,va='top',ha='left',weight='semibold',fontsize=10)

    ##Scenario guide lines
    ##axvline(.2,ls='--',color='k')
    #build_path = lambda gamma, alpha_w_inf, gamma_0=0.2: (gamma_0 + (gamma-gamma_0)*alpha_w_inf)/(gamma+1e-10)
    #plot(gamma,build_path(gamma,0.9),'k-')
    #text(1.5,build_path(1.5,0.9)+.01,'10% solar',va='baseline',weight='demibold')
    #plot(gamma,build_path(gamma,0.8),'k-')
    #text(1.5,build_path(1.5,0.8)+.01,'20% solar',va='baseline',weight='demibold')
    #plot(gamma,build_path(gamma,0.7),'k-')
    #text(1.5,build_path(1.5,0.7)+.01,'30% solar',va='baseline',weight='demibold')



    axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=1)
    xlabel(r'Gross share of electricity demand $\gamma_{'+ISO+'}$')
    ylabel(r'Wind fraction $\alpha^{'+ISO+'}_w$')

    #txtlabels = ['0-1 pp','1-5 pp','5-25 pp','>25 pp']
    txtlabels = ['0-1 pp','1-5 pp',' >5 pp']

    leg = legend(pp,txtlabels,loc='lower right',ncol=2,title='Deviation from optimal mix');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    add_duplicate_yaxis(gcf(),unit_multiplier=1.,label=r'Solar fraction $\alpha^{'+ISO+'}_s$',invert_axis=True)

    tight_layout()
    save_file_name = 'plot_country_optimal_mix_vs_gamma_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)

    ########
    # Triangle plot of optimal mix map
    ####

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.5])
    
    conv_x = lambda gamma, alpha_w: mean(L)*365*24/1e3*gamma*alpha_w
    conv_y = lambda gamma, alpha_w: mean(L)*365*24/1e3*gamma*(1-alpha_w)

    plot(conv_x(gamma[mask],alpha_w_opt[mask]),conv_y(gamma[mask],alpha_w_opt[mask]),'w-',lw=2)

    #fill_between(conv_x(gamma,amax(gamma)*ones(len(gamma))),amax(gamma)*ones(len(gamma)),color=bg_color,edgecolor=(0,0,0,0))
    fill(mean(L)*365*24/1e3*concatenate([gamma,[0]]),mean(L)*365*24/1e3*concatenate([zeros(gamma.shape),[amax(gamma)]]),color=bg_color,edgecolor=(0,0,0,0))

    pp = list(zeros(1+len(p_interval)))
    for j in arange(len(p_interval))[::-1]:
        
        
        x = concatenate([conv_x(gamma,lower_bound[j]),conv_x(gamma,upper_bound[j])[::-1]])
        y = concatenate([conv_y(gamma,lower_bound[j]),conv_y(gamma,upper_bound[j])[::-1]])
        fill(x,y,color=color[j],lw=1,edgecolor=(0,0,0,0))
        
        
        #ub = interp(conv_x(gamma,lower_bound[j]),conv_x(gamma,upper_bound[j]),conv_y(gamma,upper_bound[j]))
        #fill_between(conv_x(gamma,lower_bound[j]),conv_y(gamma,lower_bound[j]),ub,color=color[j],lw=1,edgecolor=(0,0,0,0))
        pp[j] = Rectangle((0, 0), 1, 1, color=color[j],lw=1)
        
        plot(conv_x(gamma,lower_bound[j]),conv_y(gamma,lower_bound[j]),linespec[j],color='w',lw=1)
        plot(conv_x(gamma,upper_bound[j]),conv_y(gamma,upper_bound[j]),linespec[j],color='w',lw=1)

    pp[-1] = Rectangle((0, 0), 1, 1, color=bg_color,lw=1)
    pp = list(pp)

    ## Guide lines
    plot(.25*array([0.,mean(L)*365*24/1e3]),.25*array([mean(L)*365*24/1e3,0]),'k--',lw=1)
    text(.05*mean(L)*365*24/1e3,.25*mean(L)*365*24/1e3-.025*mean(L)*365*24/1e3,r'25%',weight='semibold',fontsize=10)
    
    plot(.5*array([0.,mean(L)*365*24/1e3]),.5*array([mean(L)*365*24/1e3,0]),'k--',lw=1)
    text(.05*mean(L)*365*24/1e3,.5*mean(L)*365*24/1e3-.025*mean(L)*365*24/1e3,r'50%',weight='semibold',fontsize=10)
    
    plot(.75*array([0.,mean(L)*365*24/1e3]),.75*array([mean(L)*365*24/1e3,0]),'k--',lw=1)
    text(.05*mean(L)*365*24/1e3,.75*mean(L)*365*24/1e3-.025*mean(L)*365*24/1e3,r'75%',weight='semibold',fontsize=10)

    plot(1.*array([0.,mean(L)*365*24/1e3]),1.*array([mean(L)*365*24/1e3,0]),'k--',lw=1)
    text(.05*mean(L)*365*24/1e3,1.*mean(L)*365*24/1e3-.025*mean(L)*365*24/1e3,r'100%',weight='semibold',fontsize=10)
    

    ## Present state
    gamma_present,alpha_w_present, label_present = .22*mean(L)*365*24/1e3, 0.0, '2010'
    plot(gamma_present,alpha_w_present,'o',ms=10,mfc='k',mec='w',mew=1)
    text(gamma_present+.025*mean(L)*365*24/1e3,alpha_w_present+.025*mean(L)*365*24/1e3,label_present,va='bottom',ha='left',weight='semibold',fontsize=10)

    axis(xmin=0,xmax=mean(L)*365*24/1e3*amax(gamma),ymin=0,ymax=mean(L)*365*24/1e3*amax(gamma))
    axis('scaled')
    xlabel(r'Wind energy [TWh/yr]')
    ylabel(r'Solar energy [TWh/yr]')
    
    ax=gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    txtlabels = ['0 - 0.35 TWh','0.35 - 1.7 TWh',' >1.7 TWh']

    leg = legend(pp,txtlabels,loc='upper right',ncol=2,title='Additional surplus');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

   # add_duplicate_yaxis(gcf(),unit_multiplier=1.,label=r'Solar fraction $\alpha^{'+ISO+'}_s$',invert_axis=True)

    tight_layout()
    save_file_name = 'plot_country_optimal_mix_vs_gamma_triangle_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)

    #####
    # Constant gamma crosssections
    #####

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])
    
    gamma_cross_section = array([0.25,0.5,1.0])
    alpha_w_cross_section = linspace(0,1,21)
    
    for i in arange(len(gamma_cross_section)):
        total_surplus = mean(L)*365*24/1e3*(get_balancing(L, Gw, Gs, gamma_cross_section[i], CS=CS, alpha=alpha_w_cross_section).transpose()[0]/len(L) - (1. - amin([gamma_cross_section[i],1])))
        
        plot(alpha_w_cross_section, total_surplus,'k-',lw=1.5)
        
        text(.15,interp(.15,alpha_w_cross_section,total_surplus)+.5,r'{0:0.0f}%'.format(gamma_cross_section[i]*100),weight='semibold',fontsize=10)
    
    xlabel('Wind share of total VRE energy')
    ylabel('Total wind and solar surplus [TWh/yr]')
    axis(xmin=0,xmax=1,ymin=0)
    
    tight_layout()
    save_file_name = 'plot_country_optimal_mix_vs_gamma_cross_sections_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)

    #####
    # Decrease in average balancing vs gamma
    #####

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

    #Highlight ROI
    if not ROI==None:
        fill_betweenx([0,10],ROI[1]*ones(2),ROI[0],lw=0,color='k',alpha=.1)

    pp_full = plot([0,1],[1,0],'--',color=(.5,.5,.5))

    pp_solar = plot(gamma,get_balancing(L, Gw, Gs, gamma, CS=CS, alpha=0.)[0]/len(L),'-',color=color_solar,lw=2)
    pp_wind = plot(gamma,get_balancing(L, Gw, Gs, gamma,  CS=CS, alpha=1.)[0]/len(L),'-',color=color_wind,lw=2)

    #pp = list(zeros(len(p_interval)))
    #for j in arange(len(p_interval))[::-1]:
    #    pp_ = plot(gamma,res_load_sum_p[j],linespec[j],color='k')
    #    pp[j] = pp_[0]
        
    pp_opt = plot(gamma,res_load_sum_opt,'k-')

    #axvline(.2,ls='--',color='k')
                
    axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=1)
    ylabel('Av. residual load [av.h.l.]')
    xlabel(r'Gross share of electricity demand $\gamma_{'+ISO+'}$')

    pp = [pp_full[0],pp_opt[0],pp_wind[0],pp_solar[0]]
    #pp2.extend(pp)
    #txtlabels = ['Optimal mix',r'100% wind','100% solar','1 pp contour','5 pp contour','25 pp contour']
    txtlabels = [r'Full integration',r'Optimal mix',r'100% wind',r'100% solar']

    leg = legend(pp,txtlabels,loc='upper right',ncol=2);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    add_duplicate_yaxis(gcf(),unit_multiplier=mean(L),label='[GWh]')

    tight_layout()
    save_file_name = 'plot_country_transition_vs_gamma_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)
    
    
    #####
    # Av. local excess vs gamma
    #####

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

    #Highlight ROI
    if not ROI==None:
        fill_betweenx([0,10],ROI[1]*ones(2),ROI[0],lw=0,color='k',alpha=.1)

    balancing_min = (1.-gamma)*(gamma<=1)

    pp_solar = plot(gamma,get_balancing(L, Gw, Gs, gamma, CS=CS, alpha=0.)[0]/len(L) - balancing_min,'-',color=color_solar,lw=2)
    pp_wind = plot(gamma,get_balancing(L, Gw, Gs, gamma,  CS=CS, alpha=1.)[0]/len(L) - balancing_min,'-',color=color_wind,lw=2)

    #pp = list(zeros(len(p_interval)))
    #for j in arange(len(p_interval))[::-1]:
    #    pp_ = plot(gamma,res_load_sum_p[j],linespec[j],color='k')
    #    pp[j] = pp_[0]
        
    pp_opt = plot(gamma,res_load_sum_opt - balancing_min,'k-')

    #axvline(.2,ls='--',color='k')
                
    axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=.33)
    ylabel('Av. excess balancing [av.h.l.]')
    xlabel(r'Gross share of electricity demand $\gamma_{'+ISO+'}$')

    pp = [pp_opt[0],pp_wind[0],pp_solar[0]]
    #pp2.extend(pp)
    #txtlabels = ['Optimal mix',r'100% wind','100% solar','1 pp contour','5 pp contour','25 pp contour']
    txtlabels = [r'Optimal mix',r'100% wind',r'100% solar']

    leg = legend(pp,txtlabels,loc='upper left',ncol=1);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    add_duplicate_yaxis(gcf(),unit_multiplier=mean(L),label='[GWh]')

    tight_layout()
    save_file_name = 'plot_country_excess_vs_gamma_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)



################################################################################################
################################################################################################



#
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,11), CS=linspace(0,14,11))
#
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,11), CS=concatenate([[0],linspace(1e-2,14,10)]))
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,51), CS=concatenate([[0],linspace(1e-2,14,50)]))
#
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,11), CS=concatenate([[0],linspace(1e-2,14,10)]),alpha_w=1,txtlabel='wind-only',savelabel='wind_only')
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,51), CS=concatenate([[0],linspace(1e-2,14,50)]),alpha_w=1,txtlabel='wind-only',savelabel='wind_only')
#
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,11), CS=concatenate([[0],linspace(1e-2,14,10)]),alpha_w=0.6,txtlabel='Seasonal mix',savelabel='season_mix')
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,51), CS=concatenate([[0],linspace(1e-2,14,50)]),alpha_w=0.6,txtlabel='Seasonal mix',savelabel='season_mix')
#
# plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,51), CS=concatenate([[0],linspace(1e-2,14,50)]),alpha_w=None,txtlabel='opt. mix',savelabel='opt_mix')
#
#   Hydrogen storage efficiencies: eta_in=0.4, eta_out=0.4
#   =============
#
#   plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,51), CS=concatenate([[0],linspace(1e-2,14,50)]),alpha_w=None,txtlabel='opt. mix',savelabel='opt_mix_H2_storage', eta_in=0.4, eta_out=0.4)
#
#
def plot_value_of_storage(ISO='DK', gamma=linspace(0,1.05,11), CS=linspace(0,27,11), alpha_w=None, txtlabel='', savelabel='', eta_in=1, eta_out=1):

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    annual_TWh = mean(L)*365*24/1e3
    
    ## Baseline, no storage  
    ## XXX: I should check if this is true in all cases.
    if alpha_w=='season':
        alpha_w_path, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=NaN, returnall=True, eta_in=eta_in, eta_out=eta_out)
    elif alpha_w==None:
        alpha_w_path, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=None, returnall=True, eta_in=eta_in, eta_out=eta_out)
    else:
        ## Here, a constant alpha_w is assumed.
        res_load_sum_0 = get_balancing(L, Gw, Gs, gamma, alpha_w, CS=None)[0]/len(L)
        alpha_w_path = alpha_w*ones_like(gamma)    

    Gamma_, CS_ = meshgrid(gamma,CS)
    alpha_w_ = kron(array(ones_like(gamma),ndmin=2).transpose(),alpha_w_path)
    res_load_sum_0_ = kron(array(ones_like(gamma),ndmin=2).transpose(),res_load_sum_0)
    
    Surplus_ = zeros_like(Gamma_)
    for i in arange(len(Gamma_.flat)):
        ## XXX: In fact, res_load_sum is returned.
        Surplus_.flat[i] = get_balancing(L, Gw, Gs, Gamma_.flat[i], alpha_w_.flat[i], CS_.flat[i], eta_in=eta_in, eta_out=eta_out)[0]/len(L) + 0.001

    ## XXX: Gamma - load covered by VRES = surplus VRES (not corrected for eta_in, eta_out). (hourly average values, normalized)
    Deviation_from_target = Gamma_ - (1. - Surplus_)
    ## XXX: Deviation_from_target[0] is the values for no storage (assuming that CS includes 0).
    Storage_output = kron(array(ones_like(CS),ndmin=2).transpose(),Deviation_from_target[0]) - Deviation_from_target

    from matplotlib.colors import ListedColormap, BoundaryNorm
    color=[(0.53,0.73,0.37), (1.,.82,.20), bg_color]
    cmap = ListedColormap(color)
    norm = BoundaryNorm([0, 0.01, 0.05, 1], cmap.N)
    
    close(1); figure(1); clf()

    contourf(Gamma_*mean(L)*365*24/1e3, CS_*mean(L), Deviation_from_target,[0,.01,0.05,1.],cmap=cmap,norm=norm)

    dx = 0.02*annual_TWh
    y = amax(CS)*mean(L)/6.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)

    xlabel('Wind plus solar energy ('+txtlabel+') [TWh/yr]')
    ylabel('Storage volume [GWh]')
    
    pp = [Rectangle((0, 0), 1, 1, color=c,lw=1) for c in color]
    txtlabels = ['0 - 0.35 TWh','0.35 - 1.7 TWh',' >1.7 TWh']

    legend(pp,txtlabels,loc='upper left',ncol=2,title='Surplus after storage');

    tight_layout()
    save_file_name = 'plot_value_of_storage_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)
    
    ## Storage throughput per volume
    close(1); figure(1); clf()
    
    intervals = [0,10,50,100,150,200,300]
    pp_c = contourf(Gamma_*mean(L)*365*24/1e3, CS_*mean(L), Storage_output*365*24/CS_,intervals,cmap=cm.jet)
    
    dx = 0.02*annual_TWh
    y = amax(CS)*mean(L)/6.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)
    
    xlabel('Wind plus solar energy ('+txtlabel+') [TWh/yr]')
    ylabel('Storage volume [GWh]')    
    
    pp = [Rectangle((0, 0), 1, 1, color=pp_c_.get_facecolor()[0],lw=1) for pp_c_ in pp_c.collections]
    intervals = intervals[:len(pp)]
    txtlabels = concatenate([['{0:.0f} - {1:.0f}'.format(intervals[i],intervals[i+1]) for i in arange(len(intervals)-1)],['>{0:.0f}'.format(intervals[-1])]])

    legend(pp,txtlabels,loc='upper left',ncol=2,title=r'Cycle count [yr$^{-1}$]');
    
    tight_layout()
    save_file_name = 'plot_value_of_storage_2_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)
    



################################################################################################
################################################################################################



def get_storage_buildup_fixed_cycle_number(ISO='DK', gamma=linspace(0,1,5), N_cycles=100, gain=0.90, alpha_w=None, eta_charge=1., eta_discharge=1.):

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    ## Baseline, no storage
    if alpha_w=='season':
        alpha_w_path, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=NaN, returnall=True)
    elif alpha_w==None:
        alpha_w_path, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=None, returnall=True)
    else:
        res_load_sum_0 = get_balancing(L, Gw, Gs, gamma, alpha_w, CS=None)[0]/len(L)
        alpha_w_path = alpha_w*ones_like(gamma) 
    
    fitfunc_N_cycles = lambda gamma, alpha_w, CS, acc: get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w, CS, acc,eta_charge=eta_charge,eta_discharge=eta_discharge,returnall=True)[-1]*365*24
    
    CS_fit, P_in, P_out = zeros_like(gamma), zeros_like(gamma), zeros_like(gamma)
    for i in arange(len(gamma)):
        if gain*fitfunc_N_cycles(gamma[i],alpha_w_path[i],1e-6,acc=0)-N_cycles>0:
            CS_fit[i] = brentq(lambda CS: gain*fitfunc_N_cycles(gamma[i],alpha_w_path[i],CS,acc=0)-N_cycles, 1e-6, 365*24)
        else:
            CS_fit[i] = 1e-6 #Could be changed such that N_cycles was reduced to the maximum number of cycles at CS=1e-6.
    
    #print fitfunc_N_cycles(gamma[i],alpha_w,CS_fit,1-gain)
    
        Storage_benefit, P_in_, P_out_ = get_min_storage_cap_alt(L, Gw, Gs, gamma[i], alpha_w_path[i], CS_fit[i], 1-gain, eta_charge=eta_charge,eta_discharge=eta_discharge)
    
        P_in[i] = mean(P_in_)
        P_out[i] = mean(P_out_)
    
        print ' '
        print gamma[i], alpha_w_path[i], CS_fit[i], P_in[i], P_out[i]
    
    return CS_fit, P_in, P_out



################################################################################################
################################################################################################



def get_color_and_linestyle(alpha_w,N_cycles=None,gamma=None):

    ## Pick color
    if alpha_w==None:
        color = 'k'
        mixtxt = ' (bal. opt. mix)'
    elif alpha_w==1:
        color = color_wind
        mixtxt = ' (wind only)'
    else:
        color = 'g'
    
    ## Pick line style
    if not N_cycles==None:
        if N_cycles==50:
            ls = '-'
        elif N_cycles==100:
            ls = '--'
        else:
            ls = '-.'
    else:
        if gamma==.5:
            ls = '-'
        elif gamma==.75:
            ls = '--'
        else:
            ls = '-.'

    return color, ls, mixtxt


################################################################################################
################################################################################################



##
# plot_storage_buildup_fixed_cycle_number(gamma=linspace(0,1.05,51))
#    
def plot_storage_buildup_fixed_cycle_number(ISO='DK', gamma=linspace(0,1.05,5), N_cycles=[50,100,200,50,100], gain=0.90, alpha_w=[None,None,None,1,1], eta_charge=1., eta_discharge=1., txtlabel='', savelabel=''):

    N_cycles = array(N_cycles,ndmin=1)
    alpha_w = array(alpha_w,ndmin=1)
    
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    annual_TWh = mean(L)*365*24/1e3

    CS_fit_, P_in_, P_out_ = [], [], []
    for i in arange(len(N_cycles)):
        CS_fit, P_in, P_out = get_storage_buildup_fixed_cycle_number(ISO, gamma, N_cycles[i], gain, alpha_w[i], eta_charge, eta_discharge)
        
        CS_fit_.append(CS_fit)
        P_in_.append(P_in)
        P_out_.append(P_out)
    
    CS_fit_, P_in_, P_out_ = array(CS_fit_), array(P_in_), array(P_out_)
    print CS_fit_, P_in_, P_out_
    print CS_fit_[0], P_in_[0], P_out_[0]
    
    
    ## Plot storage energy capacity build-up
    close(1); figure(1); clf()

    for i in arange(len(N_cycles)):
        
        color, ls, mixtxt = get_color_and_linestyle(alpha_w[i],N_cycles[i])
                
        label_ = '{0:.0f}'.format(N_cycles[i]) + mixtxt
        plot(gamma*annual_TWh,mean(L)*CS_fit_[i],color=color,ls=ls,label=label_)

    
    dx = 0.02*annual_TWh
    y = 14*mean(L)/3.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)

    xlabel('Wind plus solar energy [TWh/yr]')
    ylabel('Storage volume [GWh]')

    axis(xmin=0,xmax=amax(gamma)*mean(L)*365*24/1e3,ymin=0,ymax=14*mean(L))

    legend(loc='upper left',title=r'Cycle count [yr$^{-1}$]',ncol=1)

    tight_layout()
    save_file_name = 'plot_storage_buildup_fixed_cycle_number_EnergyCap_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)

    
    ## Plot charging capacities
    close(1); figure(1); clf()

    for i in arange(len(N_cycles)):
        color, ls, mixtxt = get_color_and_linestyle(alpha_w[i],N_cycles[i])
        label_ = '{0:.0f}'.format(N_cycles[i]) + mixtxt
        
        plot(gamma*annual_TWh,1e3*mean(L)*P_in_[i],color=color,ls=ls,label=label_)

    xlabel('Wind plus solar energy [TWh/yr]')
    ylabel('Charging power [MW]')

    dx = 0.02*annual_TWh
    y = 2050/3.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)

    axis(xmin=0,xmax=amax(gamma)*annual_TWh,ymin=0,ymax=2050)

    legend(loc='upper left',title=r'Cycle count [yr$^{-1}$]',ncol=1)

    tight_layout()
    save_file_name = 'plot_storage_buildup_fixed_cycle_number_ChargingPower_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)

    ## Plot discharging capacities
    close(1); figure(1); clf()

    for i in arange(len(N_cycles)):
        color, ls, mixtxt = get_color_and_linestyle(alpha_w[i],N_cycles[i])
        label_ = '{0:.0f}'.format(N_cycles[i]) + mixtxt
        
        plot(gamma*annual_TWh,1e3*mean(L)*P_out_[i],color=color,ls=ls,label=label_)

    xlabel('Wind plus solar energy [TWh/yr]')
    ylabel('Discharging power [MW]')

    dx = 0.02*annual_TWh
    y = 2050/3.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)

    axis(xmin=0,xmax=amax(gamma)*annual_TWh,ymin=0,ymax=2050)

    legend(loc='upper left',title=r'Cycle count [yr$^{-1}$]',ncol=1)

    tight_layout()
    save_file_name = 'plot_storage_buildup_fixed_cycle_number_DischargingPower_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)
    
    ## Plot charging and discharging capacities Vs energy capacity ratio
    close(1); figure(1); clf()

    for i in arange(len(N_cycles)):
        color, ls, mixtxt = get_color_and_linestyle(alpha_w[i],N_cycles[i])
        label_ = '{0:.0f}'.format(N_cycles[i]) + mixtxt
        
        plot(gamma*annual_TWh,(1e3*mean(L)*P_in_[i])/(mean(L)*CS_fit_[i]),color=color,ls=ls,label=label_)

    dx = 0.02*annual_TWh
    y = 850/3.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)

    xlabel('Wind plus solar energy [TWh/yr]')
    ylabel('Charging ratio [MW/GWh]')

    axis(xmin=0,xmax=amax(gamma)*annual_TWh,ymin=0,ymax=850)

    legend(loc='upper right',title=r'Cycle count [yr$^{-1}$]',ncol=1)

    tight_layout()
    save_file_name = 'plot_storage_buildup_fixed_cycle_number_ChargingVsVolume_Ratio_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)

    ## Plot charging and discharging capacities Vs energy capacity ratio
    close(1); figure(1); clf()

    for i in arange(len(N_cycles)):
        color, ls, mixtxt = get_color_and_linestyle(alpha_w[i],N_cycles[i])
        label_ = '{0:.0f}'.format(N_cycles[i]) + mixtxt
        
        plot(gamma*annual_TWh,(1e3*mean(L)*P_out_[i])/(mean(L)*CS_fit_[i]),color=color,ls=ls,label=label_)

    dx = 0.02*annual_TWh
    y = 260/3.
    plot_vertical_line_and_label(0.25*annual_TWh,y,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,y,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,y,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,y,r'100%',dx)

    xlabel('Wind plus solar energy [TWh/yr]')
    ylabel('Discharging ratio [MW/GWh]')

    axis(xmin=0,xmax=amax(gamma)*annual_TWh,ymin=0,ymax=260)

    legend(loc='upper left',title=r'Cycle count [yr$^{-1}$]',ncol=1)

    tight_layout()
    save_file_name = 'plot_storage_buildup_fixed_cycle_number_DischargingVsVolume_Ratio_'+savelabel+'_'+ISO+'.pdf'
    save_figure(save_file_name)
    


################################################################################################
################################################################################################


    
def plot_value_of_storage_alt(ISO='DK', gamma=[.25,.50,.75,1.00], CS=[1,15,30,60], alpha_w=None):

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    ## Baseline, no storage
    if alpha_w==None:
        alpha_w_path, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=None, returnall=True)
    else:
        res_load_sum_0 = get_balancing(L, Gw, Gs, gamma, alpha_w, CS=None)[0]/len(L)
        alpha_w_path = alpha_w*ones_like(gamma)

    Gamma_, CS_ = meshgrid(gamma,CS)
    alpha_w_ = kron(array(ones_like(gamma),ndmin=2).transpose(),alpha_w_path)
    res_load_sum_0_ = kron(array(ones_like(gamma),ndmin=2).transpose(),res_load_sum_0)
    
    Surplus_ = zeros_like(Gamma_)
    for i in arange(len(Gamma_.flat)):
        Surplus_.flat[i] = (res_load_sum_0_.flat[i] - get_balancing(L, Gw, Gs, Gamma_.flat[i], alpha_w_.flat[i], CS_.flat[i])[0])/len(L)*mean(L*365*24)/1e3


    #Set plot options	
    matplotlib.rcParams['font.size'] = 10
    
    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

#    contourf(Gamma_, CS_, Surplus_,mean(L*365*24)/1e3*array([0,0.01,0.05,1]))
    contourf(Gamma_, CS_, Surplus_)

    xlabel('Gross share')
    ylabel('Storage volume')
    
    colorbar()

    tight_layout()
    save_file_name = 'plot_value_of_storage_'+ISO+'.pdf'
    save_figure(save_file_name)    
    


################################################################################################
################################################################################################



def plot_surplus_bar(ISO='DK',gamma_bar=array([.25,.50,.75,1.]),alpha_w=1,CS=None,N_gamma=111,returnall=True):

    surplus_bar = zeros(gamma_bar.shape)
    for i in arange(len(gamma_bar)):
        L, wind_local, solar_local, curtailment, filling, extraction = get_compare_VRES_load(ISO, gamma_bar[i], alpha_w, CS)

        surplus_bar[i] = mean(curtailment)*365*24/1e3

    gamma = linspace(amin(gamma_bar)-.1,amax(gamma_bar)+.1,N_gamma)
    surplus = zeros(gamma.shape)
    for i in arange(len(gamma)):
        L, wind_local, solar_local, curtailment, filling, extraction = get_compare_VRES_load(ISO, gamma[i], alpha_w, CS)

        surplus[i] = mean(curtailment)*365*24/1e3
    

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

    plot(gamma*100,surplus,'k-',lw=1.5,zorder=0)
    bars = bar(gamma_bar*100,surplus_bar,align='center',width=10,color='r')
    
    #Label bars:
    for i in arange(len(bars)):
        height = bars[i].get_height()
        if height<=0.4*365*24/1e3:
            height = height+0.5
        else:
            height = height/2
            
        #gca().text(bars[i].get_x()+bars[i].get_width()/2., height, '{0:0.1f}%\nof\n{1:0.1f} GW'.format(surplus_bar[i]*100/(gamma_bar[i]*mean(L))+1e-10,gamma_bar[i]*mean(L)),
        #        ha='center', va='center')
        gca().text(bars[i].get_x()+bars[i].get_width()/2., height, '{0:0.1f} TWh'.format(abs(surplus_bar[i] + 1e-6)),
                ha='center', va='center')
    
    
    #autolabel_bars(bars,gca(),heighttxt='{0:0.1f}%\nof\n'+'{0:0.1f} GW'.format(),scale_factors=100/(gamma_bar*mean(L)))

    axis(xmin=amin(gamma_bar*100)-10,xmax=amax(gamma_bar*100)+10,ymin=0,ymax=1.05*amax(surplus_bar))
    xticks(gamma_bar*100,['{0:.0f}%'.format(x) for x in gamma_bar*100],va='top')
    xlabel(r'Target share of total electricity demand'.format(mean(L)))
    ylabel('Surplus [TWh]')

    ax=gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    tight_layout()
    save_file_name = 'plot_surplus_bar_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)
    
    if returnall:
        return gamma, surplus, gamma_bar, surplus_bar


################################################################################################
################################################################################################



#   Wind vs wind-solar mix
#   plot_surplus_bar_comp(alpha_w=[1,None],CS=[None,None],color_bar=[color_wind,color_opt],hatch_bar=[False,False],textlabels=['Wind only','Balancing optimal mix'])
#
#   With and without storage
#   plot_surplus_bar_comp(ISO='DK',gamma_bar=array([.25,.50,.75,1.]),alpha_w=[1,None,1,None],CS=[None,None,12,12],N_gamma=21,bw = .10,color_bar=[color_wind,color_opt,color_wind,color_opt],hatch_bar=[False,False,True,True],textlabels=['Wind only','Bal. opt. mix','Wind only + storage','Bal. opt. mix + storage'])
#
def plot_surplus_bar_comp(ISO='DK',gamma_bar=array([.25,.50,.75,1.]),alpha_w=[1,None,1,None],CS=[None,None,12,12],N_gamma=21,bw = .10,color_bar=[color_wind,color_opt,color_wind,color_opt],hatch_bar=[False,False,True,True],textlabels=['Wind only','Balancing optimal mix','Wind only, 50 GWh storage','Balancing optimal mix, 50 GWh storage']):


    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    annual_TWh = mean(L)*365*24/1e3

    gamma_, surplus_, gamma_bar_, surplus_bar_ = zeros((len(alpha_w),N_gamma)), zeros((len(alpha_w),N_gamma)), zeros((len(alpha_w),len(gamma_bar))), zeros((len(alpha_w),len(gamma_bar)))
    for i in arange(len(alpha_w)):
    
        gamma_[i], surplus_[i], gamma_bar_[i], surplus_bar_[i] = plot_surplus_bar(ISO,gamma_bar,alpha_w[i],CS[i],N_gamma,returnall=True)
    
        print diff(surplus_[i])/((gamma_[i][1]-gamma_[i][0])*annual_TWh)

    close(1); figure(1); clf()

    fill_between(gamma_[i]*annual_TWh,amax(surplus_,axis=0),amin(surplus_,axis=0),color=(.7,.7,.7))

    plot_vertical_line_and_label(.25*annual_TWh,interp(25,gamma_[0]*100,amax(surplus_,axis=0))+.5,textlabel=r'25%',dx=.5,ls='--',color='k',lw=1)
    plot_vertical_line_and_label(.50*annual_TWh,interp(50,gamma_[0]*100,amax(surplus_,axis=0))+.5,textlabel=r'50%',dx=.5,ls='--',color='k',lw=1)
    plot_vertical_line_and_label(.75*annual_TWh,interp(75,gamma_[0]*100,amax(surplus_,axis=0))+.5,textlabel=r'75%',dx=.5,ls='--',color='k',lw=1)
    plot_vertical_line_and_label(1.00*annual_TWh,interp(100,gamma_[0]*100,amax(surplus_,axis=0))+.5,textlabel=r'100%',dx=.5,ls='--',color='k',lw=1)

    for i in arange(len(alpha_w)):

        #plot(gamma_[i]*100,surplus_[i],'k-',lw=1.5,zorder=0)
        
        if hatch_bar[i]:
            bars = bar(annual_TWh*(gamma_bar_[i] + 1.20*bw/len(alpha_w)*(i-.5*len(alpha_w)+.5)+0.1*bw/len(alpha_w)),surplus_bar_[i],align='center',width=bw*annual_TWh/len(alpha_w),color='k', fc=color_bar[i], hatch="//",label=textlabels[i])
        else:
            bars = bar(annual_TWh*(gamma_bar_[i] + 1.20*bw/len(alpha_w)*(i-.5*len(alpha_w)+.5)+0.1*bw/len(alpha_w)),surplus_bar_[i],align='center',width=bw*annual_TWh/len(alpha_w),color=color_bar[i],label=textlabels[i])

    

    #axis(xmin=0,xmax=amax(gamma_bar*annual_TWh)+.10*annual_TWh,ymin=0,ymax=1.0*amax(surplus_bar_)+1)
    axis(xmin=0,xmax=amax(gamma_bar*annual_TWh)+.10*annual_TWh,ymin=0,ymax=12.5)
    
    #xticks(gamma_bar*100,['{0:.0f}%'.format(x) for x in gamma_bar*100],va='top')
    xticks(concatenate([[0],gamma_bar])*annual_TWh,['{0:.1f}'.format(x) for x in concatenate([[0],gamma_bar])*annual_TWh],va='top')
    xlabel(r'Wind plus solar energy [TWh/yr]')
    ylabel('VRE surplus [TWh/yr]')

    ax=gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    leg = legend(loc='upper left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    tight_layout(pad=0.2)
    save_file_name = 'plot_surplus_bar_comp_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)



################################################################################################
################################################################################################



def autolabel_bars(rects,ax,heighttxt=r'{0:0.1f}',scale_factors=1):
    # attach some text labels
    for i in arange(len(rects)):
        rect = rects[i]
        if size(scale_factors)>1:
            scale_factor = scale_factors[i]
        else:
            scale_factor = scale_factors
            
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., height/2., heighttxt.format(height*scale_factor),
                ha='center', va='center',backgroundcolor='w')


################################################################################################
################################################################################################



def get_compare_VRES_load(ISO='DK', gamma=0.5, alpha_w=.5, CS=None, silent=True):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

    if alpha_w==None:
        alpha_w = get_optimal_mix_balancing(L, Gw, Gs, gamma)
        if ~silent:
            print 'Mix not specified. Optimal mix alpha_w={0} chosen'.format(alpha_w)
                    
    wind = gamma*alpha_w*Gw*mean(L)
    solar = gamma*(1-alpha_w)*Gs*mean(L)
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

    return L, wind_local, solar_local, curtailment, filling, extraction



################################################################################################
################################################################################################



# Jan. 2000:
# plot_hourly_generation(alpha_w=None,date_start=datestr2num('1-1-2000'),monday_offset=5,titletxt='Denmark, Jan. 2000',label='Jan2000_optimal')
# plot_hourly_generation(alpha_w=1,date_start=datestr2num('1-1-2000'),monday_offset=5,titletxt='Denmark, Jan. 2000',label='Jan2000_wind')
# plot_hourly_generation(alpha_w=0.,date_start=datestr2num('1-1-2000'),monday_offset=5,titletxt='Denmark, Jan. 2000',label='Jan2000_solar')
#
# plot_hourly_generation(alpha_w=1.,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Denmark, Jan. 2000',label='Jan2000_day_wind')
# plot_hourly_generation(alpha_w=None,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Denmark, Jan. 2000',label='Jan2000_day_optimal')
#
# July 2000:
# plot_hourly_generation(alpha_w=None,date_start=datestr2num('7-1-2000'),monday_offset=5,titletxt='Denmark, July 2000',label='July2000_optimal')
# plot_hourly_generation(alpha_w=1,date_start=datestr2num('7-1-2000'),monday_offset=5,titletxt='Denmark, July 2000',label='July2000_wind')
# plot_hourly_generation(alpha_w=0.,date_start=datestr2num('7-1-2000'),monday_offset=5,titletxt='Denmark, July 2000',label='July2000_solar')
#
def plot_hourly_generation(ISO='DK', gamma=0.5, alpha_w=.5, CS=None, date_start=datestr2num('1-1-2000'), N_days=30, monday_offset=5, titletxt='Denmark, Jan. 2000', label='TestFigure'):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    mask = find((t+datetime_offset>=date_start)*(t+datetime_offset<=date_start+N_days))

    if alpha_w==None:
        alpha_w = get_optimal_mix_balancing(L, Gw, Gs, gamma)
        print 'Mix not specified. Optimal mix alpha_w={0} chosen'.format(alpha_w)
                    
    wind = gamma*alpha_w*Gw*mean(L)
    solar = gamma*(1-alpha_w)*Gs*mean(L)
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

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

    fill_between(t[mask],wind_local[mask],color=color_wind,edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask],wind_local[mask],color=color_solar,edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+filling[mask],wind_local[mask]+solar_local[mask],color='b',edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+extraction[mask],wind_local[mask]+solar_local[mask],color='m',edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+filling[mask]+curtailment[mask],wind_local[mask]+solar_local[mask]+filling[mask],color='r',edgecolor=color_edge,lw=1)
    pp_wind = Rectangle((0, 0), 1, 1, facecolor=color_wind)
    pp_solar = Rectangle((0, 0), 1, 1, facecolor=color_solar)
    pp_curtailment = Rectangle((0, 0), 1, 1, facecolor='r')

    pp_load = plot(t[mask],L[mask],color='k',lw=1.5)

    axis(xmin=t[mask[0]],xmax=t[mask[-1]],ymin=0,ymax=2.05*mean(L))

    ylabel('Power [GW]')

    day_names = array([calendar.day_abbr[mod(i,7)] for i in arange(N_days)+monday_offset])
    day_names[find([d!='Mon' for d in day_names])] = ''
    xticks(t[mask[0]]+arange(N_days),day_names,rotation=-45,ha='left')

    pp = [pp_load[0],pp_wind,pp_solar,pp_curtailment]
    txtlabels = ['Load ({0})'.format(ISO),'Wind','Solar','Surplus']

    titletxt = titletxt + '\nwind/solar mix: {0}/{1}'.format(int(round(100*alpha_w)),int(round(100*(1-alpha_w))))
    leg = legend(pp,txtlabels,loc='upper left',ncol=4,title=titletxt);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    ax=gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    tight_layout()
    save_file_name = 'plot_hourly_generation_'+ISO+'_'+label+'.pdf'
    save_figure(save_file_name)


################################################################################################
################################################################################################



# plot_hourly_generation_alt(alpha_w=1.,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, wind only',label='week_wind')
# plot_hourly_generation_alt(alpha_w=None,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, optimal wind and solar mix',label='week_optimal')
#
# plot_hourly_generation_alt(alpha_w=1.,date_start=datestr2num('3-6-2000'),N_days=2,monday_offset=7,titletxt='Spring week, wind only',label='week_wind')
#
# plot_hourly_generation_alt(alpha_w=1.,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, wind only',label='week_wind_storage',CS=7.4)
# plot_hourly_generation_alt(alpha_w=None,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, optimal wind and solar mix',label='week_optimal_storage',CS=7.4)

### 2013
#
#   plot_hourly_generation_alt(alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_no_storage')
#
#   plot_hourly_generation_alt(alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_storage_no_gain',CS=2.54)
#
#   plot_hourly_generation_alt(alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_storage_gain_100',gain_storage=.99999,CS=2.54)
#
#   plot_hourly_generation_alt(alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_storage_gain_90',gain_storage=.9,CS=2.54)
#
#   plot_hourly_generation_alt(alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_storage_gain_99',gain_storage=.99,CS=2.54)

### DK paper
#
#   plot_hourly_generation_alt(alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_no_storage')
#
#
#   plot_hourly_generation_alt(alpha_w=1,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_wind_only_no_storage')
#
#   plot_hourly_generation_alt(ISO='DK-hist',alpha_w=1,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_wind_only_hist')
#
#
#   plot_hourly_generation_alt(ISO='DK-hist',alpha_w=None,gamma=.75,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_optimal_hist')
#
#   plot_hourly_generation_alt(ISO='DK-hist',alpha_w=1,gamma=.4,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_wind_only_hist_40p')
#
#   plot_hourly_generation_alt(ISO='DK-hist',alpha_w=1,gamma=.25,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt=None,label='week_wind_only_hist_25p')
#
def plot_hourly_generation_alt(ISO='DK', gamma=0.5, alpha_w=.5, CS=None, date_start=datestr2num('3-6-2000'), N_days=7, monday_offset=7, titletxt='Spring week', label='TestFigure', P_in=None, P_out=None, gain_storage=None, eta_in = 1., eta_out = 1., fancy_storage=False):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    mask = find((t+datetime_offset>=date_start)*(t+datetime_offset<=date_start+N_days))

    if alpha_w==None:
        alpha_w = get_optimal_mix_balancing(L, Gw, Gs, gamma)
        print 'Mix not specified. Optimal mix alpha_w={0} chosen'.format(alpha_w)
                    
    wind = gamma*alpha_w*Gw*mean(L)
    solar = gamma*(1-alpha_w)*Gs*mean(L)
    mismatch = (wind+solar) - L
    
    if CS!=None or CS==0:
        
        if gain_storage!=None:
            Storage_benefit, P_in_, P_out_ = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w, CS, 1-gain_storage,eta_charge=eta_in,eta_discharge=eta_out)
            P_in, P_out = P_in_[0][0], P_out_[0][0]
            print Storage_benefit, P_in, P_out
        
        mismatch_r = get_policy_2_storage_modified(mismatch/mean(L), eta_in = eta_in, eta_out = eta_out, CS = CS, P_in=P_in, P_out=P_out,fancy_storage=fancy_storage)[0]*mean(L)
        #mismatch_r = get_policy_2_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = CS)[0]
        
    else:
        mismatch_r = mismatch
    
    curtailment = get_positive(mismatch_r)
    filling = get_positive(mismatch - mismatch_r)
    extraction = get_positive(-(mismatch - mismatch_r))

    wind_local = wind - (curtailment+filling)*wind/(wind+solar+1e-10)
    solar_local = solar - (curtailment+filling)*solar/(wind+solar+1e-10)

    print 'sum(curtailment), sum(filling), sum(extraction), sum(wind_local), sum(solar_local)'
    print sum(curtailment), sum(filling), sum(extraction), sum(wind_local), sum(solar_local)

    plt.close(1); plt.figure(1); plt.clf()

    fill_between(t[mask],wind_local[mask],color=color_wind,edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask],wind_local[mask],color=color_solar,edgecolor=color_edge,lw=1)
    p = fill_between(t[mask],wind_local[mask]+solar_local[mask]+filling[mask],wind_local[mask]+solar_local[mask],color='g',edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+extraction[mask],wind_local[mask]+solar_local[mask],color='g',edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+filling[mask]+curtailment[mask],wind_local[mask]+solar_local[mask]+filling[mask],color='r',edgecolor=color_edge,lw=1)
    
    p.set_facecolors("none")

    plt.plot(t[mask],np.amax(L)*np.ones_like(t[mask]),'k--')

    from matplotlib.patches import PathPatch
    for path in p.get_paths():
        p1 = PathPatch(path, fc="none", hatch="/")
        gca().add_patch(p1)
        p1.set_zorder(p.get_zorder()-0.1)
    
    
    
    pp_wind = Rectangle((0, 0), 1, 1, facecolor=color_wind)
    pp_solar = Rectangle((0, 0), 1, 1, facecolor=color_solar)
    pp_curtailment = Rectangle((0, 0), 1, 1, facecolor='r')
    pp_storage = Rectangle((0, 0), 1, 1, facecolor='g')
    pp_filling = Rectangle((0, 0), 1, 1, facecolor='none', hatch="/")

    pp_load = plot(t[mask],L[mask],color='k',lw=1.5)

    axis(xmin=t[mask[0]],xmax=t[mask[-1]],ymin=0,ymax=9.1)#1.9*mean(L))

    ylabel('Power [GW]')

    day_names = array([calendar.day_abbr[mod(i,7)] for i in arange(N_days)+monday_offset])
    #day_names[find([d!='Mon' for d in day_names])] = ''
    xticks(t[mask[0]]+arange(N_days),day_names,rotation=-45,ha='left')

    if CS==None:
        pp = [pp_load[0],pp_wind,pp_solar,pp_curtailment]
        txtlabels = ['Load'.format(ISO),'Wind','Solar','VRE-surplus']
        leg = legend(pp,txtlabels,loc='upper left',ncol=3,title=titletxt);
    else:
        pp = [pp_load[0],pp_filling,pp_wind,pp_storage,pp_solar,pp_curtailment]
        txtlabels = ['Load'.format(ISO),'Charging','Wind','Discharging','Solar','Rem. surplus']    
        leg = legend(pp,txtlabels,loc='upper left',ncol=3,title=titletxt);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    tight_layout(pad=0.2)
    save_file_name = 'plot_hourly_generation_alt_'+ISO+'_'+label+'.pdf'
    save_figure(save_file_name)


################################################################################################
################################################################################################



# plot_monthly_summary('DK',gamma=0.75,alpha_w=None,label='optimal')
# plot_monthly_summary('DK',gamma=0.75,alpha_w=1,label='wind')
#
# plot_monthly_summary('DK',gamma=0.5,alpha_w=None,label='optimal')
# plot_monthly_summary('DK',gamma=0.5,alpha_w=1,label='wind')
# plot_monthly_summary('DK',gamma=0.5,alpha_w=0,label='solar')
#
# plot_monthly_summary('DK',gamma=0.5,alpha_w=None,label='optimal_storage',CS=2.54)
# plot_monthly_summary('DK',gamma=0.5,alpha_w=1,label='wind_storage',CS=2.54)
# plot_monthly_summary('DK',gamma=0.5,alpha_w=0,label='solar_storage',CS=2.54)
# plot_monthly_summary('DK',gamma=0.75,alpha_w=None,label='optimal_storage',CS=2.54)
# plot_monthly_summary('DK',gamma=0.75,alpha_w=1,label='wind_storage',CS=2.54)
#
# plot_monthly_summary('DK',gamma=1,alpha_w=None,label='100p_optimal')
# plot_monthly_summary('DK',gamma=1,alpha_w=1,label='100p_wind')
def plot_monthly_summary(ISO='DK', gamma=.5, alpha_w=None, CS=None, titletxt='Denmark, 2000-2007',label='TestFigure'):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    

    if alpha_w==None:
        alpha_w = get_optimal_mix_balancing(L, Gw, Gs, gamma)[0]
        optimal = True
        print 'Mix not specified. Optimal mix alpha_w={0} chosen'.format(alpha_w)
    else:
        optimal = False
    
    wind = gamma*alpha_w*Gw*mean(L)
    solar = gamma*(1-alpha_w)*Gs*mean(L)
    mismatch = (wind+solar) - L
    
    if CS!=None or CS==0:
        mismatch_r = get_policy_2_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = CS*mean(L))[0]
    else:
        mismatch_r = mismatch
    
    curtailment = get_positive(mismatch_r)
    filling = get_positive(mismatch - mismatch_r)
    extraction = get_positive(-(mismatch - mismatch_r))
    
    print "Curtailment: ", sum(curtailment)/8.
    print "Filling: ", sum(filling)/8.
    print "Extraction: ", sum(extraction)/8.
    
    wind_local = wind - (curtailment+filling)*wind/(wind+solar+1e-10)
    solar_local = solar - (curtailment+filling)*solar/(wind+solar+1e-10)
    
    t_month = array([d.month for d in num2date(t+datetime_offset)])
    
    months = arange(0,12)+1
    load_sum, wind_sum, solar_sum, storage_sum, curtail_sum, weight = zeros(12), zeros(12), zeros(12), zeros(12), zeros(12), zeros(12)
    for month in months:
        load_sum[month-1] = sum(L[find(t_month==month)])
        wind_sum[month-1] = sum(wind_local[find(t_month==month)])
        solar_sum[month-1] = sum(solar_local[find(t_month==month)])
        storage_sum[month-1] = sum(extraction[find(t_month==month)])
        curtail_sum[month-1] = sum(curtailment[find(t_month==month)])
        weight[month-1] = len(find(t_month==month))/float(len(t))

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

    #Monthly values
    pp_wind = bar(months,wind_sum*100/load_sum,color=color_wind)
    pp_solar = bar(months,solar_sum*100/load_sum,bottom=wind_sum*100/load_sum,color=color_solar)
    pp_storage = bar(months,storage_sum*100/load_sum,bottom=(wind_sum+solar_sum)*100/load_sum,color='g')
    pp_curtailment = bar(months,curtail_sum*100/load_sum,bottom=(wind_sum+solar_sum+storage_sum)*100/load_sum,color='r')
    
    pp_gross = axhline(100*gamma,color='k',ls='--',lw=1.5)
    
    #Average values
    bar(13.5,sum(wind_sum)*100/sum(load_sum),color=color_wind)
    bar(13.5,sum(solar_sum)*100/sum(load_sum),bottom=sum(wind_sum)*100/sum(load_sum),color=color_solar)
    bar(13.5,sum(storage_sum)*100/sum(load_sum),bottom=sum(wind_sum+solar_sum)*100/sum(load_sum),color='g')
    bar(13.5,sum(curtail_sum)*100/sum(load_sum),bottom=sum(wind_sum+solar_sum+storage_sum)*100/sum(load_sum),color='r')
    
    print 'Local excess: {0} pp'.format(sum(curtail_sum)*100/sum(load_sum),bottom=sum(wind_sum+solar_sum)*100/sum(load_sum))
    
    axis(xmin=.5,xmax=14.75,ymin=0,ymax=100)
    
    month_names = array([calendar.month_abbr[m] for m in months])
    xticks(concatenate([months,[13.5]]),concatenate([month_names,['Av.']]),rotation=-45,ha='left')
    
    ylabel('Percentage points (pp)')
    
    pp = [pp_gross,pp_wind[0],pp_solar[0],pp_curtailment[0]]
    txtlabels = ['Gross share','Wind','Solar','Local excess']
    titletxt = titletxt + '\nwind/solar mix: {0}/{1}'.format(int(round(100*alpha_w)),int(round(100*(1-alpha_w))))
    leg = legend(pp,txtlabels,loc='upper left',ncol=4,title=titletxt);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    tight_layout()
    save_file_name = 'plot_monthly_summary_'+ISO+'_'+label+'.pdf'
    save_figure(save_file_name)


    ### Single bar summary:
    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([3,2])

    ax=axes([0.03,0.03,.33,1])

    #Average values
    bar(0,mean(wind_sum)/mean(load_sum),color=color_wind,align='center',width=1)
    bar(0,mean(solar_sum)/mean(load_sum),bottom=mean(wind_sum)/mean(load_sum),color=color_solar,align='center',width=1)
    bar(0,mean(storage_sum)/mean(load_sum),bottom=mean(wind_sum+solar_sum)/mean(load_sum),color='g',align='center',width=1)
    bar(0,mean(curtail_sum)/mean(load_sum),bottom=mean(wind_sum+solar_sum+storage_sum)/mean(load_sum),color='r',align='center',width=1)

    axis(xmin=-.75,xmax=.75,ymin=0,ymax=1.1)

    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.tick_params(labelbottom=0, labeltop=0, labelleft=0, labelright=0)
    
    if CS==None:
        text(0.37,0.9,'Total: {0:0.1f} GWh\nSurplus: {1:.1f}%'.format(gamma*mean(L),mean(curtail_sum)/(mean(load_sum)*gamma)*100),fontsize=10,weight='semibold',ha='left',va='top',transform = gcf().transFigure,bbox=dict(boxstyle="round, pad=.75", fc="w",lw=1.5))
    else:
        text(0.37,0.9,'Total: {0:0.1f} GWh\nSurplus:      {1:.1f}%\nStorage:     -{2:.1f}%\n----------------------\nRemainder: {3:.1f}%\n\nStorage volume:\n{4:.1f} GWh'.format(gamma*mean(L),mean(curtail_sum+storage_sum)/(mean(load_sum)*gamma)*100,mean(storage_sum)/(mean(load_sum)*gamma)*100,mean(curtail_sum)/(mean(load_sum)*gamma)*100,mean(L)*CS),fontsize=10,weight='semibold',ha='left',va='top',transform = gcf().transFigure,bbox=dict(boxstyle="round, pad=.75", fc="w",lw=1.5))


    #tight_layout()
    #save_file_name = 'plot_single_bar_summary_'+ISO+'_'+label+'.pdf'
    #save_figure(save_file_name)

    ### Single bar summary (alternative):
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='sans-serif')
    matplotlib.rcParams['text.latex.preamble'] = '\usepackage{color}'
    
    
    
    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([3,2])

    ax=axes([0.03,0.03,.33,1])

    #Average values
    bar(0,mean(wind_sum)/mean(load_sum),color=color_wind,align='center',width=1)
    bar(0,mean(solar_sum)/mean(load_sum),bottom=mean(wind_sum)/mean(load_sum),color=color_solar,align='center',width=1)
    bar(0,mean(storage_sum)/mean(load_sum),bottom=mean(wind_sum+solar_sum)/mean(load_sum),color='g',align='center',width=1)
    bar(0,mean(curtail_sum)/mean(load_sum),bottom=mean(wind_sum+solar_sum+storage_sum)/mean(load_sum),color='r',align='center',width=1)

    axis(xmin=-.75,xmax=.75,ymin=0,ymax=1.1)

    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.tick_params(labelbottom=0, labeltop=0, labelleft=0, labelright=0)
    

    if optimal == True:
       tabletitle = 'Optimal mix' 
    elif alpha_w==1:
        tabletitle = r'Wind-only'
    elif alpha_w==0:
        tabletitle = r'Solar-only'       
    else:
       tabletitle = '{0:0.0f}\% wind, {1:0.0f}\% solar'.format(100*alpha_w,(1-alpha_w)*100)
    
    print tabletitle
    
    if CS==None:
       
        #tabletext = r'\begin{tabular}{ l r } \multicolumn{2}{c}{\bf RES-E budget (DK)} \\[.5ex]' + \
        
        tabletext = r'\begin{tabular}{ l r } \multicolumn{2}{c}{\bf ' + tabletitle + r'} \\[.5ex]' + \
            r'Total VRE & ' + '{0:0.1f}'.format(gamma*mean(L)*24*365/1e3) + r'\phantom{)} TWh \\ ' + \
            r'Consumed & (' + '{0:0.1f}'.format(sum(wind_sum+solar_sum)/len(L)*24*365/1e3) + r') TWh  \\ \hline ' + \
            r'Surplus &  ' + '{0:0.1f}'.format(sum(curtail_sum)/len(L)*24*365/1e3) + r'\phantom{)} TWh \\[1ex]' + \
            r'\multicolumn{2}{l}{\it No storage.} \\' + \
            r'\end{tabular}'
        text(0.37,0.9,tabletext,fontsize=10,weight='semibold',ha='left',va='top',transform = gcf().transFigure,bbox=dict(boxstyle="round, pad=.75", fc="w",lw=1.5))
        
    else:
        
        tabletext = r'\begin{tabular}{ l r } \multicolumn{2}{c}{\bf ' + tabletitle + r'} \\[.5ex]' + \
            r'Total VRE & ' + '{0:0.1f}'.format(gamma*mean(L)*24*365/1e3) + r'\phantom{)} TWh \\ ' + \
            r'Consumed & (' + '{0:0.1f}'.format(sum(wind_sum+solar_sum)/len(L)*24*365/1e3) + r') TWh  \\ ' + \
            r'Storage & (' + '{0:0.1f}'.format(sum(storage_sum)/len(L)*24*365/1e3) + r') TWh  \\ \hline ' + \
            r'Surplus &  ' + '{0:0.1f}'.format(sum(curtail_sum)/len(L)*24*365/1e3) + r'\phantom{)} TWh \\[1ex]' + \
            r'\multicolumn{2}{l}{\it Storage volume:} \\' + \
            r'\multicolumn{2}{l}{'+'{0:0.0f}'.format(mean(L)*CS)+r' GWh} \\' + \
            r'\end{tabular}'


        text(0.37,0.9,tabletext,fontsize=10,weight='semibold',ha='left',va='top',transform = gcf().transFigure,bbox=dict(boxstyle="round, pad=.75", fc="w",lw=1.5))       
        
        #r'Storage & (' + '{0:0.1f}'.format(sum(storage_sum)/len(L)*24*365/1e3) + r') TWh  \\ \hline ' + \
        
        #text(0.37,0.9,'Total: {0:0.1f} TWh\nSurplus:      {1:.1f} TWh\nStorage:     -{2:.1f} TWh\n----------------------\nRemainder: {3:.1f}%\n\nStorage volume:\n{4:.1f} GWh'.format(gamma*mean(L),mean(curtail_sum+storage_sum)/(mean(load_sum)*gamma)*100,mean(storage_sum)/(mean(load_sum)*gamma)*100,mean(curtail_sum)/(mean(load_sum)*gamma)*100,mean(L)*CS),fontsize=10,weight='semibold',ha='left',va='top',transform = gcf().transFigure,bbox=dict(boxstyle="round, pad=.75", fc="w",lw=1.5))


    #tight_layout()
    save_file_name = 'plot_single_bar_summary_alt_'+ISO+'_{0:0.0f}p_'.format(gamma*100)+label+'.pdf'
    save_figure(save_file_name)


    ### Single bar summary legend
    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,1])
    
    axis('off')
    
    pp_wind = Rectangle((0, 0), 1, 1, facecolor=color_wind)

    pp_solar = Rectangle((0, 0), 1, 1, facecolor=color_solar)
    pp_curtailment = Rectangle((0, 0), 1, 1, facecolor='r')
    pp_storage = Rectangle((0, 0), 1, 1, facecolor='g')    
    
    pp = [pp_wind,pp_solar,pp_storage,pp_curtailment]
    txtlabels = ['Wind','Solar','Storage output','Surplus - Storage output']
    leg = legend(pp,txtlabels,loc='upper left',ncol=2,frameon=False);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    save_file_name = 'plot_single_bar_summary_legend.pdf'
    save_figure(save_file_name)
    
    
 
################################################################################################
################################################################################################

   

###############
### Tools: ####
###############

def get_optimal_mix_storage(L, GW, GS, gamma=1., p_interval=0.01, returnall=False):

	L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

	l = weighed_sum(L)
	Gw = weighed_sum(GW)	
	Gs = weighed_sum(GS)

	mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
	E_H = lambda alpha_w: get_policy_2_storage(mismatch(alpha_w))[1]
	
	alpha_w_opt = fmin(E_H,0.5,disp=True)
	
	if returnall:
	
		balancing_fixed_E_H = lambda alpha_w: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w),storage_capacity = E_H(alpha_w_opt))[0])) - p_interval*sum(l)
	
		if balancing_fixed_E_H(0)>0:
			if balancing_fixed_E_H(1)>0:
				alpha_w_opt_1p_interval = array([brentq(balancing_fixed_E_H, 0, alpha_w_opt),brentq(balancing_fixed_E_H, alpha_w_opt, 1)])
			else:
				alpha_w_opt_1p_interval = array([brentq(balancing_fixed_E_H, 0, alpha_w_opt),1.])
		else:	
			if balancing_fixed_E_H(1)>0:
				alpha_w_opt_1p_interval = array([0.,brentq(balancing_fixed_E_H, alpha_w_opt, 1)])
			else:
				alpha_w_opt_1p_interval = array([0.,1.])
	
		print balancing_fixed_E_H(0), balancing_fixed_E_H(1)
		print alpha_w_opt, alpha_w_opt_1p_interval
	
	
		#Returns: alpha_w_opt, alpha_w_opt_1p_interval
		return alpha_w_opt, alpha_w_opt_1p_interval, E_H(alpha_w_opt)
	else:
		return E_H(alpha_w_opt), alpha_w_opt



################################################################################################
################################################################################################



def get_balancing(L, GW, GS, gamma=1, alpha=1., CS=None, returnall=False, eta_in=1., eta_out=1.):

	L, GW, GS, gamma, alpha = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2), array(gamma,ndmin=1), array(alpha,ndmin=1)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

	l = weighed_sum(L)
	Gw = weighed_sum(GW)	
	Gs = weighed_sum(GS)
	
	mismatch = lambda alpha_w, gamma: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
	
	if CS==None:
		res_load_sum = lambda alpha_w, gamma: sum(get_positive(-mismatch(alpha_w,gamma)))
	else:
		res_load_sum = lambda alpha_w, gamma: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w,gamma),eta_in,eta_out,storage_capacity = CS)[0]))
	
	Gamma, Alpha = meshgrid(gamma,alpha)
	
	Res_load_sum = zeros(Gamma.shape)
	for i in arange(size(Gamma)):
		Res_load_sum.flat[i] = res_load_sum(Alpha.flat[i],Gamma.flat[i])
	
	if returnall:
		return Res_load_sum, Gamma, Alpha
	else:
		return Res_load_sum


################################################################################################
################################################################################################



def get_mismatch(L, GW, GS, gamma=1, alpha_w=1.,CS=None):

    L, GW, GS, gamma, alpha_w = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2), array(gamma,ndmin=1), array(alpha_w,ndmin=1)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.

    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w, gamma: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l

    if CS==None:
        return mismatch(alpha_w,gamma)
    else:
        return get_policy_2_storage(mismatch(alpha_w,gamma),storage_capacity = CS)[0]


################################################################################################
################################################################################################



def get_optimal_path_balancing(L, GW, GS, gamma=linspace(0,1,5), p_interval=0.01, CS=None, returnall=False, normalized=True, eta_in=1, eta_out=1):
    """Wraper for get_optimal_mix_balancing(). This function allows gamma to be an array."""

    gamma = array(gamma,ndmin=1)
    
    if returnall==True:
        
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = get_optimal_mix_balancing(L, GW, GS, gamma[0], p_interval, CS, returnall, normalized, eta_in=eta_in, eta_out=eta_out)
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = expand2array(alpha_w_opt,gamma), expand2array(alpha_w_opt_1p_interval,gamma), expand2array(res_load_sum_opt,gamma), expand2array(mismatch_opt,gamma), expand2array(res_load_sum_1p,gamma)
        
        for i in arange(1,len(gamma)):
            alpha_w_opt[i], alpha_w_opt_1p_interval[i], res_load_sum_opt[i], mismatch_opt[i], res_load_sum_1p[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval, CS, returnall, normalized, eta_in=eta_in, eta_out=eta_out)
        
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p

    else:
        alpha_w_opt = zeros(gamma.shape)
        for i in arange(len(gamma)):
            alpha_w_opt[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval, CS, returnall, normalized, eta_in=eta_in, eta_out=eta_out)
        
        return alpha_w_opt


################################################################################################
################################################################################################


        
def get_optimal_mix_balancing(L, GW, GS, gamma=1., p_interval=0.01, CS=None, returnall=False, normalized=True, DefaultWind=True, eta_in=1., eta_out=1.):
    """ Old version of get_optimal_mix()"""

    if returnall:
        ## The returned value of CS is ommitted for compatability reasons.
        return get_optimal_mix(L, GW, GS, gamma, p_interval, CS, returnall, normalized, DefaultWind, eta_in, eta_out)[:-1]
    else:
        return get_optimal_mix(L, GW, GS, gamma, p_interval, CS, returnall, normalized, DefaultWind, eta_in, eta_out)


################################################################################################
################################################################################################



def get_optimal_mix(L, GW, GS, gamma=1., p_interval=0.01, CS=None, returnall=False, normalized=True, DefaultWind=True, eta_in=1., eta_out=1.):

    L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l

    if CS==None:
        ## Balancing optimal mix
        res_load_sum = lambda alpha_w: sum(get_positive(-mismatch(alpha_w))) - alpha_w*0.001*sign(DefaultWind-.5)
        alpha_w_opt = fmin(res_load_sum,0.5,disp=False)
        CS = 0.
        
    elif isnan(CS):
        ## Seasonal optimal mix: In this case, the seasonal optimal mix is found by minimizing CS instead of the summed residual load.
        CS_ = lambda alpha_w: get_policy_2_storage(mismatch(alpha_w),eta_in,eta_out)[1]
        alpha_w_opt = fmin(CS_,0.5,disp=False)
        CS = CS_(alpha_w_opt)
        res_load_sum = lambda alpha_w: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w),eta_in,eta_out,storage_capacity = CS)[0])) - alpha_w*0.001*sign(DefaultWind-.5)
    
    else:
        ## Optimal mix given a specific storage energy capacity.
        res_load_sum = lambda alpha_w: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w),eta_in,eta_out,storage_capacity = CS)[0])) - alpha_w*0.001*sign(DefaultWind-.5)
        alpha_w_opt = fmin(res_load_sum,0.5,disp=False)

    ## Correct for nonsense alpha_w since fmin does not include bounds.
    if alpha_w_opt>1.:
        alpha_w_opt = 1.
    elif alpha_w_opt<0.:
        alpha_w_opt = 0.

    if normalized:
        if CS==None:
            mismatch_opt = mismatch(alpha_w_opt)
        else:
            mismatch_opt = get_policy_2_storage(mismatch(alpha_w_opt),eta_in,eta_out,storage_capacity = CS)[0]
    else:
        if CS==None:
            mismatch_opt = mismatch(alpha_w_opt)*mean(sum(L,axis=0))
        else:
            mismatch_opt = get_policy_2_storage(mismatch(alpha_w_opt),eta_in,eta_out,storage_capacity = CS)[0]*mean(sum(L,axis=0))
            
    res_load_sum_opt = res_load_sum(alpha_w_opt)

    if returnall:
        res_load_sum_1p_interval = lambda alpha_w: res_load_sum(alpha_w)-(res_load_sum(alpha_w_opt)+p_interval*sum(l))
        
        if sign(res_load_sum_1p_interval(0))!=sign(res_load_sum_1p_interval(alpha_w_opt)):
            lower_bound = brentq(res_load_sum_1p_interval, 0, alpha_w_opt)
        else:
            lower_bound = 0
        
        if sign(res_load_sum_1p_interval(1))!=sign(res_load_sum_1p_interval(alpha_w_opt)):
            upper_bound = brentq(res_load_sum_1p_interval, alpha_w_opt, 1)
        else:
            upper_bound = 1
        
        alpha_w_opt_1p_interval = array([lower_bound,upper_bound])
        res_load_sum_1p = amax([res_load_sum(lower_bound),res_load_sum(upper_bound)])
        
        #Returns: alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p, CS
    else:
        return alpha_w_opt


################################################################################################
################################################################################################



#
# get_min_storage_cap(L, Gw, Gs, gamma=1, alpha_w=1., CS=6)
#
def get_min_storage_cap(L, GW, GS, gamma=1, alpha_w=1., CS=None, acc=1e-4):
    """Finds upper and lower limits of the storage in and output capacities. acc sets the relative increase in balancing energy that is acceptable.
    WARNING! This function is slow to evaluate."""


    L, GW, GS, gamma, alpha_w = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2), array(gamma,ndmin=1), array(alpha_w,ndmin=1)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    if len(alpha_w)==1:
        alpha_w = alpha_w*ones_like(gamma)
    
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w, gamma: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
    mismatch_ = lambda alpha_w, gamma, P_in, P_out: amax([amin([mismatch(alpha_w, gamma),P_in*ones_like(l)],axis=0),-P_out*ones_like(l)],axis=0)

    res_load_sum = lambda alpha_w, gamma, P_in, P_out: sum(get_positive(-get_policy_2_storage(mismatch_(alpha_w,gamma,P_in,P_out),storage_capacity = CS)[0]) + get_positive(-mismatch(alpha_w, gamma)-P_out))
    
    Res_load_sum = zeros_like(gamma)
    P_in_ = zeros((len(gamma),2))
    P_out_ = zeros((len(gamma),2))
    for i in arange(len(gamma)):
        Res_load_sum[i] = res_load_sum(alpha_w[i],gamma[i],1e3,1e3) 
        print i   
        print Res_load_sum[i]
        print res_load_sum(alpha_w[i], gamma[i], 1e-8, 1e3)
        print res_load_sum(alpha_w[i], gamma[i], 1e-8, 1e3) - Res_load_sum[i] > acc*Res_load_sum[i]
        print ' ' 


        if res_load_sum(alpha_w[i], gamma[i], 1e-8, 1e3) - Res_load_sum[i] > acc*Res_load_sum[i]:

            P_in_0 = brentq(lambda P_in: res_load_sum(alpha_w[i], gamma[i], P_in, 1e3)-(1+acc/2.)*Res_load_sum[i], 1e-8, 10)
            P_out_0 = brentq(lambda P_out: res_load_sum(alpha_w[i], gamma[i], P_in_0, P_out)-(1+acc)*Res_load_sum[i], 1e-8, 10)

            P_out_1 = brentq(lambda P_out: res_load_sum(alpha_w[i], gamma[i], 1e3, P_out)-(1+acc/2.)*Res_load_sum[i], 1e-8, 10)
            P_in_1 = brentq(lambda P_in: res_load_sum(alpha_w[i], gamma[i], P_in, P_out_1)-(1+acc)*Res_load_sum[i], 1e-8, 10)

            P_in_[i] = [P_in_0, P_in_1]
            P_out_[i] = [P_out_0, P_out_1]
            
        else:
            P_in_[i] = [0,0]
            P_out_[i] = [0,0]

    return Res_load_sum, P_in_, P_out_



################################################################################################
################################################################################################



def get_policy_2_storage_modified(mismatch, eta_in=1., eta_out=1., CS=NaN, P_in=None, P_out=None, fancy_storage=False):
    """ This function behaves like Morten's get_policy_2_storage(). However, it allows limiting the charging and discharging capacities of the storage. """
    
    if P_in==None:
        ## No constraints on charging
        P_in = amax(mismatch)
    if P_out==None:
        ## No constraints on discharging
        P_out = -amin(mismatch)
        
    ## The mismatch is cut by the maximum charging and discharging capacities.
    mismatch_ = lambda P_in, P_out: amax([amin([mismatch,P_in*ones_like(mismatch)],axis=0),-P_out*ones_like(mismatch)],axis=0) 
    
    ## The cut mismatch is feed into Mortens code to produce the cut reduced mismatch.     
    if fancy_storage:
        mismatch_r_, CS_used = get_storage(mismatch_(P_in,P_out),eta_in,eta_out,CS)   
    else:
        mismatch_r_, CS_used = get_policy_2_storage(mismatch_(P_in,P_out),eta_in,eta_out,CS)   
        
    ## Calculate the real reduced mismatch.
    mismatch_r = mismatch_r_ + (mismatch - mismatch_(P_in,P_out))
        
    ## REturns the reduced mismatch and the storage energy capacity that whas actually used (0<=CS_used<=CS).
    return mismatch_r, CS_used


################################################################################################
################################################################################################




#
# get_min_storage_cap_alt(L, Gw, Gs, gamma=1, alpha_w=1., CS=6)
# Storage_benefit, P_in_, P_out_ = get_min_storage_cap_alt(L, GW, GS, gamma=1, alpha_w=1., CS=None, acc=1e-2,returnall=False,eta_charge=1.,eta_discharge=1.):
def get_min_storage_cap_alt(L, GW, GS, gamma=1, alpha_w=1., CS=None, acc=1e-2,returnall=False,eta_charge=1.,eta_discharge=1.):
    """Finds upper and lower limits of the storage in and output capacities. acc sets the relative increase in balancing energy that is acceptable.
    WARNING! This function is slow to evaluate."""

    L, GW, GS, gamma, alpha_w = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2), array(gamma,ndmin=1), array(alpha_w,ndmin=1)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    if len(alpha_w)==1:
        alpha_w = alpha_w*ones_like(gamma)
    
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w, gamma: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
    mismatch_ = lambda alpha_w, gamma, P_in, P_out: amax([amin([mismatch(alpha_w, gamma),P_in*ones_like(l)],axis=0),-P_out*ones_like(l)],axis=0)

    ## Summed residual load after storage.
    res_load_sum = lambda alpha_w, gamma, P_in, P_out: sum(get_positive(-get_policy_2_storage(mismatch_(alpha_w,gamma,P_in,P_out),storage_capacity = CS)[0]) + get_positive(-mismatch(alpha_w, gamma)-P_out))
    
    ## Difference between residual load before and after storage.
    storage_benefit = lambda alpha_w, gamma, P_in, P_out: sum(get_positive(-mismatch(alpha_w, gamma))) - res_load_sum(alpha_w,gamma,P_in,P_out)
    
    Storage_benefit, E_surplus, E_residual, E_discharge, E_charge, N_cycles = zeros_like(gamma), zeros_like(gamma), zeros_like(gamma), zeros_like(gamma), zeros_like(gamma), zeros_like(gamma)
    P_in_ = zeros((len(gamma),2))
    P_out_ = zeros((len(gamma),2))
    for i in arange(len(gamma)):
        Storage_benefit[i] = storage_benefit(alpha_w[i],gamma[i],1e3,1e3) 

        ## Values to be returned
        E_surplus[i] = mean(get_positive(mismatch(alpha_w[i],gamma[i]))) ## Mean hourly surplus before storage.
        E_residual[i] = mean(get_positive(-mismatch(alpha_w[i],gamma[i]))) ## Mean hourly residual load before storage. 
        E_discharge[i] = (1.-acc)*Storage_benefit[i]/len(l) ## Mean hourly output of the storage. (almost same as 'Storage_benefit')
        E_charge[i] = E_discharge[i]/eta_discharge/eta_charge ## Mean hourly input of the storage.
        N_cycles[i] = (E_discharge[i]/eta_discharge)/CS ## Average number of accumulated storage cycles per hour.

        print i,
        sys.stdout.flush()    

        if acc==0:  #Intended for quick calculation where the charging and discharging powers are not relevant.
            P_in_[i] = [NaN,NaN]
            P_out_[i] = [NaN,NaN]
        #elif Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 1e-8, 1e3) > acc*Storage_benefit[i]:
        elif ((1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 1e-8, 10) > 0)*((1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 10, 1e-8) > 0):

            P_in_0 = brentq(lambda P_in: (1-acc/2.)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], P_in, 1e3), 1e-8, 10)
            P_out_0 = brentq(lambda P_out: (1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], P_in_0, P_out), 1e-8, 10)

            P_out_1 = brentq(lambda P_out: (1-acc/2.)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 1e3, P_out), 1e-8, 10)
            P_in_1 = brentq(lambda P_in: (1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], P_in, P_out_1), 1e-8, 10)

            P_in_[i] = [P_in_0, P_in_1]
            P_out_[i] = [P_out_0, P_out_1]
            
        else:
            P_in_[i] = [0,0]
            P_out_[i] = [0,0]

    if returnall==True:
        return Storage_benefit, P_in_, P_out_, E_surplus, E_residual, E_discharge, E_charge, N_cycles
    else:
        return Storage_benefit, P_in_, P_out_


################################################################################################
################################################################################################




#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=6,acc=1e-8,txtlabel='optimal mix',savelabel='acc1e-8_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=6,acc=1e-6,txtlabel='optimal mix',savelabel='acc1e-6_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=6,acc=1e-3,txtlabel='optimal mix',savelabel='acc1e-3_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=6,acc=1e-2,txtlabel='optimal mix',savelabel='acc1e-2_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=6,acc=1e-1,txtlabel='optimal mix',savelabel='acc1e-1_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=6,acc=5e-1,txtlabel='optimal mix',savelabel='acc5e-1_OptMix')
#
#
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=.254,acc=1e-8,txtlabel='optimal mix',savelabel='acc1e-8_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=.254,acc=1e-3,txtlabel='optimal mix',savelabel='acc1e-3_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=.254,acc=1e-2,txtlabel='optimal mix',savelabel='acc1e-2_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=.254,acc=1e-1,txtlabel='optimal mix',savelabel='acc1e-1_OptMix')
#  plot_min_storage_cap(gamma=linspace(0,1.05,21),alpha_w=None,CS=.254,acc=5e-1,txtlabel='optimal mix',savelabel='acc5e-1_OptMix')




def plot_min_storage_cap(ISO='DK', gamma=linspace(0,1.05,11), alpha_w=1., CS=6, acc=1e-6, txtlabel='', savelabel=''):

    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    ## Use optimal path? If yes, the optimal path WITHOUT storage is used.
    if alpha_w==None:
        alpha_w, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=None, returnall=True)
    
    # Function
    Res_load_sum, P_in_, P_out_ = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w, CS, acc)
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10
    
    close(1); figure(1); clf()
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])
    
    plot(gamma*mean(L)*365*24/1e3,mean(L)*P_in_.transpose()[0],'k-',label='Filling (upper)')
    plot(gamma*mean(L)*365*24/1e3,mean(L)*P_in_.transpose()[1],'r-',label='Filling (lower)')
    plot(gamma*mean(L)*365*24/1e3,mean(L)*P_out_.transpose()[0],'r--',label='Extraction (upper)')
    plot(gamma*mean(L)*365*24/1e3,mean(L)*P_out_.transpose()[1],'k--',label='Extraction (lower)') 
    
    axvline(.25*mean(L)*365*24/1e3,ls='--',color='k',lw=1)
    axvline(.5*mean(L)*365*24/1e3,ls='--',color='k',lw=1)
    axvline(.75*mean(L)*365*24/1e3,ls='--',color='k',lw=1)
    axvline(1.*mean(L)*365*24/1e3,ls='--',color='k',lw=1)

    text(.25*mean(L)*365*24/1e3-.3,.05*amax(mean(L)),r'25%',weight='semibold',fontsize=10,ha='right')
    text(.5*mean(L)*365*24/1e3-.3,.05*amax(mean(L)),r'50%',weight='semibold',fontsize=10,ha='right')
    text(.75*mean(L)*365*24/1e3-.3,.05*amax(mean(L)),r'75%',weight='semibold',fontsize=10,ha='right')
    text(1.0*mean(L)*365*24/1e3-.3,.05*amax(mean(L)),r'100%',weight='semibold',fontsize=10,ha='right')
    
    xlabel('Wind plus solar energy ('+txtlabel+') [TWh]')
    ylabel('Power [GW]')

    axis(xmin=0,xmax=amax(gamma*mean(L)*365*24/1e3),ymin=0,ymax=mean(L)*1.5)

    leg = legend(loc='upper left',ncol=2,title='Storage power capacities (vol. {0:.0f} GWh)'.format(CS*mean(L)));
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    tight_layout()
    save_file_name = 'plot_min_storage_cap_'+savelabel+'_'+ '10xCS{0:0.0f}h'.format(10*CS) +'_'+ISO+'.pdf'
    save_figure(save_file_name)


################################################################################################
################################################################################################


    
#   2013:
#
#   plot_min_storage_cap_summary(gamma=np.linspace(0,1.05,31))
#
#   plot_min_storage_cap_summary(gamma=np.linspace(0,1.05,31),CS=.254)
#
def plot_min_storage_cap_summary(ISO='DK', gamma=np.linspace(0,1.05,5), gain=[.99,.9], linestyle=['-','--'], CS=2.54, txtlabel='', savelabel=''):
    """Compares wind only and optimal mix scenarios"""

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    annual_TWh = mean(L)*365*24/1e3
    
    ## Use optimal path, the optimal path WITHOUT storage is used.
    alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=None, returnall=True)
    alpha_w_wind = 1.
    
    
    Res_load_sum_opt, Res_load_sum_wind = zeros((len(gain),len(gamma))), zeros((len(gain),len(gamma)))
    P_in_opt_, P_out_opt_ = zeros((len(gain),len(gamma))), zeros((len(gain),len(gamma)))
    P_in_wind_, P_out_wind_ = zeros((len(gain),len(gamma))), zeros((len(gain),len(gamma)))
    
    for i in arange(len(gain)):
        
        acc = 1. - gain[i]
        
        #Optimal mix
        print L, Gw, Gs, gamma, alpha_w_opt, CS, acc
        Res_load_sum_opt[i], P_in_opt, P_out_opt = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w_opt, CS, acc)
        P_in_opt_[i], P_out_opt_[i] = mean(P_in_opt.transpose(),axis=0), mean(P_out_opt.transpose(),axis=0)
    
        #Wind only
        Res_load_sum_wind[i], P_in_wind, P_out_wind = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w_wind, CS, acc)
        P_in_wind_[i], P_out_wind_[i] = mean(P_in_wind.transpose(),axis=0), mean(P_out_wind.transpose(),axis=0)
      
    max_y = 1e3*1.2*mean(L)*amax(concatenate([P_in_opt_,P_out_opt_,P_in_wind_,P_out_wind_]))
      
    ### Chargeing storage                      
    close(1); figure(1); clf()

    #plot(gamma,balancing_power_mean,'-',color=(.5,.5,.5),lw=2.5,label='Mean')
    for i in arange(len(gain)):
        plot(gamma*mean(L)*365*24/1e3,1e3*mean(L)*P_in_opt_[i],ls=linestyle[i],lw=1,color='k',label=r'{0:.0f}% (bal. opt. mix)'.format(gain[i]*100))

    #plot(gamma*mean(L)*365*24/1e3,Res_load_sum_opt[0],'-',color=color_wind,lw=3, label='Wind')
    for i in arange(len(gain)):
        plot(gamma*mean(L)*365*24/1e3,1e3*mean(L)*P_in_wind_[i],ls=linestyle[i],lw=1.5,color=color_wind,label=r'{0:.0f}% (wind only)'.format(gain[i]*100))

    dx = 0.02*annual_TWh
    plot_vertical_line_and_label(0.25*annual_TWh,max_y/6.,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,max_y/6.,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,max_y/6.,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,max_y/6.,r'100%',dx)

    axis(xmin=0,xmax=amax(gamma*mean(L)*365*24/1e3),ymin=0,ymax=max_y)
    
    ylabel('Charging power [MW]')
    xlabel('Wind plus solar energy [TWh/yr]')
    
    legend(loc='upper left',ncol=1,title='gain');

    #add_duplicate_yaxis(gcf(),unit_multiplier=1./mean(L),label='[av.l.h.]')
    
    tight_layout()
    save_file_name = 'plot_min_storage_cap_summary_charging_'+savelabel+'_'+ '10xCS{0:0.0f}h'.format(10*CS) +'_'+ISO+'.pdf'
    save_figure(save_file_name)
    
    ### Dischargeing storage                      
    close(1); figure(1); clf()
    
    #plot(gamma,excess_power_mean,'-',color=(.5,.5,.5),lw=2.5,label='Mean')
    for i in arange(len(gain)):
        plot(gamma*mean(L)*365*24/1e3,1e3*mean(L)*P_out_opt_[i],ls=linestyle[i],lw=1,color='k',label=r'{0:.0f}% (bal. opt. mix)'.format(gain[i]*100))    
    
    #plot(gamma,excess_power_mean_wind,'-',color=color_wind,lw=3, label='Wind')
    for i in arange(len(gain)):
        plot(gamma*mean(L)*365*24/1e3,1e3*mean(L)*P_out_wind_[i],ls=linestyle[i],lw=1.5,color=color_wind,label=r'{0:.0f}% (wind only)'.format(gain[i]*100))
    
    dx = 0.02*annual_TWh
    plot_vertical_line_and_label(0.25*annual_TWh,max_y/6.,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,max_y/6.,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,max_y/6.,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,max_y/6.,r'100%',dx)       
                
    axis(xmin=0,xmax=amax(gamma*mean(L)*365*24/1e3),ymin=0,ymax=max_y)

    ylabel('Discharging power [MW]')
    xlabel('Wind plus solar energy [TWh/yr]')

    legend(loc='upper left',ncol=1,title='gain');

    #add_duplicate_yaxis(gcf(),unit_multiplier=1./mean(L),label='[av.l.h.]')

    tight_layout()
    save_file_name = 'plot_min_storage_cap_summary_discharging_'+savelabel+'_'+ '10xCS{0:0.0f}h'.format(10*CS) +'_'+ISO+'.pdf'
    save_figure(save_file_name)    


################################################################################################
################################################################################################



#
#   plot_min_storage_cap_fixed_gamma_summary( gamma=.75)
#   plot_min_storage_cap_fixed_gamma_summary( gamma=.75,gain=array([0,.9,.99,1-1e-8]))

def plot_min_storage_cap_fixed_gamma_summary( gamma=.5, gain=linspace(0,1-1e-8,111), CS=2.54, linestyle=['-','--','-.',':'], txtlabel='', savelabel='',ISO='DK'):
    
    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    ## Use optimal path, the optimal path WITHOUT storage is used.
    alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_0, mismatch_opt, res_load_sum_1p = get_optimal_path_balancing(L, Gw, Gs, gamma, CS=None, returnall=True)
    alpha_w_wind = 1.
    
    Res_load_sum_opt, Res_load_sum_wind = np.zeros(gain.shape), np.zeros(gain.shape)
    P_in_opt, P_out_opt = np.zeros((len(gain),2)), np.zeros((len(gain),2))
    P_in_wind, P_out_wind = np.zeros((len(gain),2)), np.zeros((len(gain),2))
    
    for i in arange(len(gain)):
        
        acc = 1. - gain[i]
        
        Res_load_sum_opt[i], P_in_opt[i], P_out_opt[i] = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w_opt, CS, acc)
        Res_load_sum_wind[i], P_in_wind[i], P_out_wind[i] = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w_wind, CS, acc)
    
    max_y = 1e3*1.2*mean(L)*amax(concatenate([P_in_opt,P_out_opt,P_in_wind,P_out_wind]))
    
    plt.close(1); plt.figure(1); plt.clf()
    
    plt.plot(100*gain,1e3*mean(L)*mean(P_in_opt.transpose(),axis=0),'k-',label='Charging (bal. opt. mix)')
    plt.plot(100*gain,1e3*mean(L)*mean(P_out_opt.transpose(),axis=0),'k--',label='Discharging (bal. opt. mix)')
    
    plt.plot(100*gain,1e3*mean(L)*mean(P_in_wind.transpose(),axis=0),'-',lw=1.5,color=color_wind,label='Charging (wind only)')
    plt.plot(100*gain,1e3*mean(L)*mean(P_out_wind.transpose(),axis=0),'--',lw=1.5,color=color_wind,label='Discharging (wind only)')
    
    plt.xlabel(r'gain [%]')
    plt.ylabel('Power [MW]')
    
    plt.axvline(90,ls='--',color='k',lw=1,alpha=.75)
    plt.text(90-2,max_y*2/3.,r'90%',weight='semibold',ha='right')
    
    plt.axis(xmin=0,xmax=100,ymin=0)
    
    plt.legend(loc='upper left',ncol=1);
        
    plt.tight_layout()
    
    save_file_name = 'plot_min_storage_cap_fixed_gamma_summary_'+'100x_gamma_{0:.0f}'.format(100*gamma)+'_'+'10xCS{0:0.0f}h'.format(10*CS)+'_'+ISO+'.pdf'

    save_figure(save_file_name)  
    
    print 'Gain: ' + str(gain)
    print 'P_in_opt: ' + str(mean(L)*mean(P_in_opt.transpose(),axis=0))
    print 'P_out_opt: ' + str(mean(L)*mean(P_out_opt.transpose(),axis=0))
    print 'P_in_wind: ' + str(mean(L)*mean(P_in_wind.transpose(),axis=0))
    print 'P_out_wind: ' + str(mean(L)*mean(P_out_wind.transpose(),axis=0))



################################################################################################
################################################################################################



#  get_storage_summary_table(ISO='DK',gamma=[.5,.75,1.],CS=array([0.0,3])/3.94329, alpha_w=[1.], storage_gain=.90)
#
def get_storage_summary_table(ISO='DK',gamma=[.5,.75,1.],CS=array([0.1,1.,10,30,nan])/3.94329, alpha_w=None, storage_gain=.90):
    """ All numbers are per year. CS=NaN gives results for seasonal storage at the relevant mix. alpha_w=None is balancing optimal mix, alpha_w=NaN is seasonal optimal mix."""
    
    ## Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)    
    
    gamma_,CS_ = meshgrid(gamma,CS)
    
    for i in arange(len(gamma_.flat)):
    
        E_surplus, E_residual, N_cycles, E_charge, E_discharge, P_charge, P_discharge, alpha_w_used, CS_used = get_storage_summary(ISO,gamma_.flat[i],alpha_w,CS_.flat[i],storage_gain)
    
        print ' '
    
        print 'Gamma: {0:0.2f}, alpha_W: {1:0.2f}, CS: {2:0.2f} GWh'.format(gamma_.flat[i],alpha_w_used[0],CS_used*mean(L))
        print '-------------------------------------'
        print 'E_surplus: {0:0.3f} TWh/yr, E_residual: {1:0.3f} TWh/yr'.format(E_surplus[0]*mean(L)*365*24/1e3, E_residual[0]*mean(L)*365*24/1e3)
        print 'E_charge: {0:0.3f} TWh/yr, E_discharge: {1:0.3f} TWh/yr, N_cycles: {2:0.1f} cycles/yr'.format(E_charge[0]*mean(L)*365*24/1e3, E_discharge[0]*mean(L)*365*24/1e3, N_cycles[0]*365*24)
        print 'P_charge: {0:0.0f} MW, P_discharge: {1:0.0f} MW'.format(mean(P_charge[0])*mean(L)*1e3,mean(P_discharge[0])*mean(L)*1e3)
        
        print ' '
        
        #print E_surplus, E_residual, N_cycles, E_charge, E_discharge, P_charge, P_discharge, alpha_w, CS
        #print ' '


################################################################################################
################################################################################################

        
def get_storage_summary(ISO='DK', gamma=1., alpha_w=None, CS=nan, storage_gain=.99,country_data=None):   
    """CS=NaN gives results for seasonal storage at the relevant mix. alpha_w=None is balancing optimal mix, alpha_w=NaN is seasonal optimal mix. country_data=(t, L, Gw, Gs, datetime_offset, datalabel)"""
    
    ## Load data
    if country_data==None:
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    else:
        t, L, Gw, Gs, datetime_offset, datalabel = country_data
    
    ## Find mix if it is not specified directly as a fraction.
    if alpha_w==None:
        ## Identify balancing optimal mix
        if isnan(CS):
            alpha_w = get_optimal_mix(L, Gw, Gs, gamma, CS=None)
            CS = get_policy_2_storage(gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - L/mean(L))[1]
        else:
            alpha_w = get_optimal_mix(L, Gw, Gs, gamma, CS=None)
    elif isnan(alpha_w):
        ## Identify seasonal optimal mix and find the seasonal optimal storage energy capacity
        alpha_w, alpha_w_1p_interval, res_load_sum, mismatch, res_load_sum_1p, CS_pp = get_optimal_mix(L, Gw, Gs, gamma, CS=NaN, returnall=True)
    
    if isnan(CS):
        CS = get_policy_2_storage(gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - L/mean(L))[1]
    
    ## Calculate optimal storage dynamics
    Storage_benefit, P_charge, P_discharge, E_surplus, E_residual, E_discharge, E_charge, N_cycles = get_min_storage_cap_alt(L, Gw, Gs, gamma, alpha_w, CS, 1.-storage_gain,returnall=True)
    
    ## Average hourly values are returned.
    return E_surplus, E_residual, N_cycles, E_charge, E_discharge, P_charge, P_discharge, alpha_w, CS


################################################################################################
################################################################################################


##
#
# plot_seasonal_storage_singularity(gamma=linspace(0,1.22,111),alpha_w=[1,None])
#
#        
def plot_seasonal_storage_singularity(ISO='DK', gamma=linspace(0,1.22,5), alpha_w=[1,None,nan],textlabel=['Wind-only','Opt. mix','Seasonal opt. mix'],line_color=[color_wind,'k','k'],ls=['-','-','--']):
    """alpha_w: None=balancing optimal mix, NaN=Seasonal optimal mix"""
    
    ## Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    annual_TWh = mean(L)*365*24/1e3

    CS = zeros((len(alpha_w),len(gamma)))
    
    for i in arange(len(alpha_w)):
        for j in arange(len(gamma)):
            E_surplus, E_residual, N_cycles, E_charge, E_discharge, P_charge, P_discharge, alpha_w_, CS[i][j] = get_storage_summary(ISO, gamma[j], alpha_w[i], CS=nan, storage_gain=1.00, country_data=(t, L, Gw, Gs, datetime_offset, datalabel))    
      
    ### Plot storage singularity                      
    close(1); figure(1); clf()
      
    for i in arange(len(alpha_w)):
        plot(gamma*annual_TWh,CS[i]*mean(L),ls=ls[i],color=line_color[i],label=textlabel[i],lw=1.5)  
    
    dx = 0.02*annual_TWh
    plot_vertical_line_and_label(0.25*annual_TWh,11500/6.,r'25%',dx)
    plot_vertical_line_and_label(0.5*annual_TWh,11500/6.,r'50%',dx)
    plot_vertical_line_and_label(.75*annual_TWh,11500/6.,r'75%',dx)
    plot_vertical_line_and_label(1.*annual_TWh,11500/6.,r'100%',dx)
    
    xlabel('Wind plus solar energy [TWh/yr]')
    ylabel('Seasonal volume [GWh]')
      
    legend(loc='upper left')  
      
    axis(xmin=0,xmax=amax(gamma*annual_TWh),ymin=0,ymax=11500)  
      
    tight_layout()
    save_file_name = 'plot_seasonal_storage_singularity'+'_'+ISO+'.pdf'
    save_figure(save_file_name)  


################################################################################################
################################################################################################


##
# plot_storage_balancing_synergy(CS=linspace(0,13,51))
#
def plot_storage_balancing_synergy(ISO='DK',CS=linspace(0,13,5),gamma=[.5,.75,1.], alpha_w=[1,None],textlabel=['Wind only','Bal. opt. mix','Seasonal opt. mix']):

    ## Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    annual_GWh = mean(L)*365*24

    E_residual, E_discharge = zeros((len(gamma),len(alpha_w),len(CS))), zeros((len(gamma),len(alpha_w),len(CS)))
    E_seasonal_storage = zeros((len(gamma),len(alpha_w)))
    
    for k in arange(len(gamma)):
        for i in arange(len(alpha_w)):
            E_seasonal_storage[k][i] = get_storage_summary(ISO, gamma[k], alpha_w[i], CS=NaN, storage_gain=1.00, country_data=(t, L, Gw, Gs, datetime_offset, datalabel))[4]
            
            for j in arange(len(CS)):
                E_surplus, E_residual[k][i][j], N_cycles, E_charge, E_discharge[k][i][j], P_charge, P_discharge, alpha_w_, CS_ = get_storage_summary(ISO, gamma[k], alpha_w[i], CS=CS[j], storage_gain=1.00, country_data=(t, L, Gw, Gs, datetime_offset, datalabel))   
      
    ### Plot storage singularity                      
    close(1); figure(1); clf()
    
    for k in arange(len(gamma)):
        for i in arange(len(alpha_w)):
            
            color, ls, mixtxt = get_color_and_linestyle(alpha_w[i],gamma=gamma[k])
            
            textlabel = r'{0:0.0f}%'.format(100*gamma[k]) + mixtxt
            
            plot(CS*mean(L),100*E_discharge[k][i]/E_seasonal_storage[k][i],ls=ls,color=color,label=textlabel,lw=1.5)
   
    #axhline(100,color='k',ls='-',alpha=.75)
   
    xlabel('Storage volume [GWh]')
    ylabel(r'Stored share of surplus [%]')
      
    axis(xmin=0,xmax=amax(CS*mean(L)),ymin=0,ymax=101)
    
    leg = legend(loc='upper left',title='VRE gross share',ncol=2)  
    frame = leg.get_frame()
    frame.set_alpha(0.75)
          
                      
    tight_layout()
    save_file_name = 'plot_storage_balancing_synergy'+'_'+ISO+'.pdf'
    save_figure(save_file_name)
    
    

def plot_vertical_line_and_label(x,y,textlabel=None,dx=None,ls='--',color='k',lw=1,ymax=1,ymin=0,alpha=.75):

    if dx==None:
        dx=0.05*x

    axvline(x,ls=ls,color=color,lw=lw,ymin=ymin,ymax=ymax,alpha=alpha)
    text(x-dx,y,textlabel,weight='semibold',ha='right')
    

################################################################################################
################################################################################################


    
      
######
# Convinient access to ISET country data.
# Main function:  get_ISET_country_data()
###

def get_ISET_country_data(ISO='DK',path='./data/',silent=True):
    """
    Returns data for a specific country in the ISET data set. 
    
    The country is specified using ISO two letter names. For a list of names use get_ISET_country_names().
    
    The function retrives/stores data from/in ./data as default. If the data is not available the function
     attempts to download it from pepsi. This requires ssh access: 
     ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk 
     
    Filenames are: ISET_country_ISO.npz
    
    Returns
    -------
    t: array
        hours numbered from 0 to N-1
    L: array
        load in MW
    Gw: array
        normalized wind power generation
    Gs: array
        normalized solar power generation
    datetime_offset: scalar
        num2date(datetime_offset) yields the date and hour of t=0
    datalabel: string
        country ISO
        
    Example: Returns the Danish load, wind and solar power generation
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('DK')
    
    """
    
    if not valid_ISO(ISO):
        sys.exit("Error (43nlksd): No such country ISO ({0}). For a list of names use get_ISET_country_names().".format(ISO))
    
    try:
        a = get_ISET_country_data.cache
    except AttributeError:
        get_ISET_country_data.cache = dict()
    
    try:
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data.cache[ISO];
    except KeyError:
        filename = 'ISET_country_' + ISO + '.npz'
        
        try:
            #Load the data file if it exists:
            npzfile = np.load(path + filename)
            if ~silent:
                print 'Loaded file: ', path + filename
                sys.stdout.flush()
            
        except IOError:
            print 'Datafile does not exist:', path + filename
            print 'Trying to download data from pepsi...'
            sys.stdout.flush()
            try: 
                #Normalized data: t, l, Gw, Gs, datetime_offset, datalabels = get_data_countries(schema='norm_agg_avg_1hour_pdata_caps_eu2020',localhost=True);
                t, L, GW, GS, datetime_offset, datalabels = get_data_countries(localhost=True);
            except:
                sys.exit("Error (sdf3dz1): Could not connect to pepsi. Setup ssh access first: ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk")
                
            #Save all country files:
            for i in arange(len(datalabels)):
                ISO_ = ISET2ISO_country_codes(datalabels[i])
                filename_ =  'ISET_country_' + ISO_ + '.npz'
                np.savez(path + filename_,t=t, L=L[i], Gw=GW[i]/mean(GW[i]), Gs=GS[i]/mean(GS[i]), datetime_offset=datetime_offset, datalabel=ISO_)
                if ~silent:
                    print 'Saved file: ', path + filename_
                    sys.stdout.flush()
             
            #Load the relevant file now that it has been created:       
            npzfile = np.load(path + filename)
            if ~silent:
                print 'Loaded file: ', path + filename
                sys.stdout.flush()
        
        t, L, Gw, Gs, datetime_offset, datalabel = npzfile['t'], npzfile['L'], npzfile['Gw'], npzfile['Gs'], npzfile['datetime_offset'], npzfile['datalabel']
        npzfile.close()
    
        get_ISET_country_data.cache[ISO] = (t, L, Gw, Gs, datetime_offset, datalabel)
    
    return t, L, Gw, Gs, datetime_offset, datalabel

def valid_ISO(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    if ISO=='DK-hist':
        return 1
    else:
        return (ISO in table['ISO'])
    
def ISO2ISET_country_codes(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISET = table['ISET'][find(table['ISO']==ISO)][0]
    
    return ISET

def ISET2ISO_country_codes(ISET='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISO = table['ISO'][find(table['ISET']==ISET)][0]
    
    return ISO

def get_ISET_country_names(filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)

    print 'Index\tISO\tISET\tName'
    print '====================================='
    for i in arange(len(table)):
        print '{0}\t{1}\t{2}\t{3}'.format(i, table[i]['ISO'], table[i]['ISET'], table[i]['name'])
    sys.stdout.flush()