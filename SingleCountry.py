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
from MortenStorage import get_policy_2_storage

#Specific functions, not sure if all are actually used
from scipy.optimize import brentq
from scipy.optimize import fmin
from scipy.optimize import leastsq
from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles
from mpl_toolkits.axes_grid1 import make_axes_locatable
import calendar #Used to place day labels on plots

##Colors, should be moved to a color module
color_wind = (0.5,0.7,1.)
color_solar = (1.,.8,0.)
bg_color = (.75,.0,.0)
color_edge = (.4,.4,.4)

def plot_storage_balancing_synergy(ISO='DK', gamma=[.5,.75,1.,1.25,1.5], CS=linspace(0,24,5)):

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
    
    conv_x = lambda gamma, alpha_w: mean(L)*gamma*alpha_w
    conv_y = lambda gamma, alpha_w: mean(L)*gamma*(1-alpha_w)

    plot(conv_x(gamma[mask],alpha_w_opt[mask]),conv_y(gamma[mask],alpha_w_opt[mask]),'w-',lw=2)

    #fill_between(conv_x(gamma,amax(gamma)*ones(len(gamma))),amax(gamma)*ones(len(gamma)),color=bg_color,edgecolor=(0,0,0,0))
    fill(mean(L)*concatenate([gamma,[0]]),mean(L)*concatenate([zeros(gamma.shape),[amax(gamma)]]),color=bg_color,edgecolor=(0,0,0,0))

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
    plot(.2*array([0.,mean(L)]),.2*array([mean(L),0]),'k--',lw=1)
    text(.05*mean(L),.2*mean(L)-.025*mean(L),r'20%',weight='semibold',fontsize=10)
    
    plot(.5*array([0.,mean(L)]),.5*array([mean(L),0]),'k--',lw=1)
    text(.05*mean(L),.5*mean(L)-.025*mean(L),r'50%',weight='semibold',fontsize=10)
    
    plot(1.*array([0.,mean(L)]),1.*array([mean(L),0]),'k--',lw=1)
    text(.05*mean(L),1.*mean(L)-.025*mean(L),r'100%',weight='semibold',fontsize=10)
    

    ## Present state
    gamma_present,alpha_w_present, label_present = .22*mean(L), 0.0, '2010'
    plot(gamma_present,alpha_w_present,'o',ms=10,mfc='k',mec='w',mew=1)
    text(gamma_present+.025*mean(L),alpha_w_present+.025*mean(L),label_present,va='bottom',ha='left',weight='semibold',fontsize=10)

    axis(xmin=0,xmax=mean(L)*amax(gamma),ymin=0,ymax=mean(L)*amax(gamma))
    axis('scaled')
    xlabel(r'Av. wind power [GWh/h]')
    ylabel(r'Av. solar PV power [GWh/h]')
    
    ax=gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    txtlabels = ['0-1 pp','1-5 pp',' >5 pp']

    leg = legend(pp,txtlabels,loc='upper right',ncol=2,title='Deviation from optimal mix');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

   # add_duplicate_yaxis(gcf(),unit_multiplier=1.,label=r'Solar fraction $\alpha^{'+ISO+'}_s$',invert_axis=True)

    tight_layout()
    save_file_name = 'plot_country_optimal_mix_vs_gamma_triangle_'+ISO+'_CS_'+str(CS)+'.pdf'
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


def plot_surplus_bar(ISO='DK',gamma_bar=array([.25,.50,.75,1.]),alpha_w=1,CS=None):

    surplus_bar = zeros(gamma_bar.shape)
    for i in arange(len(gamma_bar)):
        L, wind_local, solar_local, curtailment, filling, extraction = get_compare_VRES_load(ISO, gamma_bar[i], alpha_w, CS)

        surplus_bar[i] = mean(curtailment)

    gamma = linspace(amin(gamma_bar)-.1,amax(gamma_bar)+.1,111)
    surplus = zeros(gamma.shape)
    for i in arange(len(gamma)):
        L, wind_local, solar_local, curtailment, filling, extraction = get_compare_VRES_load(ISO, gamma[i], alpha_w, CS)

        surplus[i] = mean(curtailment)
    

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
        if height<=0.4:
            height = height+0.1
        else:
            height = height/2
            
        gca().text(bars[i].get_x()+bars[i].get_width()/2., height, '{0:0.1f}%\nof\n{1:0.1f} GW'.format(surplus_bar[i]*100/(gamma_bar[i]*mean(L))+1e-10,gamma_bar[i]*mean(L)),
                ha='center', va='center')
    
    
    #autolabel_bars(bars,gca(),heighttxt='{0:0.1f}%\nof\n'+'{0:0.1f} GW'.format(),scale_factors=100/(gamma_bar*mean(L)))

    axis(xmin=amin(gamma_bar*100)-10,xmax=amax(gamma_bar*100)+10,ymin=0,ymax=1.05*amax(surplus_bar))
    xticks(gamma_bar*100,['{0:.0f}%'.format(x) for x in gamma_bar*100],va='top')
    xlabel(r'Target share of total electricity demand'.format(mean(L)))
    ylabel('Renewable surplus [GW]')

    ax=gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    tight_layout()
    save_file_name = 'plot_surplus_bar_'+ISO+'_CS_'+str(CS)+'.pdf'
    save_figure(save_file_name)

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

# plot_hourly_generation_alt(alpha_w=1.,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, wind only',label='week_wind')
# plot_hourly_generation_alt(alpha_w=None,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, optimal wind and solar mix',label='week_optimal')
#
#
# plot_hourly_generation_alt(alpha_w=1.,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, wind only',label='week_wind_storage',CS=7.4)
# plot_hourly_generation_alt(alpha_w=None,date_start=datestr2num('3-6-2000'),N_days=7,monday_offset=7,titletxt='Spring week, optimal wind and solar mix',label='week_optimal_storage',CS=7.4)

def plot_hourly_generation_alt(ISO='DK', gamma=0.5, alpha_w=.5, CS=None, date_start=datestr2num('1-1-2000'), N_days=30, monday_offset=5, titletxt='Denmark, Jan. 2000', label='TestFigure'):

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
    p = fill_between(t[mask],wind_local[mask]+solar_local[mask]+filling[mask],wind_local[mask]+solar_local[mask],color='g',edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+extraction[mask],wind_local[mask]+solar_local[mask],color='g',edgecolor=color_edge,lw=1)
    fill_between(t[mask],wind_local[mask]+solar_local[mask]+filling[mask]+curtailment[mask],wind_local[mask]+solar_local[mask]+filling[mask],color='r',edgecolor=color_edge,lw=1)
    
    p.set_facecolors("none")

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

    axis(xmin=t[mask[0]],xmax=t[mask[-1]],ymin=0,ymax=1.9*mean(L))

    ylabel('Power [GW]')

    day_names = array([calendar.day_abbr[mod(i,7)] for i in arange(N_days)+monday_offset])
    #day_names[find([d!='Mon' for d in day_names])] = ''
    xticks(t[mask[0]]+arange(N_days),day_names,rotation=-45,ha='left')

    if CS==None:
        pp = [pp_load[0],pp_wind,pp_solar,pp_curtailment]
        txtlabels = ['Load ({0})'.format(ISO),'Wind','Solar','Surplus']
        leg = legend(pp,txtlabels,loc='upper left',ncol=5,title=titletxt);
    else:
        pp = [pp_load[0],pp_filling,pp_wind,pp_curtailment,pp_solar,pp_storage]
        txtlabels = ['Load ({0})'.format(ISO),'Stored surplus','Wind','Remainder surplus','Solar','Storage']    
        leg = legend(pp,txtlabels,loc='upper left',ncol=4,title=titletxt);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    tight_layout()
    save_file_name = 'plot_hourly_generation_alt_'+ISO+'_'+label+'.pdf'
    save_figure(save_file_name)

# plot_monthly_summary('DK',gamma=0.5,alpha_w=None,label='optimal')
# plot_monthly_summary('DK',gamma=0.5,alpha_w=1,label='wind')
# plot_monthly_summary('DK',gamma=0.5,alpha_w=0,label='solar')
#
# plot_monthly_summary('DK',gamma=0.5,alpha_w=None,label='optimal_storage',CS=7.6)
# plot_monthly_summary('DK',gamma=0.5,alpha_w=1,label='wind_storage',CS=7.6)
# plot_monthly_summary('DK',gamma=0.5,alpha_w=0,label='solar_storage',CS=7.6)
#
# plot_monthly_summary('DK',gamma=1,alpha_w=None,label='100p_optimal')
# plot_monthly_summary('DK',gamma=1,alpha_w=1,label='100p_wind')
def plot_monthly_summary(ISO='DK', gamma=.5, alpha_w=None, CS=None, titletxt='Denmark, 2000-2007',label='TestFigure'):

    #Load data
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)

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
    save_file_name = 'plot_single_bar_summary_'+ISO+'_'+label+'.pdf'
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
    txtlabels = ['Wind','Solar','Storage','Surplus/Remainder']
    leg = legend(pp,txtlabels,loc='upper left',ncol=4,frameon=False);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    save_file_name = 'plot_single_bar_summary_legend.pdf'
    save_figure(save_file_name)
    

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


def get_balancing(L, GW, GS, gamma=1, alpha=1., CS=None,returnall=False):

	L, GW, GS, gamma, alpha = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2), array(gamma,ndmin=1), array(alpha,ndmin=1)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

	l = weighed_sum(L)
	Gw = weighed_sum(GW)	
	Gs = weighed_sum(GS)
	
	mismatch = lambda alpha_w, gamma: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
	
	if CS==None:
		res_load_sum = lambda alpha_w, gamma: sum(get_positive(-mismatch(alpha_w,gamma)))
	else:
		res_load_sum = lambda alpha_w, gamma: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w,gamma),storage_capacity = CS)[0]))
	
	Gamma, Alpha = meshgrid(gamma,alpha)
	
	Res_load_sum = zeros(Gamma.shape)
	for i in arange(size(Gamma)):
		Res_load_sum.flat[i] = res_load_sum(Alpha.flat[i],Gamma.flat[i])
	
	if returnall:
		return Res_load_sum, Gamma, Alpha
	else:
		return Res_load_sum

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

def get_optimal_path_balancing(L, GW, GS, gamma=linspace(0,1,5), p_interval=0.01, CS=None, returnall=False, normalized=True):
    """Wraper for get_optimal_mix_balancing(). This function allows gamma to be an array."""

    gamma = array(gamma,ndmin=1)
    
    if returnall==True:
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = get_optimal_mix_balancing(L, GW, GS, gamma[0], p_interval, CS, returnall, normalized)
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = expand2array(alpha_w_opt,gamma), expand2array(alpha_w_opt_1p_interval,gamma), expand2array(res_load_sum_opt,gamma), expand2array(mismatch_opt,gamma), expand2array(res_load_sum_1p,gamma)
        
        for i in arange(1,len(gamma)):
            alpha_w_opt[i], alpha_w_opt_1p_interval[i], res_load_sum_opt[i], mismatch_opt[i], res_load_sum_1p[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval, CS, returnall, normalized)
        
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p

    else:
        alpha_w_opt = zeros(gamma.shape)
        for i in arange(len(gamma)):
            alpha_w_opt[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval, CS, returnall, normalized)
        
        return alpha_w_opt
        

def get_optimal_mix_balancing(L, GW, GS, gamma=1., p_interval=0.01, CS=None, returnall=False, normalized=True, DefaultWind=True):

    L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l

    if CS==None:
        res_load_sum = lambda alpha_w: sum(get_positive(-mismatch(alpha_w))) - alpha_w*0.001*sign(DefaultWind-.5)
    else:
        res_load_sum = lambda alpha_w: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w),storage_capacity = CS)[0])) - alpha_w*0.001*sign(DefaultWind-.5)
    
    alpha_w_opt = fmin(res_load_sum,0.5,disp=False)

    if alpha_w_opt>1.:
        alpha_w_opt = 1.
    elif alpha_w_opt<0.:
        alpha_w_opt = 0.

    if normalized:
        if CS==None:
            mismatch_opt = mismatch(alpha_w_opt)
        else:
            mismatch_opt = get_policy_2_storage(mismatch(alpha_w_opt),storage_capacity = CS)[0]
    else:
        if CS==None:
            mismatch_opt = mismatch(alpha_w_opt)*mean(sum(L,axis=0))
        else:
            mismatch_opt = get_policy_2_storage(mismatch(alpha_w_opt),storage_capacity = CS)[0]*mean(sum(L,axis=0))
            
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
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p
    else:
        return alpha_w_opt







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
    
    return t, L, Gw, Gs, datetime_offset, datalabel

def valid_ISO(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
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