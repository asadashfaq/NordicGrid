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

#Custom functions
from Database_v1 import get_data_countries, get_data_regions
from MortenStorage import get_policy_2_storage

def generate_data_files():
    """Works from anywhere  if you setup an ssh tunnel first: ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk"""
    # t, L, GW, GS, datetime_offset, datalabels = get_data_countries(localhost=True);
    # t, l, Gw, Gs, datetime_offset, datalabels = get_data_countries(schema='norm_agg_avg_1hour_pdata_caps_eu2020',localhost=True);
    
    return 5
    
def ISET2ISO_country_codes(filename='ISET2ISO_country_codes.npy',path='./settings/'):

    return np.load(path+filename)
    
###
# DK: i=7; plot_country_optimal_mix_vs_gamma(L[i], GW[i], GS[i], gamma=linspace(0,2.05,31))
# DK, 6h storage:  i=7; plot_country_optimal_mix_vs_gamma(L[i], GW[i], GS[i], gamma=linspace(0,2.05,31),CS=6)
#
#(1.,.53,.20)
def plot_country_optimal_mix_vs_gamma(L, GW, GS, gamma=linspace(0,2,11),p_interval=[0.01,0.05,0.25], CS=None,linespec=['--','-.',':'],color=[(0.53,0.73,0.37),(1.,.82,.20),(.90,.27,.20)]):
	
	L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.

	bg_color = (.75,.0,.0)

	alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, res_load_sum_1p = zeros(len(gamma)), zeros((len(gamma),2)), zeros(len(gamma)), zeros(len(gamma))
	lower_bound, upper_bound, res_load_sum_p = zeros((len(p_interval),len(gamma))), zeros((len(p_interval),len(gamma))), zeros((len(p_interval),len(res_load_sum_1p)))
	for j in arange(len(p_interval)):
		for i in arange(len(gamma)):
			alpha_w_opt[i], alpha_w_opt_1p_interval[i], res_load_sum_opt[i], mismatch_opt, res_load_sum_1p[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval=p_interval[j], CS=CS, returnall=True, normalized=True)

		lower_bound[j] = alpha_w_opt_1p_interval.transpose()[0]
		upper_bound[j] = alpha_w_opt_1p_interval.transpose()[1]
		res_load_sum_p[j] = res_load_sum_1p/len(L[0])
	
	#mask = find((lower_bound[0]!=0) + (upper_bound[0]!=1))
	mask = arange(len(gamma) - amax([argmin(lower_bound[0][::-1]),argmax(upper_bound[0][::-1])]),len(gamma))
	
	res_load_sum_opt = res_load_sum_opt/len(L[0])
	
	#Set plot options	
	matplotlib.rcParams['font.size'] = 10

	figure(1); clf()

	gcf().set_dpi(300)
	gcf().set_size_inches([4.75,6])

	#Upper panel
	ax1 = axes([.11,.565,.885,.42])
	plot(gamma[mask],alpha_w_opt[mask],'w-',lw=2)
	
	fill_between(gamma,ones(len(gamma)),color=bg_color,edgecolor=(0,0,0,0))
	
	pp = list(zeros(1+len(p_interval)))
	for j in arange(len(p_interval))[::-1]:
		fill_between(gamma,lower_bound[j],upper_bound[j],color=color[j],lw=1,edgecolor=(0,0,0,0))
		pp[j] = Rectangle((0, 0), 1, 1, color=color[j],lw=1)
		
		plot(gamma,lower_bound[j],linespec[j],color='w',lw=2)
		plot(gamma,upper_bound[j],linespec[j],color='w',lw=2)

	axvline(.2,ls='--',color='k')
	#plot(gamma,(0.2 + (gamma-0.2)*0.75)/gamma,'k--')

	pp[-1] = Rectangle((0, 0), 1, 1, color=bg_color,lw=1)
	pp = tuple(pp)
	
	axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=1)
	xlabel(r'RES power generation coefficient $\gamma$')
	ylabel(r'Wind fraction $\alpha_W$')
	
	txtlabels = ('0-1 pp','1-5 pp','5-25 pp','>25 pp')
	
	leg = legend(pp,txtlabels,loc='lower right',ncol=2);
	ltext  = leg.get_texts();
	setp(ltext, fontsize='small')    # the legend text fontsize
	
	
	#Lower panel
	ax2 = axes([.11,.065,.885,.42])
	
	plot([0,1],[1,0],'-',color='k',alpha=.3)
	
	pp_solar = plot(gamma,get_balancing(L, GW, GS, gamma, CS=CS, alpha=0.)[0]/len(L[0]),'-',color=color_solar,lw=2)
	pp_wind = plot(gamma,get_balancing(L, GW, GS, gamma,  CS=CS, alpha=1.)[0]/len(L[0]),'-',color=color_wind,lw=2)
	
	pp = list(zeros(len(p_interval)))
	for j in arange(len(p_interval))[::-1]:
		pp_ = plot(gamma,res_load_sum_p[j],linespec[j],color='k')
		pp[j] = pp_
	pp = tuple(pp)
	
	pp_opt = plot(gamma,res_load_sum_opt,'k-')

	axvline(.2,ls='--',color='k')
	
				
	axis(xmin=0,xmax=amax(gamma),ymin=0,ymax=1)
	ylabel('Av. residual load [av.h.l.]')
	xlabel(r'RES power generation coefficient $\gamma$')
	
	pp = (pp_opt,pp_wind,pp_solar) + pp
	txtlabels = ('Optimal mix','Wind only','Solar only','1 pp contour','5 pp contour','25 pp contour')
	
	leg = legend(pp,txtlabels,loc='upper right',ncol=2);
	ltext  = leg.get_texts();
	setp(ltext, fontsize='small')    # the legend text fontsize
	
	save_file_name = 'plot_country_optimal_mix_vs_gamma_'+'CS_'+str(CS)+'.png'
	save_figure(figname=save_file_name, fig=gcf(), path='./', dpi=300)
