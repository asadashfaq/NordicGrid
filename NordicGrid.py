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
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Custom modules
sys.path.append( './zdcpf/' ) #This can be done in a more fancy way using __init__.py or some such.
from zdcpf import *
from shortcuts import *

#Colors. To be placed somewhere else
colors_countries = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers.
color_balancing = '#2B2825' #Conspicuous Creep from COLOURlovers.
color_RES = '#CACF43'
color_import = '#FCF8BC'
color_export = '#0B8C8F'
color_curtailment = '#D6156C'

color_load = (.9,.5,.5,.3)
color_slow = (0.2,0.,0.,1)
color_medium = (.5,0,0,1)
color_fast = (0.8,0,0,1)
color_wind = (0.5,0.7,1.,1)
color_solar = (1.,0.8,0,1)
color_edge = (.2,.2,.2)
#color_RES = (0.4,.65,0.25)


#This function will later be replaced by some fancy save/load thing.
#
# year, data_nodes, data_flows = get_nodes_and_flows_vs_year(lapse=50*24)
#
#
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

#
# plot_generation_summary_vs_year(year,data_nodes,lapse=50*24)
#
def plot_generation_summary_vs_year(year,data,lapse=50*24):

    N = data[0]
    region_name = ['Norway','Sweden','Denmark West','Denmark East','Germany North']

    for node_id in arange(len(N)):
        #Calculate averages to be displayed
        balancing_av, import_av, RES_local_av, curtailment_av, export_av = zeros(year.shape), zeros(year.shape), zeros(year.shape), zeros(year.shape), zeros(year.shape)
        for i in arange(len(year)):
            
            #Power used locally
            balancing_av[i] = data[i][node_id].get_localBalancing()[:lapse].mean()/data[i][node_id].mean  #Wrong. Balancing can be exported too!!!!
            RES_local_av[i] = data[i][node_id].get_localRES()[:lapse].mean()/data[i][node_id].mean
            import_av[i] = data[i][node_id].get_import()[:lapse].mean()/data[i][node_id].mean
            
            #Power not used locally
            export_av[i] = data[i][node_id].get_export()[:lapse].mean()/data[i][node_id].mean
            curtailment_av[i] = data[i][node_id].curtailment[:lapse].mean()/data[i][node_id].mean



        #Set plot options	
        matplotlib.rcParams['font.size'] = 10

        close(1);figure(1); clf()

        gcf().set_dpi(300)
        gcf().set_size_inches([5.25,3.5])

        ax1 = axes()

        pp_balancing = fill_between(year,balancing_av,color=color_balancing,edgecolor='k',lw=.5)
        pp_VRES = fill_between(year,RES_local_av+balancing_av,balancing_av,color=color_RES,edgecolor='k',lw=.5)
        pp_import = fill_between(year,import_av+RES_local_av+balancing_av,RES_local_av+balancing_av,color=color_import,edgecolor='k',lw=.5)
        pp_export = fill_between(year,export_av+import_av+RES_local_av+balancing_av,import_av+RES_local_av+balancing_av,color=color_export,edgecolor='k',lw=.5)
        pp_curtailment = fill_between(year,curtailment_av+export_av+import_av+RES_local_av+balancing_av,export_av+import_av+RES_local_av+balancing_av,color=color_curtailment,edgecolor='k',lw=.5)

        pp_balancing = Rectangle((0, 0), 1, 1, facecolor=color_balancing)
        pp_VRES = Rectangle((0, 0), 1, 1, facecolor=color_RES)
        pp_import = Rectangle((0, 0), 1, 1, facecolor=color_import)
        pp_export = Rectangle((0, 0), 1, 1, facecolor=color_export)
        pp_curtailment = Rectangle((0, 0), 1, 1, facecolor=color_curtailment)

        axis(ymin=0, ymax =2.05, xmin=amin(year), xmax=amax(year))
        yticks(arange(0,2.05,.5))
        xlabel('Reference year')
        ylabel('Power [av.l.h.]')

        pp = [pp_balancing,pp_VRES,pp_import,pp_export,pp_curtailment]
        pp_text = ['Balancing','VRES','Import','Export','Curtailment']
        leg = legend(pp[::-1],pp_text[::-1],title=region_name[node_id],loc='upper left')
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        divider = make_axes_locatable(plt.gca())
        ax2 = divider.append_axes("right", "0%", pad="0%")
        #ax2 = axes(ax1.get_position())

        #ax2 = twinx()
        ax2.yaxis.set_ticks_position('right')
        ax2.yaxis.set_label_position('right')
        axis(ymin=0,ymax=2.05)
        ylabel('[GW]')
        xticks([0],[''])
        yticks(ax1.get_yticks(),around(ax1.get_yticks()*N[node_id].mean/1e3,1))

        tight_layout(pad=.3)
        savename = 'plot_generation_summary_vs_year_Stacked_' + N[node_id].name.replace(' ','_') + '.png'
        save_figure(savename)
    
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