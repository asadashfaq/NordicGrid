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

#Custom modules
sys.path.append( './EuropeanGridR/' ) #This can be done in a more fancy way using __init__.py or some such.
from aures import solve as aures_solve #R's magic flow solver.
from auresc import Nodes
from colourutils import track_imports
from shortcuts import *
from SingleCountry import get_ISET_country_data, get_balancing

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


def test_new_flow_calc(lapse=None):

    N = Nodes(admat='./settings/admat_2011.txt',path='./data/',files=['ISET_NordicGrid_DE-N.npz', 'ISET_NordicGrid_DK-W.npz', 'ISET_NordicGrid_SE.npz','ISET_NordicGrid_DK-E.npz', 'ISET_NordicGrid_NO.npz'])
    
    N,F = aures_solve(N, lapse=lapse)
    
    return N, F

#
#   year, data_nodes, data_flows = get_nodes_and_flows_vs_year(year=linspace(1985,2053,5),lapse=365*24,path_nodes='./output_data/test/')
#
#   year, data_nodes, data_flows = get_nodes_and_flows_vs_year(year=linspace(1985,2053,21)[2:14],lapse=None,path_nodes='./output_data/test/')
#
def get_nodes_and_flows_vs_year(year=linspace(1985,2053,21), add_color=False,path_nodes='./output_data/',admat='./settings/admat_2011.txt',path_data='./data/',files_data=['ISET_NordicGrid_NO.npz', 'ISET_NordicGrid_SE.npz', 'ISET_NordicGrid_DK-W.npz', 'ISET_NordicGrid_DK-E.npz', 'ISET_NordicGrid_DE-N.npz'], path_settings='./settings/',copper=0, h0=None, b=1.0, lapse=None, squaremin=False, maxb=True):
    
    """ (almost) Updated to use EuropeanGridR."""

    #Initialize
    N = Nodes(admat=admat,path=path_data,files=files_data)
    
    Gamma = get_basepath_gamma(year)
    dirList = os.listdir(path_nodes)
    
    data_nodes, data_flows = [], []
    for i in arange(len(year)):
        print 'Year: {0:.0f}'.format(year[i])
        sys.stdout.flush()

        nodes_filename = 'nodes_year_'+str(year[i])
        flows_filename = 'flows_year_'+str(year[i])

        if nodes_filename+'.npz' in dirList:
            N._load_nodes_(nodes_filename+'.npz',full_load=True,path=path_nodes)
            F = load(path_nodes+flows_filename+'.npy')
            
            print 'Loaded nodes file: '+path_nodes+nodes_filename
            print 'Loaded flows file: '+path_nodes+flows_filename
            sys.stdout.flush()
            
            if add_color:
                #N.add_colored_import(F,lapse=lapse)
                N = track_imports(N,F,admat=admat,lapse=lapse)
                N.save_nodes(nodes_filename,path=path_nodes,)
                print 'Added color and resaved.'
                sys.stdout.flush()
            
        else:
            #Apply year to nodes (Dummy gamma's)
            N.set_gammas(Gamma.transpose()[i])
            N.set_alphas([1.,1.,1.,1.,.9])
            
            print 'Gamma: ' + str(Gamma.transpose()[i])
            sys.stdout.flush()
            
            #Calculate flows
            
            N,F = aures_solve(N, path=path_settings, copper=copper, h0=h0, b=b, lapse=lapse, squaremin=squaremin, maxb=maxb)
            
            #### DOES NOT WORK YET!!!
            if add_color:
                #N.add_colored_import(F,lapse=lapse)
                N = track_imports(N,F,admat=admat,lapse=lapse)
                print 'Added color.'
            
            N.save_nodes(nodes_filename,path=path_nodes)
            np.save(path_nodes+flows_filename,F)
        
        #append data
        data_nodes.append(deepcopy(N))
        data_flows.append(deepcopy(F))

    return year, data_nodes, data_flows

# plot_generation_summary_vs_year(year, data_nodes,lapse=365*24,datalabel='Test')
#
def plot_generation_summary_vs_year(year,data,lapse=50*24,datalabel=None):
    """ (almost) Updated to use EuropeanGridR."""
    
    if datalabel != None:
        datalabel = datalabel + '_'

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

        add_duplicate_yaxis(gcf(),unit_multiplier=mean(N[node_id].load)/1e3,label='[GW]')

        tight_layout(pad=.5)
        savename = 'plot_generation_summary_vs_year_Stacked_' + datalabel + region_name[node_id].replace(' ','_') + '.pdf'
        save_figure(savename)


##New version
#
#   plot_generation_summary_vs_year_2(year,data_nodes,lapse=None)
#
def plot_generation_summary_vs_year_2(year,data_nodes,node_id=[2,3],lapse=50*24,datalabel='Denmark'):
    """ (almost) Updated to use EuropeanGridR. NOTE: If balancing is shared, this function may need an update."""

    N = data_nodes[0]
    region_name = ['Norway','Sweden','Denmark West','Denmark East','Germany North']
    node_id = list(array(node_id,ndmin=1))
    
    gross_share, RES_local_av, balancing_av, curtailment_av, import_av, export_av = zeros(len(year)), zeros(len(year)), zeros(len(year)), zeros(len(year)), zeros(len(year)), zeros(len(year))
    
    for i in arange(len(year)):
        load_sum, gross_share_, RES_local_av_, import_share_, balancing_av_, import_av_, export_av_, curtailment_av_ = 0., 0., 0., 0., 0., 0., 0., 0.
        for node_id_ in node_id:
            node = data_nodes[i][node_id_]
            
            load_sum = load_sum + node.mean
            gross_share_ = gross_share_ + node.gamma*node.mean
            
            ## Power used locally ##
            
            ## Own RES used locally + import from friends.
            RES_local_av_ = RES_local_av_ + node.get_localRES()[:lapse].mean() + sum(node.get_colored_import()[find([id in node_id for id in arange(len(N))])],axis=0)[:lapse].mean()
            
            ## Own balancing
            ## NOTE: If balancing is shared, this may have to change.
            balancing_av_ = balancing_av_ + node.get_localBalancing()[:lapse].mean()
            
            ## Import from not-friends.
            import_av_ = import_av_ + sum(node.get_colored_import()[find([id not in node_id for id in arange(len(N))])],axis=0)[:lapse].mean()
            
            ## Power not used locally ##
            
            ## Own RES curtailed.
            curtailment_av_ = curtailment_av_ + node.curtailment[:lapse].mean()
            
            ## Own export - export to friends.
            export_av_ = export_av_ + (node.get_export() - sum(node.get_colored_import()[find([id in node_id for id in arange(len(N))])],axis=0))[:lapse].mean()
            
        #Normalize all
        gross_share[i] = gross_share_/load_sum
        RES_local_av[i] = RES_local_av_/load_sum #
        balancing_av[i] = balancing_av_/load_sum #
        curtailment_av[i] = curtailment_av_/load_sum #
        import_av[i] = import_av_/load_sum
        export_av[i] = export_av_/load_sum

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1);figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])

    ax = axes()

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

    axis(ymin=0, ymax =1.55, xmin=amin(year), xmax=amax(year))
    yticks(arange(0,1.55,.25))
    xlabel('Reference year')
    ylabel('Power [av.l.h.]')

    pp = [pp_balancing,pp_VRES,pp_import,pp_export,pp_curtailment]
    pp_text = ['Balancing','VRES','Import','Export','Curtailment']
    leg = legend(pp[::-1],pp_text[::-1],title=datalabel,loc='lower left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    #Annotations
    year_export_1pp = interp(0.01,export_av,year)
    if year_export_1pp<amax(year):
        gross_export_1pp = interp(year_export_1pp,year,gross_share)
        ax.annotate('{0:.0f} ({1:.0f}%):\nExport: 0.01 av.l.h.'.format(year_export_1pp,100*gross_export_1pp), xy=(year_export_1pp, 1),  xycoords='data',xytext=(-90, 12), textcoords='offset points',bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),fontsize='small')

    year_import_1pp = interp(0.01,import_av,year)
    if year_import_1pp<amax(year):
        gross_import_1pp = interp(year_import_1pp,year,gross_share)
        ax.annotate('{0:.0f} ({1:.0f}%):\nImport: 0.01 av.l.h.'.format(year_import_1pp,100*gross_import_1pp), xy=(year_import_1pp, 1.0),  xycoords='data',xytext=(-80, 40), textcoords='offset points',bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),fontsize='small')

    year_curt_1pp = interp(0.01,curtailment_av,year)
    if year_curt_1pp<amax(year):
        gross_curt_1pp = interp(year_curt_1pp,year,gross_share)
        offset_curt_1pp = 1 + interp(year_curt_1pp,year,export_av)
        ax.annotate('{0:.0f} ({1:.0f}%):\nCurtailment: 0.01 av.l.h.'.format(year_curt_1pp,100*gross_curt_1pp), xy=(year_curt_1pp, offset_curt_1pp),  xycoords='data',xytext=(0, 45), textcoords='offset points',bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),fontsize='small')


    add_duplicate_yaxis(gcf(),unit_multiplier=load_sum/1e3,label='[GW (2007)]')

    tight_layout(pad=.5)
    savename = 'plot_generation_summary_vs_year_Stacked_2_' + str(datalabel) + '.pdf'
    save_figure(savename)

#    
# plot_colored_import_export(year, data_nodes, datalabel='2011')
# plot_colored_import_export(year_cu, data_nodes_cu, datalabel='copper')
#    
def plot_colored_import_export(year, data, colors=colors_countries, lapse=None, datalabel=''):
    """ (almost) Updated to use EuropeanGridR. NOTE: If balancing is shared, this function may need an update."""
    
    region_names = ['NO','SE','DK-W','DK-E','DE-N']
    Nodes = data[0]

    if lapse==None:
        lapse=Nodes[0].mismatch.shape[0]

    for node_id in arange(len(Nodes)):

        ## Calculate average import
        colored_import_av = zeros((len(year),len(Nodes)))
        for	i in arange(len(year)):
            colored_import_av[i] = mean(data[i][node_id].get_colored_import().transpose()[:lapse], axis=0)/data[i][node_id].mean

        ## Calculate export
        colored_export_av = zeros((len(year),len(Nodes)))
        for	i in arange(len(year)):
            for j in arange(len(Nodes)):
                colored_export_av[i][j] = mean(data[i][j].get_colored_import()[node_id][:lapse])/data[i][node_id].mean

        close(1); figure(1); clf()
        
        ## Set non-standard figure size.
        gcf().set_size_inches([5.25,3.5])

        ## Plot import ##
        subplot(211)

        pp = []; pp_text=[]
        
        ## Plot first country.
        fill_between(year,cumsum(colored_import_av,axis=1).transpose()[0],color=colors[0],edgecolor='k',lw=.5)
        
        ## Add legend entry.
        pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[0]))
        pp_text.append(region_names[0])
        
        ## Plot all other countries.
        for i in arange(1,len(Nodes)):
            fill_between(year,cumsum(colored_import_av,axis=1).transpose()[i],cumsum(colored_import_av,axis=1).transpose()[i-1],label=region_names[i],color=colors[i],edgecolor='k',lw=.5)
            
            ## Add legend entry.
            pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[i]))
            pp_text.append(region_names[i])

        ## Remove own id from legend.
        pp.pop(node_id); pp_text.pop(node_id)
        
        ## Plot legend (in reverse order.)
        leg = legend(tuple(pp[::-1]),tuple(pp_text[::-1]),loc='upper left',title=region_names[node_id]+' - import')
        
        ## Set ticks
        x = amax(cumsum(colored_import_av,axis=1).transpose()[len(Nodes)-1])
        dx = 10**floor(log10(x))
        if 5*dx<x:
            yticks(arange(0,2*x,2.5*dx))
            ymax_=1.05*arange(0,2*x,2.5*dx)[argmax(arange(0,2*x,2.5*dx)>=x)]
        else:
            yticks(arange(0,2*x,dx))
            ymax_=1.05*arange(0,2*x,dx)[argmax(arange(0,2*x,dx)>=x)]
            
        axis(ymin=0, ymax=ymax_, xmin=1995, xmax=amax(year))
        xlabel('Reference year')
        ylabel('Power [av.l.h.]')
        
        ## Include real units on secondary y-axis.
        add_duplicate_yaxis(gcf(),unit_multiplier=Nodes[node_id].mean,label='[MW]',tickFormatStr='%.0f')

        ## Plot export ##
        subplot(212)
        
        pp = []; pp_text=[]
        
        ## Plot first country.
        fill_between(year,cumsum(colored_export_av,axis=1).transpose()[0],color=colors[0],edgecolor='k',lw=.5)
        
        ## Add legend entry.
        pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[0]))
        pp_text.append(region_names[0])
        
        ## Plot all other countries.
        for i in arange(1,len(Nodes)):
            fill_between(year,cumsum(colored_export_av,axis=1).transpose()[i],cumsum(colored_export_av,axis=1).transpose()[i-1],label=region_names[i],color=colors[i],edgecolor='k',lw=.5)
            
            ## Add legend entry.
            pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[i]))
            pp_text.append(region_names[i])

        ## Remove own id from legend.
        pp.pop(node_id); pp_text.pop(node_id)
        
        ## Plot legend (in reverse order.)
        leg = legend(tuple(pp[::-1]),tuple(pp_text[::-1]),loc='upper left',title=region_names[node_id]+' - export')
        
        ## Set ticks
        x = amax(cumsum(colored_export_av,axis=1).transpose()[len(Nodes)-1])
        dx = 10**floor(log10(x))
        if 5*dx<x:
            yticks(arange(0,2*x,2.5*dx))
            ymax_=1.05*arange(0,2*x,2.5*dx)[argmax(arange(0,2*x,2.5*dx)>=x)]
        else:
            yticks(arange(0,2*x,dx))
            ymax_=1.05*arange(0,2*x,dx)[argmax(arange(0,2*x,dx)>=x)]
            
        axis(ymin=0, ymax=ymax_, xmin=1995, xmax=amax(year))
        xlabel('Reference year')
        ylabel('Power [av.l.h.]')
        
        ## Include real units on secondary y-axis.
        add_duplicate_yaxis(gcf(),unit_multiplier=Nodes[node_id].mean,label='[MW]',tickFormatStr='%.0f')
        
        tight_layout(pad=.5)
        savename = 'plot_colored_import_export_' + region_names[node_id] + '_' + datalabel + '.pdf'
        save_figure(savename)

#    
# plot_colored_import_export(year, data_nodes, lapse=None, datalabel='Test_2')
#    
# FOR NINGLING, MAY 2012:
# plot_colored_import_export(year, data_nodes, datalabel='2011_mod')
# plot_colored_import_export(year_ref, data_nodes_ref, datalabel='2011_ref')
#
def plot_colored_import_export_alt(year, data, colors=colors_countries, lapse=None, datalabel=''):
    region_names = ['NO','SE','DK-W','DK-E','DE-N']
    Nodes = data[0]

    if lapse==None:
        lapse=Nodes[0].mismatch.shape[0]

    for node_id in arange(len(Nodes)):

        ## Calculate average import
        colored_import_av = zeros((len(year),len(Nodes)))
        for	i in arange(len(year)):
            colored_import_av[i] = mean(data[i][node_id].get_colored_import().transpose()[:lapse], axis=0)/data[i][node_id].mean

        ## Calculate export
        colored_export_av = zeros((len(year),len(Nodes)))
        for	i in arange(len(year)):
            for j in arange(len(Nodes)):
                colored_export_av[i][j] = mean(data[i][j].get_colored_import()[node_id])/data[i][node_id].mean

        ## Set plot options	##
        matplotlib.rcParams['font.size'] = 10

        figure(1); clf()

        gcf().set_dpi(300)
        gcf().set_size_inches([5.25,3.5])

        ## Plot import figure ##
        subplot(211)

        envelope = cumsum(colored_import_av,axis=1).transpose()[-1]

        pp = []; pp_text=[]
        fill_between(year,cumsum(colored_import_av,axis=1).transpose()[0]/envelope,color=colors[0],edgecolor='k',lw=.5)
        pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[0]))
        pp_text.append(region_names[0])
        for i in arange(1,len(Nodes)):
            fill_between(year,cumsum(colored_import_av,axis=1).transpose()[i]/envelope,cumsum(colored_import_av,axis=1).transpose()[i-1]/envelope,label=region_names[i],color=colors[i],edgecolor='k',lw=.5)
            
            ## Add legend entry.
            pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[i]))
            pp_text.append(region_names[i])

        ## Plot lines
        plot(year,envelope,'k-',lw=3)

        pp.pop(node_id); pp_text.pop(node_id)
        leg = legend(tuple(pp[::-1]),tuple(pp_text[::-1]),loc='upper left',title=region_names[node_id]+' - import')
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        #x = amax(cumsum(colored_import_av,axis=1).transpose()[len(Nodes)-1])
        #dx = 10**floor(log10(x))
        #if 5*dx<x:
        #    yticks(arange(0,2*x,2.5*dx))
        #    ymax_=1.05*arange(0,2*x,2.5*dx)[argmax(arange(0,2*x,2.5*dx)>=x)]
        #else:
        #    yticks(arange(0,2*x,dx))
        #    ymax_=1.05*arange(0,2*x,dx)[argmax(arange(0,2*x,dx)>=x)]
            
        axis(ymin=0, xmin=1995, xmax=amax(year))
        xlabel('Reference year')
        ylabel('Power [av.l.h.]')
        
        #majorFormatter = matplotlib.ticker.ScalarFormatter(useMathText=True)#FormatStrFormatter('%.1f')
        #majorFormatter.set_powerlimits((0, 0))
        #gca().yaxis.set_major_formatter(majorFormatter)
        
        add_duplicate_yaxis(gcf(),unit_multiplier=Nodes[node_id].mean,label='[MW]',tickFormatStr='%.0f')

        subplot(212)
        
        pp = []; pp_text=[]
        fill_between(year,cumsum(colored_export_av,axis=1).transpose()[0],color=colors[0],edgecolor='k',lw=.5)
        pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[0]))
        pp_text.append(region_names[0])
        for i in arange(1,len(Nodes)):
            fill_between(year,cumsum(colored_export_av,axis=1).transpose()[i],cumsum(colored_export_av,axis=1).transpose()[i-1],label=region_names[i],color=colors[i],edgecolor='k',lw=.5)
            pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[i]))
            pp_text.append(region_names[i])

        pp.pop(node_id); pp_text.pop(node_id)
        leg = legend(tuple(pp[::-1]),tuple(pp_text[::-1]),loc='upper left',title=region_names[node_id]+' - export')
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        x = amax(cumsum(colored_export_av,axis=1).transpose()[len(Nodes)-1])
        dx = 10**floor(log10(x))
        if 5*dx<x:
            yticks(arange(0,2*x,2.5*dx))
            ymax_=1.05*arange(0,2*x,2.5*dx)[argmax(arange(0,2*x,2.5*dx)>=x)]
        else:
            yticks(arange(0,2*x,dx))
            ymax_=1.05*arange(0,2*x,dx)[argmax(arange(0,2*x,dx)>=x)]
            
        axis(ymin=0, ymax=ymax_, xmin=1995, xmax=amax(year))
        xlabel('Reference year')
        ylabel('Power [av.l.h.]')
        
        #majorFormatter = matplotlib.ticker.ScalarFormatter(useMathText=True)#FormatStrFormatter('%.1f')
        #majorFormatter.set_powerlimits((0, 0))
        #gca().yaxis.set_major_formatter(majorFormatter)
        
        add_duplicate_yaxis(gcf(),unit_multiplier=mean(Nodes[node_id].load),label='[MW]',tickFormatStr='%.0f')
        
        tight_layout(pad=.5)
        savename = 'plot_generation_summary_vs_year_Colored_import_alt_' + region_names[node_id] + datalabel + '.pdf'
        save_figure(savename)

##
# plot_gross_net_total_share(year,data_nodes,node_id=[2,3],label='2011')
# plot_gross_net_total_share(year_cu,data_nodes_cu,node_id=[2,3],label='copper')
#
def plot_gross_net_total_share(year,data_nodes,node_id=2,lapse=None,ROI=[.2,.5],label=None):
    
    node_id = list(array(node_id,ndmin=1))
    
    gross_share, net_res_load, res_load, net_plus_export = zeros(len(year)), zeros(len(year)), zeros(len(year)), zeros(len(year))
    
    for i in arange(len(year)):
        gross_share_, net_share_, import_share_, localBalancing_, net_plus_export_, load_sum = 0., 0., 0., 0., 0., 0.
        for node_id_ in node_id:
            gross_share_ = gross_share_ + data_nodes[i][node_id_].get_wind()[:lapse].mean() + data_nodes[i][node_id_].get_solar()[:lapse].mean()
            net_share_ = net_share_ + data_nodes[i][node_id_].get_localRES()[:lapse].mean()
            
            localBalancing_ = localBalancing_ + data_nodes[i][node_id_].get_localBalancing()[:lapse].mean()
            load_sum = load_sum + data_nodes[i][node_id_].mean
            
            net_plus_export_ = net_plus_export_ + data_nodes[i][node_id_].get_wind()[:lapse].mean() + data_nodes[i][node_id_].get_solar()[:lapse].mean() - data_nodes[i][node_id_].curtailment[:lapse].mean()
        
        gross_share[i] = gross_share_/load_sum
        net_res_load[i] = 1 - net_share_/load_sum
        res_load[i] = localBalancing_/load_sum
        net_plus_export[i] = 1- net_plus_export_/load_sum
    
    gross_share, net_res_load, res_load, net_plus_export = concatenate([[0],gross_share]), concatenate([[1],net_res_load]), concatenate([[1],res_load]), concatenate([[1],net_plus_export])
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1);figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    #Reference line, wind only, DK as one region:
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('DK')
    gamma = interp(linspace(0,amax(gross_share),111),gross_share,gross_share)
    
    plot(gamma,1-gamma,ls='-',color=(.5,.5,.5),label='Local limit') #Full integration of local VRES, no import.
    
    plot(gamma,get_balancing(L, Gw, Gs, gamma, alpha=.8)[0]/len(L),'-',color=(0.53,0.73,0.37),lw=2,label='80/20 mix') #"Optimal mix", actually alpha_w=0.8
    plot(gamma,get_balancing(L, Gw, Gs, gamma, alpha=1.)[0]/len(L),'-',color=color_wind,lw=2,label='Wind-only') #Wind-only
    
    plot(gross_share,net_res_load,ls=':',color='k',label='Residual load (net)') #Residual load before import. 
    plot(gross_share,res_load,ls='-',color='k',lw=1.5,label='Residual load') #Local residual load after import.
    plot(gross_share,net_plus_export,ls='--',color='k',label='Net + export') #Global use of local VRES.
    
    axhline(1.0,color='k',ls='--')
    
    #Highlight ROI
    if not ROI==None:
        fill_betweenx([0,10],ROI[1]*ones(2),ROI[0],lw=0,color='k',alpha=.1)
    
    axis(xmin=0,ymin=0,ymax=1.0)
    
    xlabel(r'Gross share of electricity demand $\gamma_{DK}$')
    ylabel(r'Av. residual load [av.h.l.]')
    
    leg = legend(loc='upper right')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
        
    add_duplicate_yaxis(gcf(),unit_multiplier=load_sum/1e3,label='[GW]',tickFormatStr='%.1f')

    tight_layout()
    save_figure('plot_gross_net_total_share_'+label+'.pdf')
    

def plot_flow_statistics(data_flows,i_year = 10,i_link = 4,savelabel='test'):

    
    flows = data_flows[i_year]
    
    close(1); figure(1); clf();
    
    bins = linspace(1.1*amin(flows[i_link]),1.1*amax(flows[i_link]),100)
    h_data, bins, patch = hist(flows[i_link],bins,normed=True)
    
    
    axis(xmin=1.1*amin(flows[i_link]),xmax=1.1*amax(flows[i_link]),ymin=0, ymax=1.1*amax(h_data[find(abs(bins[1:])>=2*abs(bins[0]-bins[1]))]))

    save_figure('plot_flow_statistics_'+savelabel+'.pdf')

def plot_surplus_statistics(year,data_nodes,i_year = [10],i_node = 2,savelabel='test',n_bins=21,xmax=4500,ymax=500,xline=[3790,2790]):

    
    
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])
    
    surplus = []
    legendlabels = []
    for i in i_year:
        surplus.append(get_positive(data_nodes[i][i_node].mismatch))
        legendlabels.append('{0:0.0f}'.format(year[i]))
    
    surplus = array(surplus).transpose()
    
    bins = linspace(0,1.01*amax(surplus),n_bins)
    weights = 365*24*ones(surplus.shape)/len(surplus)#/abs(bins[0]-bins[1])
    
    h_data, bins, patch = hist(surplus,bins,normed=False,weights=weights,label=legendlabels)
    
    for x in xline:
        axvline(x,ls='--',color='k')
        text(x+abs(bins[0]-bins[1])*.25,0.8*ymax,'{0:0.0f} MW'.format(x))
    
    xlabel('Power [MW]')
    ylabel('Frequency [hours/yr]')
    axis(xmin=0,xmax=xmax,ymin=0, ymax=ymax)
    

    leg = legend(title='Surplus',loc='upper left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')  

    tight_layout()
    save_figure('plot_surplus_statistics_'+savelabel+'.pdf')

    
# FOR NINGLING, MAY 2012:
#
# DK-W -> DE-N
# plot_compare_flow_statistics(year,data_flows,data_flows_ref,i_year = i,i_link = 5,linkname='DK-W -> DE-N',savelabel='DKW2DEN_ref_mod',legendlabels=['Modified layout','2011 layout'])
#
# DK-W -> DK-E
# plot_compare_flow_statistics(year,data_flows,data_flows_ref,i_year = i,i_link = 4,linkname='DK-W -> DK-E',savelabel='DKW2DKE_ref_mod',legendlabels=['Modified layout','2011 layout'])
#
# NO -> DK-W
# plot_compare_flow_statistics(year,data_flows,data_flows_ref,i_year = i,i_link = 1,linkname='NO -> DK-W',savelabel='NO2DKW_ref_mod',legendlabels=['Modified layout','2011 layout'])
#
# SE -> DK-W
# plot_compare_flow_statistics(year,data_flows,data_flows_ref,i_year = i,i_link = 2,linkname='SE -> DK-W',savelabel='SE2DKW_ref_mod',legendlabels=['Modified layout','2011 layout'])
#
#
def plot_compare_flow_statistics(year,data_flows_1,data_flows_2,i_year = 15,i_link = 5,legendlabels=['1','2'],linkname='DK-W -> DE-N',savelabel='test',N_bins=21,ymax=510,xmin=-1010,xmax=1610):
      
    flows_1 = data_flows_1[i_year]
    flows_2 = data_flows_2[i_year]
    
    flows = array([flows_1[i_link],flows_2[i_link]]).transpose()
    
    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])
    
    #bins = linspace(1.01*amin(flows),1.01*amax(flows),N_bins)
    bins = linspace(xmin,xmax,N_bins)
    weights = 365*24*ones(flows.shape)/len(flows_1[i_link])#/abs(bins[0]-bins[1])
    h_data, bins, patch = hist(flows,bins,weights=weights,normed=False,label=legendlabels)
    
    if ymax==None:
        ymax = 1.5*amax(h_data[0][amin(find(bins>=abs(bins[0]-bins[1])))])
    axis(xmin=xmin,xmax=xmax,ymin=0, ymax=ymax)
    # axis(xmin=1.1*amin(flows_1[i_link]),xmax=1.1*amax(flows_1[i_link]),ymin=0, ymax=1.1*amax(h_data[0][find(abs(bins[1:])>=2*abs(bins[0]-bins[1]))]))

    xlabel('Power [MW]')
    ylabel('Frequency [hours/yr]')
    
    leg = legend(title='{0:0.0f}:'.format(year[i_year]) + '\n' + linkname)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    tight_layout()
    save_figure('plot_compare_flow_statistics_'+'{0:0.0f}_'.format(year[i_year])+savelabel+'.pdf')
      
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