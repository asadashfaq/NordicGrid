#! /usr/bin/env python
from pylab import *
from scipy import *
import scipy.optimize as optimize
from numpy import concatenate as conc
import numpy as np
from cvxopt import matrix,solvers
from time import time
import os
from copy import deepcopy
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

class node:
    def __init__(self,fileName,ID,setup,name):
        self.id = ID
        data = map(np.double,np.load(fileName))
        self.gamma = float(setup[ID][0])
        self.alpha = float(setup[ID][1])  #Alpha should be expanded to a vector.  completalpha() can be applied in update()
        self.load = data[0]
        self.nhours = len(self.load)
        self.normwind = data[1]
        self.normsolar = data[2]
        self.mean = np.mean(self.load)
        self.balancing = np.zeros(self.nhours)
        self.curtailment = np.zeros(self.nhours)
        self.baltot = []
        self.curtot = []
        self.name = name
        self.mismatch = None
        self.colored_import = None #Set using self.set_colored_i_import()
		
        self._update_()
    
    def _update_(self):
        self.mismatch=(self.getwind()+self.getsolar())-self.load
    
    def getimport(self):
        """Returns import power time series in units of MW."""
        return get_positive(get_positive(-self.mismatch) - self.balancing) #Balancing is exported if it exceeds the local residual load.
        
    def getexport(self):	
        """Returns export power time series in units of MW."""
        return get_positive(self.mismatch) - self.curtailment + get_positive(self.balancing - get_positive(-self.mismatch))
	
    def getlocalRES(self):
        """Returns the local use of RES power time series in units of MW."""
        return self.getwind() + self.getsolar() - self.curtailment  - self.getexport()
		
    def getlocalBalancing(self):
        """Returns the local use of balancing power time series in units of MW."""
        return get_positive(-self.mismatch) - self.getimport()
		
    def getwind(self):
        """Returns wind power time series in units of MW."""
        return self.mean*self.gamma*self.alpha*self.normwind
	
    def getsolar(self):
        """Returns solar power time series in units of MW."""
        return self.mean*self.gamma*(1.-self.alpha)*self.normsolar
			
    def setgamma(self,gamma):
        self.gamma = gamma
        self._update_()
    
    def setalpha(self,alpha):
        self.alpha=alpha
        self._update_()

    def set_colored_import_i(self,i,colored_import_i):
        if self.colored_import == None:
            self.colored_import = zeros((len(colored_import_i),len(self.load)))
			
        self.colored_import.transpose()[i] = colored_import_i



###########        This function generates matrices, or     ###################
###########        at least their basic structure.          ###################

def generatemat(constraintsfile='lala.txt'):
    K=np.genfromtxt("incidence.txt",delimiter=",",dtype='d')
    Nnodes=np.size(K,0)
    Nlinks=np.size(K,1)
    # These numbers include the dummy node and link
    # With this info, we create the P matrix, sized
    P1=np.eye(Nlinks+2*Nnodes)  # because a row is needed for each flow, and two for each node
    P=conc((P1[:Nlinks],P1[-2*Nnodes:]*1e-6))  # and the bal/cur part has dif. coeffs
    # Then we make the q vector, whose values will be changed all the time
    q=np.zeros(Nlinks+2*Nnodes)  # q has the same size and structure as the solution 'x'
    # Then we build the equality constraint matrix A
    # The stucture is more or less [ K | -I | I ]
    A1=conc((K,-np.eye(Nnodes)),axis=1)
    A=conc((A1,np.eye(Nnodes)),axis=1)
    # See documentation for why first row is cut
    A=np.delete(A,np.s_[0],axis=0)
    # b vector will be defined by the mismatches, in MAIN
    # Finally, the inequality matrix and vector, G and h.
    # Refer to doc to understand what the hell I'm doing, as I build G...
    g1=np.eye(Nlinks)
    G1=g1
    for i in range(Nlinks-1):
        i+=i
        G1=np.insert(G1,i+1,-G1[i],axis=0)
    G1=conc((G1,-G1[-1:]))
    # G1 is ready, now we make G2
    G2=np.zeros((2*Nlinks,2*Nnodes))
    # G3 is built as [ 0 | -I | 0 ]
    g3=conc((np.zeros((Nnodes,Nlinks)),-np.eye(Nnodes)),axis=1)
    G3=conc((g3,np.zeros((Nnodes,Nnodes))),axis=1)
    g4=np.eye(Nnodes)
    G4=g4
    for i in range(Nnodes-1):
        i+=i
        G4=np.insert(G4,i+1,-G4[i],axis=0)
    G4=conc((G4,-G4[-1:]))
    G5=conc((np.zeros((2*Nnodes,Nlinks+Nnodes)),G4),axis=1)
    G=conc((G1,G2),axis=1)
    G=conc((G,G3))
    G=conc((G,G5))
    # That was crazy! Now, the h vector is partly made from the constraints.txt file
    h=np.genfromtxt("constraints.txt",dtype='d')
    # but also added some extra zeros for the rest of the constraints.
    h=np.append(h,np.zeros(3*Nnodes))
    # And that's it!
    return P,q,G,h,A


def dcpowerflow(P,q,G,h,A,b):
    sol=solvers.qp(P,q,G,h,A,b)
    return sol['x']

# This generates the Nodal time series from the binary files and the setup.txt file.
def makenodes():
    setup=np.genfromtxt('setupnodes.txt',delimiter=',')
    N=node('N.npy',0,setup,'Norway')
    S=node('S.npy',1,setup,'Sweden')
    DKv=node('DKW.npy',2,setup,'Denmark West')
    DKo=node('DKE.npy',3,setup,'Denmark East')
    D=node('DN.npy',4,setup,'Germany North')
    return [N,S,DKv,DKo,D]

def setconstraints(h,Nlinks,value):
    # default 'value' is np.genfromtxt("constraints.txt",delimiter=',',dtype='d')
    # 'value' can be a single double or a vector
    if np.size(value)==1:
        h[2:2*Nlinks+2]=value
        return h
    if np.size(value)!=2*Nlinks:
        print "Wrong constraint vector size. ", np.size(value,0)," were received, ",2*Nlinks," were expected."
        return h
    h[2:2*Nlinks+2]=value
    return h

def setgammas(Nodes,value):
    # to change a single node's gamma, just write XX.setgamma(yy)
    # 'value' can be a single number or a vector
    if np.size(value)==1:
        for i in Nodes: i.setgamma(value)
        return Nodes
    if np.size(value)!=np.size(Nodes):
        print "Wrong gamma vector size. ", np.size(value,0)," were received, ",np.size(Nodes)," were expected."
        return Nodes
    for i in Nodes:
        i.setgamma(value[i.id])
    return Nodes
    
def setalphas(Nodes,value):
    # to change a single node's alpha, just write XX.setalpha(yy)
    # 'value' can be a single number or a vector
    if np.size(value)==1:
        for i in Nodes: i.setalpha(value)
        return Nodes
    if np.size(value)!=np.size(Nodes):
        print "Wrong alpha vector size. ", np.size(value,0)," were received, ",np.size(Nodes)," were expected."
        return Nodes
    for i in Nodes:
        i.setalpha(value[i.id])
    return Nodes  

def runtimeseries(Nodes,F,P,q,G,h,A,coop,lapse=100*24):
    """lapse: 0-70128  """
	
    if lapse==None:
        lapse=Nodes[0].mismatch.shape[0]
	
    Nlinks=np.size(F,0)
    Nnodes=np.size(A,0)
    start=time()
    b=matrix([[0,0,0,0,0]],tc='d')
    P_b=P[Nlinks+2:Nlinks+Nnodes+2,:]*1e6
    for t in range(lapse):
        for i in Nodes:
            b[i.id]=i.mismatch[t]
            # from default, both curtailment and balancing have a minimum of 0.
            # in order to prevent export of curtailment, max curtailment is set to b
            h[2*Nlinks+Nnodes+5+2*i.id]=0
            if b[i.id]>0:
                h[2*Nlinks+Nnodes+5+2*i.id]=b[i.id]
        # then, we set the values of q_b and q_r for bal and cur, according to doc.
        # for Gorm's inequalities, we need f,L,delta
        f=P[0,0]
        L=Nnodes-1
        d=np.array(b)
        excess=np.dot(d.T,d>0)[0][0]
        deficit=abs(np.dot(d.T,d<0)[0][0])
        delta=min(excess,deficit)
        q_r=L*f*2*delta*0.5
        q_b=L*f*2*delta+q_r*(1.5)
        q[Nlinks+2:Nlinks+Nnodes+2]=q_b
        q[Nlinks+Nnodes+2:]=q_r
        if coop==1:
            P[Nlinks+2:Nlinks+Nnodes+2,:]=P_b*L*f*deficit*.99
        opt=dcpowerflow(P,q,G,h,A,b)   ########### Save relevant solution as flows
        for j in range(Nlinks):
            F[j][t]=opt[j+1]           
        for k in Nodes:                ########### Save balancing at each node
            k.balancing[t]=opt[2+Nlinks+k.id]
            k.curtailment[t]=opt[3+Nlinks+Nnodes+k.id]  
        end=time()
        if (np.mod(t,547)==0) and t>0:
            print "Elapsed time is ",round(end-start)," seconds. t = ",t," out of ",lapse
            sys.stdout.flush()
    for i in Nodes:
        i.baltot.append(sum(i.balancing))
        i.curtot.append(sum(i.curtailment))
    end=time()
    print "Calculation took ",round(end-start)," seconds."
    sys.stdout.flush()
    return Nodes,F
          
#
# Nodes, F, FF, C = main(lapse=100*24)
#
def main(lapse=None):

	#####################################################################################
	###################################        MAIN       ###############################
	#####################################################################################
	Nodes=makenodes()
	solvers.options['show_progress']=False
	coop=0
	###################### Building CVXOPT Matrices ###############################
	###############################################################################
	P,q,G,h,A = generatemat()
	Nnodes=np.size(np.genfromtxt("incidence.txt",delimiter=','),0)-1
	Nlinks=np.size(np.genfromtxt("incidence.txt",delimiter=','),1)-1
	F=np.zeros((Nlinks,70128))
	# they must also be in cvxopt native "matrix" format (not the same as numpy's!)
	P=matrix(P,tc='d')
	q=matrix(q,tc='d')
	G=matrix(G,tc='d')
	h=matrix(h,tc='d')
	A=matrix(A,tc='d')
	##################### Running through the time series #########################
	###############################################################################
	
	Nodes = setgammas(Nodes,[1,1,1,1,1]) #setgammas(Nodes,[.1,.13,.63,.31,.4])
	Nodes = setalphas(Nodes,[1.,1.,1.,1.,.9])
	
	Nodes,F=runtimeseries(Nodes,F,P,q,G,h,A,coop=0,lapse=lapse)
	
	t = 100
	export_ = array([Nodes[i].getexport()[t] for i in arange(len(Nodes))])
	import_ = array([Nodes[i].getimport()[t] for i in arange(len(Nodes))])

	print 'Flow:: ', floor(F.transpose()[t])
	print 'Import: ', floor(import_)
	print 'Export: ', floor(export_)
	
	FF, C = get_colored_flow(F.transpose()[t], export_)
	
	print C
	
	
	#####################      An example sweeping        ##########################
	# lets run the series three times, for three different global gammas
	#gammas=[0.25,0.5,0.75]
	#for i in gammas:
	#   Nodes=setgammas(Nodes,i)
	#    Nodes,F=runtimeseries(Nodes,F,P,q,G,h,A,coop=0)

	return Nodes, F, FF, C

#Test lal delete me

#####
# Plots
#####

#
# plot_ts(Nodes[2],F)
#
def plot_ts(node,F):

    flow_DKW = F[1] + F[2] - F[4] - F[5]
    flow_in_DKW = flow_DKW*(flow_DKW>0)
    flow_out_DKW = -flow_DKW*(flow_DKW<0)

    t=arange(len(node.mismatch))

    figure(1); clf()

    subplot(211)
    title(node.name)

    balancing = node.balancing
    transmission_import = -node.mismatch*(node.mismatch<0) - node.balancing
    transmission_export = node.mismatch*(node.mismatch>0) - node.curtailment

    RES_local = (node.getwind() + node.getsolar()) - node.curtailment  - transmission_export

    fill_between(t,balancing,color=(.5,0.,0.),lw=0)
    fill_between(t,RES_local + balancing, balancing, color=(.0,.0,.5),lw=0)
    fill_between(t, transmission_import + RES_local + balancing, RES_local + balancing,color=(.0,.5,.0),lw=0)

    pp_balancing = Rectangle((0, 0), 1, 1, facecolor=(.5,0.,0.))
    pp_RES_local = Rectangle((0, 0), 1, 1, facecolor=(.0,.0,.5))
    pp_trans_in = Rectangle((0, 0), 1, 1, facecolor=(.0,.5,.0))

    plot(t,node.load,'k-')	
                
    axis(xmin=0,xmax=30*24,ymin=0,ymax=node.mean*2.)

    pp=(pp_balancing,pp_RES_local,pp_trans_in)
    pptxt = ('Local balancing','Local RES','Import RES')

    leg = legend(pp,pptxt);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    subplot(212)

    plot(t,transmission_import,'r-',label='Import')
    plot(t,transmission_export,'g-',label='Export')	
    plot(t,node.curtailment,'y-',label='Curtailment',lw=2)	
        
    axis(xmin=0,xmax=30*24,ymin=0,ymax=node.mean)

    legend()	

    save_figure('plot_ts.png')

#
# data = plot_generation_summary_vs_year(year=linspace(1990,2050,21),node_id=3,lapse=None)
#
def plot_generation_summary_vs_year(year=linspace(1985,2053,21),node_id=3,lapse=50*24,data=None):
    """ Very complicated function. Needs cleaning."""


    if data==None:
        Nodes=makenodes()
        solvers.options['show_progress']=False
        coop=0
        ###################### Building CVXOPT Matrices ###############################
        ###############################################################################
        P,q,G,h,A = generatemat()
        Nnodes=np.size(np.genfromtxt("incidence.txt",delimiter=','),0)-1
        Nlinks=np.size(np.genfromtxt("incidence.txt",delimiter=','),1)-1
        F=np.zeros((Nlinks,70128))
        # they must also be in cvxopt native "matrix" format (not the same as numpy's!)
        P=matrix(P,tc='d')
        q=matrix(q,tc='d')
        G=matrix(G,tc='d')
        h=matrix(h,tc='d')
        A=matrix(A,tc='d')

        #Calculate flows etc.
        Gamma = get_basepath_gamma(year)
        data = []
        for i in arange(len(year)):
            print 'Year: {0:.0f}'.format(year[i])
            sys.stdout.flush()
        
            #Apply year to nodes (Dummy gamma's)
            Nodes = setgammas(Nodes,Gamma.transpose()[i])
            Nodes = setalphas(Nodes,[1.,1.,1.,1.,.9])

            print 'Gamma: ' + str(Gamma.transpose()[i])
            sys.stdout.flush()
            Nodes,F=runtimeseries(Nodes,F,P,q,G,h,A,coop=0,lapse=lapse)

            add_colored_import(Nodes, F, node_id=None, lapse=lapse)
            #for node in Nodes:
            #    node.colored_import=zeros(node.nhours)

            data.append(deepcopy(Nodes))
    else:
        Nodes = data[0]

    #Calculate averages to be displayed
    balancing_av, import_av, RES_local_av, curtailment_av, export_av = zeros(year.shape), zeros(year.shape), zeros(year.shape), zeros(year.shape), zeros(year.shape)
    colored_import_av = zeros((len(year),len(Nodes)))
    for i in arange(len(year)):
        print data[i][node_id].name, data[i][node_id].gamma, data[i][node_id].balancing.mean()
        sys.stdout.flush()
        
        #Power used locally
        balancing_av[i] = data[i][node_id].getlocalBalancing()[:lapse].mean()/data[i][node_id].mean  #Wrong. Balancing can be exported too!!!!
        RES_local_av[i] = data[i][node_id].getlocalRES()[:lapse].mean()/data[i][node_id].mean
        import_av[i] = data[i][node_id].getimport()[:lapse].mean()/data[i][node_id].mean
        
        #Power not used locally
        export_av[i] = data[i][node_id].getexport()[:lapse].mean()/data[i][node_id].mean
        curtailment_av[i] = data[i][node_id].curtailment[:lapse].mean()/data[i][node_id].mean
        
        #Colored import
        colored_import_av[i] = mean(data[i][node_id].colored_import.transpose()[:lapse], axis=0)/data[i][node_id].mean
        
        
    figure(1); clf()

    title(Nodes[node_id].name)

    plot(year,balancing_av,label='Balancing')
    plot(year,RES_local_av,label='Loacal RES')
    plot(year,import_av,label='Import')

    plot(year,export_av,label='Export')
    plot(year,curtailment_av,label='Curtailment')

    axis(ymin=0, ymax =1.5, xmin=amin(year), xmax=amax(year))
    xlabel('Year')
    ylabel('Power [local av.l.h.]')

    leg = legend()
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    savename = 'plot_generation_summary_vs_year_' + Nodes[node_id].name.replace(' ','_') + '.png'
    save_figure(savename)

    figure(1); clf()

    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])

    ax1 = axes()

    fill_between(year,balancing_av,color=color_balancing,edgecolor='k',lw=.5)
    fill_between(year,RES_local_av+balancing_av,balancing_av,color=color_RES,edgecolor='k',lw=.5)
    fill_between(year,import_av+RES_local_av+balancing_av,RES_local_av+balancing_av,color=color_import,edgecolor='k',lw=.5)
    fill_between(year,export_av+import_av+RES_local_av+balancing_av,import_av+RES_local_av+balancing_av,color=color_export,edgecolor='k',lw=.5)
    fill_between(year,curtailment_av+export_av+import_av+RES_local_av+balancing_av,export_av+import_av+RES_local_av+balancing_av,color=color_curtailment,edgecolor='k',lw=.5)

    axis(ymin=0, ymax =2.05, xmin=amin(year), xmax=amax(year))
    yticks(arange(0,2.05,.5))
    xlabel('Reference year')
    ylabel('Power [av.l.h.]')

    divider = make_axes_locatable(plt.gca())
    ax2 = divider.append_axes("right", "0%", pad="0%")
    #ax2 = axes(ax1.get_position())

    #ax2 = twinx()
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    axis(ymin=0,ymax=2.05)
    ylabel('[GW]')
    xticks([0],[''])
    yticks(ax1.get_yticks(),around(ax1.get_yticks()*Nodes[node_id].mean/1e3,1))

    tight_layout(pad=.3)
    savename = 'plot_generation_summary_vs_year_Stacked_' + Nodes[node_id].name.replace(' ','_') + '.png'
    save_figure(savename)

    figure(1); clf()

    colors = ['k','b','g','r','m']

    fill_between(year,cumsum(colored_import_av,axis=1).transpose()[0],label=Nodes[node_id].name,color=colors[0],lw=0)
    for i in arange(1,len(Nodes)):
        fill_between(year,cumsum(colored_import_av,axis=1).transpose()[i],cumsum(colored_import_av,axis=1).transpose()[i-1],label=Nodes[i].name,color=colors[i],lw=0)

    legend()

    axis(ymin=0, xmin=amin(year), xmax=amax(year))

    savename = 'plot_generation_summary_vs_year_Colored_import_' + Nodes[node_id].name.replace(' ','_') + '.png'
    save_figure(savename)

    return data

def plot_colored_import_export(year, data, colors=colors_countries, lapse=None):

    Nodes = data[0]

    if lapse==None:
        lapse=Nodes[0].mismatch.shape[0]

    for node_id in arange(len(Nodes)):

        #Calculate average import
        colored_import_av = zeros((len(year),len(Nodes)))
        for	i in arange(len(year)):
            colored_import_av[i] = mean(data[i][node_id].colored_import.transpose()[:lapse], axis=0)#/data[i][node_id].mean

        #Calculate export
        colored_export_av = zeros((len(year),len(Nodes)))
        for	i in arange(len(year)):
            for j in arange(len(Nodes)):
                colored_export_av[i][j] = mean(data[i][j].colored_import[node_id])#/data[i][node_id].mean

        #Set plot options	
        matplotlib.rcParams['font.size'] = 10

        figure(1); clf()

        gcf().set_dpi(300)
        gcf().set_size_inches([5.25,3.5])

        ax1 = axes([.11,.565,.885,.42])
        #subplot(211)
        #title(Nodes[node_id].name)

        pp = []; pp_text=[]
        fill_between(year,cumsum(colored_import_av,axis=1).transpose()[0],color=colors[0],edgecolor='k',lw=.5)
        pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[0]))
        pp_text.append(Nodes[0].name)
        for i in arange(1,len(Nodes)):
            fill_between(year,cumsum(colored_import_av,axis=1).transpose()[i],cumsum(colored_import_av,axis=1).transpose()[i-1],label=Nodes[i].name,color=colors[i],edgecolor='k',lw=.5)
            pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[i]))
            pp_text.append(Nodes[i].name)

        leg = legend(tuple(pp),tuple(pp_text),loc='upper left')
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        axis(ymin=0, xmin=1995, xmax=amax(year))
        #xlabel('Reference year')
        ylabel('Power [MW]')

        #subplot(212)
        ax2 = axes([.11,.065,.885,.42])
        pp = []; pp_text=[]
        fill_between(year,cumsum(colored_export_av,axis=1).transpose()[0],color=colors[0],edgecolor='k',lw=.5)
        pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[0]))
        pp_text.append(Nodes[0].name)
        for i in arange(1,len(Nodes)):
            fill_between(year,cumsum(colored_export_av,axis=1).transpose()[i],cumsum(colored_export_av,axis=1).transpose()[i-1],label=Nodes[i].name,color=colors[i],edgecolor='k',lw=.5)
            pp.append(Rectangle((0, 0), 1, 1, facecolor=colors[i]))
            pp_text.append(Nodes[i].name)

        leg = legend(tuple(pp),tuple(pp_text),loc='upper left')
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        axis(ymin=0, xmin=1995, xmax=amax(year))
        #xlabel('Reference year')
        ylabel('Power [MW]')

        savename = 'plot_generation_summary_vs_year_Colored_import_' + Nodes[node_id].name.replace(' ','_') + '.png'
        save_figure(savename)


def get_colored_flow(flow, export, incidence_matrix='incidence.txt'):
	"""flow: vector of flows at time t. export: vector of export at each node at time t."""

	if type(incidence_matrix)==str:
		K = np.genfromtxt(incidence_matrix,delimiter=",",dtype='d') #Incidence matrix.
		K = K[1:].transpose()[1:].transpose() #Remove dummy row/column.
	else:
		K = incidence_matrix
	
	Kf = K*kron(ones((K.shape[0],1)),-flow) #Modified incidence matrix that has positive values for the flow into a node.
	
	FF = array(mat(abs(K))*mat(Kf).transpose()) #Flow matrix with positive values indicating the positive flow into a node.
	FF = get_positive(floor(FF))
	
	#"Column sum" = 1
	for i in arange(FF.shape[1]):
		sum_ = (FF.transpose()[i].sum() + export[i])
		FF.transpose()[i] = FF.transpose()[i]/sum_
		export[i] = export[i]/sum_
	
	#Calculate color matrix	
	try:	
		C = -mat(diag(export))*inv(mat(FF)-mat(eye(FF.shape[0])))	
	except LinAlgError:
		print "Error (dfkln387c): Singular matrix"
		print mat(FF)-mat(eye(FF.shape[0]))

		C = zeros(FF.shape)
	
	return array(FF), array(C)

def add_colored_import(Nodes, F, node_id=None, incidence_matrix='incidence.txt', lapse=None):
	
	if lapse==None:
		lapse=Nodes[0].mismatch.shape[0]
	
	if type(incidence_matrix)==str:
		K = np.genfromtxt(incidence_matrix,delimiter=",",dtype='d') #Incidence matrix.
		K = K[1:].transpose()[1:].transpose() #Remove dummy row/column.
	else:
		K = incidence_matrix
			
	for t in arange(lapse):
	
		export_ = array([Nodes[i].getexport()[t] for i in arange(len(Nodes))])
		import_ = array([Nodes[i].getimport()[t] for i in arange(len(Nodes))])
	
		FF, C = get_colored_flow(F.transpose()[t], copy(export_), incidence_matrix=K)
	
		CC = C*kron(ones((K.shape[0],1)),import_)
	
		#Update Node(s)
		if node_id == None:
			for node_id_ in arange(len(Nodes)):
				Nodes[node_id_].set_colored_import_i(t,CC.transpose()[node_id_])
		else:
			Nodes[node_id].set_colored_import_i(t,CC.transpose()[node_id])


### Local utilities
def get_basepath_gamma(year,filename='basepath_gamma.npy'):

    print "Loading: {0}. Warning columns not pre-labeled!!".format(filename)
    data = np.load(filename)

    Gamma = zeros((len(data)-1,len(year)))
    for i in arange(len(Gamma)):
        Gamma[i] = interp(year,data[0],data[i+1])

    return Gamma

def generate_basepath_gamma_alpha(txtfile='basepath_wind_solar.csv',year0=1980,year_hist=2009,plot_on=False):

    print "Loading: {0}. Warning columns not pre-labeled!!".format(txtfile)
    txttitles = ['Year','Norway (wind)','Norway (solar)','Sweden (wind)','Sweden (solar)','Denmark West (wind)','Denmark West (solar)','Denmark East (wind)','Denmark East (solar)','Germany North (wind)','Germany North (solar)']
    txtlabels = ['year','N-wind','N-solar','S-wind','S-solar','DKW-wind','DKW-solar','DKE-wind','DKE-solar','DE-wind','DE-solar']
    data = np.genfromtxt(txtfile,delimiter=',',skip_header=1)
    
    year = array(data.transpose()[0])
    
    if plot_on==True:
        p_historical = year<=year_hist
    else:
        p_historical = None
    
    gamma, alpha_w = [], []
    for i in arange(1,data.shape[1]-1,2):
        wind = data.transpose()[i]
        solar = data.transpose()[i+1] 

        if ~all(isnan(wind)):
            i_data = find(~isnan(wind))
            gamma_wind = get_logistic_fit(year[i_data],wind[i_data],year0=year0,year=year,plot_on=plot_on,p_historical=p_historical[i_data],txtlabel=txtlabels[i],txttitle=txttitles[i])
        else:
            gamma_wind = zeros(wind.shape)
            
        if ~all(isnan(solar)):
            i_data = find(~isnan(solar))
            gamma_solar = get_logistic_fit(year[i_data],solar[i_data],year=year,plot_on=plot_on,p_historical=p_historical[i_data],txtlabel=txtlabels[i+1],txttitle=txttitles[i+1])
        else:
            gamma_solar = zeros(wind.shape)

        gamma.append(gamma_wind+gamma_solar)
        alpha_w.append(gamma_wind/(gamma_wind+gamma_solar))

    if plot_on==True:
        weight = array([13672.325571167074,16629.057925672601,2312.6684067162068,1630.6221299198942,18730.927464722623])
        plot_basepath_gamma_alpha(year,array(gamma),array(alpha_w),weight)

    np.save('basepath_gamma',concatenate([array(year,ndmin=2),array(gamma)]))
    print 'Saved file: basepath_gamma.npy'
    np.save('basepath_alpha_w',concatenate([array(year,ndmin=2),array(alpha_w)]))
    print 'Saved file: basepath_alpha_w.npy'


def plot_basepath_gamma_alpha(year,gamma,alpha_w,weight,txtlabels=None):

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    close(1); figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])    
    
    pp = []
    for i in arange(len(gamma)):
        pp_ = plot(year,gamma[i],'-',lw=1.5,color=colors_countries[i])
        pp.extend(pp_)
    
    gamma_mean = sum(gamma*kron(array(weight,ndmin=2).transpose(),ones(gamma.shape[1]))/sum(weight),axis=0)
    gamma_mean_DK = sum(gamma[2:4]*kron(array(weight[2:4],ndmin=2).transpose(),ones(gamma[2:4].shape[1]))/sum(weight[2:4]),axis=0)
    
    pp_mean = plot(year,gamma_mean,'k-',lw=2)
    pp.extend(pp_mean)
    
    pp_mean_DK = plot(year,gamma_mean_DK,'k--',lw=2)
    pp.extend(pp_mean_DK)
    
    
    axis(xmin=amin(year),xmax=2053,ymin=0,ymax=1.3)
    xlabel('Reference year')
    ylabel(r'Gross share of electricity demand ($\gamma_X$)')
    
    pp_text = ['Norway','Sweden','Denmark West','Denmark East','Germany North','Region (mean)','Denmark (mean)']
    
    leg = legend(pp,pp_text,loc='upper left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    tight_layout(pad=.2)
    save_figure('plot_basepath_gamma_vs_year.pdf')
    
    close(1); figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([6.5,4.3])    
    
    pp = []
    for i in arange(len(gamma)):
        pp_ = plot(gamma_mean_DK,gamma[i],'-',lw=1.5,color=colors_countries[i])
        pp.extend(pp_)
    
    pp_mean = plot(gamma_mean_DK,gamma_mean,'k-',lw=2)
    pp.extend(pp_mean)
    
    pp_mean_DK = plot(gamma_mean_DK,gamma_mean_DK,'k--',lw=2)
    pp.extend(pp_mean_DK)
    
    
    axis(xmin=0,xmax=1.025,ymin=0,ymax=1.3)
    xlabel(r'Danish gross share of electricity demand ($\gamma_{DK}$)')
    ylabel(r'Gross share of electricity demand ($\gamma_X$)')
    
    pp_text = ['Norway','Sweden','Denmark West','Denmark East','Germany North','Region (mean)','Denmark (mean)']
    
    leg = legend(pp,pp_text,loc='upper left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    tight_layout(pad=.2)
    save_figure('plot_basepath_gamma_vs_gamma_DK.pdf')

def get_logistic_fit(p_year,p_gamma,year0=1980,year=None,plot_on=False,p_historical=None,txtlabel=None,txttitle=None):
    
    if year==None:
        year = arange(amin(p_year),amax(p_year),1)
    
    fitfunc = lambda p, x: p[0]*abs(p[3])*exp(p[1]*(x-year0))/(abs(p[3])+p[0]*(exp(p[1]*(x-year0))-1.))
    errfunc = lambda p, x, y, weight: (fitfunc(p, x) - y)/weight # Distance to the target function

    # p = [P_0, r, year0, K]
    p_0 = [amin(p_gamma), .01, year0, amax(p_gamma)] # Initial guess for the parameters
    p_weight = ones(p_year.shape)
    p_fit, success = optimize.leastsq(errfunc, p_0[:], args=(p_year, p_gamma, p_weight))

    if plot_on==True:
        plot_logistic_fit(year,fitfunc(p_fit,year),p_year,p_gamma,p_historical,txtlabel,txttitle)

    return fitfunc(p_fit,year)


def plot_logistic_fit(year,gamma_fit,p_year,p_gamma,p_historical=None,txtlabel=None,txttitle=None):

    if p_historical==None:
        p_historical = ones(year.shape)

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

##
# Uncomment the line below the function. Need a better method for storing/retriving data.
#
def load_pickled_data(filename='data_200111122.p',path='./data/'):

    data=pickle.load(open(path+filename,'rb'))

    return data

# data = load_pickled_data('data_200111122.p')

### To be file utilities

def get_positive(x):
	
	return x*(x>0.)  #Possibly it has to be x>1e-10.    
    
def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()
    
    
    
    
    
    

