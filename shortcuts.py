#
#  shortcuts.py
#  General utilities that I like to use.
#
#  Created by Gorm Bruun Andresen on 24/11/2011.
#  Copyright (c) 2011 Department of Engineering, University of Aarhus. All rights reserved.
#

#Standard modules
from pylab import *
from scipy import *
import os

#Special functions
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_positive(x):
    """Replaces all negative values with zeros."""
    return x*(x>0.)  #Possibly it has to be x>1e-10.    
    
def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
    """Wraper for savefig(). Saves figure to path+filename and prints a meassage to the screen."""
    
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()

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