
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""

article = False

import os

import numpy

# Plot customization
import matplotlib

# Use LaTeX for all text
matplotlib.rcParams['text.usetex'] = True

# Markers and line widths
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['lines.markersize'] = 8

# Font Sizes
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16

# DPI of output images
matplotlib.rcParams['savefig.dpi'] = 100

import matplotlib.pyplot as plt
import datetime

from clawpack.visclaw import colormaps
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata

import clawpack.geoclaw.surge.plot as surge

try:
    from setplotfg import setplotfg
except:
    setplotfg = None

def setplot(plotdata):
    r"""Setplot function for surge plotting"""
    

    plotdata.clearfigures()  # clear any old figures,axes,items data
    # plotdata.format = 'binary'

    fig_num_counter = surge.figure_counter()

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir,'claw.data'))
    amrdata = amrclaw.AmrclawInputData(clawdata)
    amrdata.read(os.path.join(plotdata.outdir,'amr.data'))
    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir,'geoclaw.data'))
    surge_data = geodata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir,'surge.data'))

    surge.pressure_field = surge_data.pressure_index - 1

    surface_limits = [-1.0, 1.0]
    speed_max = 1.25
    speed_limits = [0.0, speed_max]
    pressure_limits = [surge_data.ambient_pressure / 100.0, 
                       surge_data.ambient_pressure * 1.5 / 100.0]

    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    #
    #  Surface
    #
    plotfigure = plotdata.new_plotfigure(name='Surface - Entire Domain', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = [clawdata.lower[0], clawdata.upper[0]]
    plotaxes.ylimits = [clawdata.lower[1],clawdata.upper[1]]

    surge.add_surface_elevation(plotaxes, plot_type="pcolor", bounds=surface_limits)
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0]

    #
    #  Water Speed
    #
    plotfigure = plotdata.new_plotfigure(name='Currents - Entire Domain',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = [clawdata.lower[0], clawdata.upper[0]]
    plotaxes.ylimits = [clawdata.lower[1],clawdata.upper[1]]

    surge.add_speed(plotaxes, plot_type="pcolor", bounds=speed_limits)
    plotaxes.plotitem_dict['speed'].amr_patchedges_show = [0]

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True
    # plotfigure.kwargs['figsize'] = (16,10)

    def gauge_after_axes(cd):

        axes = plt.gca()

        # Fix up plot
        # axes.set_title('Station %s' % cd.gaugeno)
        # axes.set_xlabel('Days relative to landfall')
        # axes.set_ylabel('Surface (m)')
        # axes.set_xlim([-2,1])
        # axes.set_ylim([-1,5])
        # axes.set_xticks([-2,-1,0,1])
        # axes.set_xticklabels([r"$-2$",r"$-1$",r"$0$",r"$1$"])
        axes.grid(True)
        # axes.legend()

        plt.hold(False)


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = [-2,1]
    # plotaxes.xlabel = "Days from landfall"
    # plotaxes.ylabel = "Surface (m)"
    # plotaxes.ylimits = [-1,5]
    plotaxes.title = 'Surface'
    plotaxes.afteraxes = gauge_after_axes

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    # ================
    #  Pressure Field
    # ================
    plotfigure = plotdata.new_plotfigure(name="Pressure", 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Pressure"
    plotaxes.scaled = True
    plotaxes.xlimits = [clawdata.lower[0], clawdata.upper[0]]
    plotaxes.ylimits = [clawdata.lower[1], clawdata.upper[1]]

    surge.add_pressure(plotaxes, plot_type="pcolor", bounds=pressure_limits)
    plotaxes.plotitem_dict['pressure'].amr_patchedges_show = [0]
    plotaxes.plotitem_dict['pressure'].amr_celledges_show = [0]

    # ==============================================================
    #  Debugging Plots, only really work if using interactive plots
    # ==============================================================
    #
    # Water Velocity Components
    #
    plotfigure = plotdata.new_plotfigure(name='Velocity Components - Entire Domain',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = False

    # X-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(121)"
    plotaxes.title = 'Velocity, X-Component'
    plotaxes.scaled = True

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surge.water_u
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1]
    plotitem.pcolor_cmax = speed_limits[1]
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]

    surge.add_land(plotaxes)

    # Y-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(122)"
    plotaxes.title = 'Velocity, Y-Component'
    plotaxes.scaled = True

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surge.water_v
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1]
    plotitem.pcolor_cmax = speed_limits[1]
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    
    surge.add_land(plotaxes)

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'            # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

