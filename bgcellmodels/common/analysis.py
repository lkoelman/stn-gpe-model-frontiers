"""
Common plotting/recording functions for cell experiments.

@author Lucas Koelman
@date 26/10/2016
"""

import neuron
h = neuron.h
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import collections
import re

import logging
logger = logging.getLogger('analysis')

################################################################################
# Recording & Analysis functions
################################################################################

def rec_currents_activations(traceSpecs, sec_tag, sec_loc, ion_species=None, 
                             currents=True, chan_states=True, ion_conc=True):
    """
    Specify trace recordings for all ionic currents
    and activation variables 

    @param traceSpecs   collections.OrderedDict of trace specifications

    @param  sec_tag     tag for the section in the section dictionary passed
                        to analysis.recordTraces. Note that the first two
                        characters must be unique

    @param  sec_loc     location in the section to record the variable
                        (location maps to a segment where the var is recorded)

    @effect             for each ionic current, insert a trace specification
                        for the current, open fractions, activation, 
                        and inactivation variables


    EXAMPLE USAGE

        secs = {'soma': soma, 'dend': dendsec}

        traceSpecs = collections.OrderedDict()
        traceSpecs['V_soma'] = {'sec':'soma', 'loc':0.5, 'var':'v'}
        traceSpecs['t_global'] = {'var':'t'}

        rec_currents_activations(traceSpecs, 'soma', 0.5)
        rec_currents_activations(traceSpecs, 'dend', 0.9, ion_species=['ca','k'])

        recordStep = 0.025
        recData = analysis.recordTraces(secs, traceSpecs, recordStep)
        h.init()
        h.run()

        figs_soma, cursors_soma = plot_currents_activations(recData, recordStep, sec_tag='soma')
        figs_dend, cursors_dend = plot_currents_activations(recData, recordStep, sec_tag='dend')

    """
    if ion_species is None:
        ion_species = ['na', 'k', 'ca', 'nonspecific']

    # Derive suffix for traces from section tag
    suf = '_' + sec_tag[0:2]
    ts = traceSpecs

    if 'na' in ion_species:
        if currents:
            # Na currents
            ts['I_NaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Na','var':'ina'} # transient sodium
            ts['I_NaP'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'NaL','var':'inaL'} # persistent sodium
        if chan_states:
            # Na Channel open fractions
            ts['O_NaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Na','var':'o'}
            # NOTE: leak currents such as I_NaL and gpas_STh are always open
            # Na channel activated fractions
            ts['A_NaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Na','var':'m'}
            # Na channel inactivated fractions
            ts['B_NaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Na','var':'h'}

    if 'k' in ion_species:
        if currents:
            # K currents
            ts['I_KDR'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'KDR','var':'ik'}
            ts['I_Kv3'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Kv31','var':'ik'}
            ts['I_KCa'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'sKCa','var':'isKCa'}
        if chan_states:
            # K channel open fractions
            ts['O_KDR'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'KDR','var':'n'}
            ts['O_Kv3'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Kv31','var':'p'}
            ts['O_KCa'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'sKCa','var':'w'}
            # K channels activated fractions - same as open fractions (single state variable)
            ts['A_KDR'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'KDR','var':'n'}
            ts['A_Kv3'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Kv31','var':'p'}
            ts['A_KCa'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'sKCa','var':'w'}
            # K channels inactivated fractions - not present (single state variable)

    if 'ca' in ion_species:
        # Ions
        ts['C_CaL_cai'+suf] = {'sec':sec_tag,'loc':sec_loc,'var':'cai'} # intracellular calcium concentration
        ts['C_CaT_cai'+suf] = {'sec':sec_tag,'loc':sec_loc,'var':'cai'} # intracellular calcium concentration
        if currents:
            # Ca currents
            ts['I_CaL'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'iLCa'}
            ts['I_CaN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'iNCa'}
            ts['I_CaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'CaT','var':'iCaT'}
        if chan_states:
            # Ca channel open fractions
            ts['O_CaL'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'o_L'}
            ts['O_CaN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'o_N'}
            ts['O_CaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'CaT','var':'o'}
            # Ca channels activated fractions
            ts['A_CaL'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'q'} # shared activation var for L/N
            ts['A_CaN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'q'} # shared activation var for L/N
            ts['A_CaT'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'CaT','var':'r'}
            # Ca channels inactivated fractions
            ts['B_CaL'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'h'} # Ca-depentdent inactivation
            ts['B_CaN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'HVA','var':'u'} # V-dependent inactivation
            ts['B_CaTf'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'CaT','var':'s'} # fast inactivation
            ts['B_CaTs'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'CaT','var':'d'} # slow inactivation

    if 'nonspecific' in ion_species:
        if currents:
            # Nonspecific currents
            ts['I_HCN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Ih','var':'ih'}
        if chan_states:
            # Nonspecific channel open fractions
            ts['O_HCN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Ih','var':'f'}
            # Nonspecific channels activated fractions - same as open fractions (single state variable)
            ts['A_HCN'+suf] = {'sec':sec_tag,'loc':sec_loc,'mech':'Ih','var':'f'}
            # Nonspecific channels activated fractions - not present (single state variable)

def plot_currents_activations(recData, recordStep, timeRange=None, sec_tag=None):
    """
    Plot currents and (in)activation variable for each ionic
    current in the same axis. Ionic currents are grouped per ion
    in one figure, and the x-axes are synchronized for zooming
    and panning.

    @param recData      traces recorded using a trace specification provided
                        by rec_currents_activations()

    @return             tuple figs, cursors where figs is a list of figures
                        that were created and cursors a list of cursors

    EXAMPLE USAGE:      see function rec_currents_activations()
    """
    figs = []
    cursors = []

    # Plot activations-currents on same axis per current
    ions_chans = [('NaT', 'NaP', 'HCN'), ('KDR', 'Kv3', 'KCa'), ('CaL', 'CaN', 'CaT')]
    for iontraces in ions_chans: # one figure for each ion
        fig, axrows = plt.subplots(len(iontraces), 1, sharex=True) # one plot for each channel
        for i, trace_abbr in enumerate(iontraces):
            # Find traces that are marked with the tracee abbreviation (e.g. 'CaT')
            if sec_tag is None:
                pat = re.compile(r'^[A-Z]_' +trace_abbr) # find 'char+_+abbr' at beginning of word
            else:
                sec_suffix = sec_tag[0:2]
                pat = re.compile(r'^[A-Z]_' + trace_abbr + r'\w+' + sec_suffix + r'$')
            chanFilter = lambda x: re.search(pat, x) # variables plotted on left axis match this filter
            twinFilter = lambda x: x.startswith('I_') # vars plotted on right axis match this filter
            # Plot traces that match pattern
            cumulPlotTraces(recData, recordStep, showFig=False, 
                                fig=None, ax1=axrows[i], yRange=(-0.1,1.1), timeRange=timeRange,
                                includeFilter=chanFilter, twinFilter=twinFilter)
            # Add figure interaction
            cursor = matplotlib.widgets.MultiCursor(fig.canvas, fig.axes, 
                        color='r', lw=1, horizOn=False, vertOn=True)
        figs.append(fig)
        cursors.append(cursor)

    return figs, cursors


def recordTraces(secs, traceSpecs, recordStep=0.05, duration=None, recData=None):
    """
    Record the given traces from section

    @param      (optional) recData : Collections.OrderedDict()
                Existing dicionary containing recording vectors

    @returns    recData : Collections.OrderedDict(str -> Hoc.Vector)
                Ordered dictionary containing recorded trace data after running simulation

    @see        Based on NetPyne's Cell.RecordTraces in file:
                https://github.com/Neurosim-lab/netpyne/blob/master/netpyne/cell.py

    USAGE
    -----

        secs = {'soma': soma, 'dend': dends[1], 'izhpp': izh, 'synpp': syn}

        traceSpecs = {
            'V_izhi': {'pointp': 'izh'},
            'g_syn': {'pointp': 'synpp', 'var': 'g'}
            'V_soma':{'sec':'soma','loc':0.5,'var':'v'},
            'GP_RT_cai':{'sec':'soma','loc':0.5,'var':'cai'},
            'GP_RT_ainf':{'sec':'soma','loc':0.5,'mech':'gpRT','var':'a_inf'}, 
            'GP_RT_r':{'sec':'soma','loc':0.5,'mech':'gpRT','var':'r'},
            'STN_r':{'sec':'soma','loc':0.5,'mech':'stn','var':'r'},
            'STN_p':{'sec':'soma','loc':0.5,'mech':'stn','var':'p'},
            'STN_q':{'sec':'soma','loc':0.5,'mech':'stn','var':'q'},
        }

    WARNING: For multithreaded execution, section _must_ have POINT_PROCESS
             associated with it to identify the tread. Hence, a PointProcessMark
             (see ppmark.mod) will be inserted if no PP present.
    """
    if recData is None:
        recData = collections.OrderedDict() # empty dict for storing recording vectors
    pp_markers = []

    for trace, spec in traceSpecs.iteritems():
        if trace in recData:
            logger.warning("Trace named {} already exists in data dictionary. Overwriting.".format(trace))

        ptr = None
        pp = None
        if 'loc' in spec or 'seg' in spec:

            if 'seg' in spec:
                rec_hobj = secs[spec['seg']]
            else:
                rec_hobj = secs[spec['sec']]

            if isinstance(rec_hobj, neuron.nrn.Segment):
                seg = rec_hobj
                sec = seg.sec
                if 'loc' in spec:
                    seg = sec(spec['loc'])
            else:
                sec = rec_hobj
                seg = sec(spec['loc'])
            
            # Get pointer/reference to variable to record
            if 'mech' in spec:  # eg. soma(0.5).hh._ref_gna
                ptr = seg.__getattribute__(spec['mech']).__getattribute__('_ref_'+spec['var'])
            else:
                # No mechanism. E.g. soma(0.5)._ref_v
                ptr = seg.__getattribute__('_ref_'+spec['var'])
            
            # find a POINT_PROCESS in segment to improve efficiency
            # seg_pps = seg.point_processes()
            # if len(seg_pps) > 0:
            #     pp = seg_pps[0]
            # else:
            pp = h.PointProcessMark(seg)
            pp_markers.append(pp)
        
        elif 'pointp' in spec: # POINT_PROCESS objects
            if spec['pointp'] in secs:
                pp = secs[spec['pointp']]
                ptr = pp.__getattribute__('_ref_'+spec['var'])
        
        elif 'netcon' in spec: # NetCon and objects that send events
            # Dont make ptr
            nc = secs[spec['netcon']]
            vec = h.Vector()
            nc.record(vec)
            recData[trace] = vec

        else: # global vars, e.g. h.t
            ptr = h.__getattribute__('_ref_'+spec['var'])

        if ptr:  # if pointer has been created, then setup recording
            if duration is not None:
                recData[trace] = h.Vector(duration/recordStep+1).resize(0)
            else:
                recData[trace] = h.Vector()

            if pp is not None:
                recData[trace].record(pp, ptr, recordStep)
            else:
                recData[trace].record(ptr, recordStep)

    return recData, pp_markers


def to_numpy(vec):
    """
    Convert vector/array from any type to numpy array.
    """
    if isinstance(vec, neuron.hoc.HocObject):
        return vec.as_numpy()
    elif isinstance(vec, np.ndarray):
        return vec
    else:
        return np.array(vec)


def plotTraces(traceData, recordStep, timeRange=None, oneFigPer='cell', 
                includeTraces=None, excludeTraces=None, labelTime=False,
                showFig=True, colorList=None, lineList=None, yRange=None,
                traceSharex=False, showGrid=True, title=None, traceXforms=None,
                fontSize=10, labelRotation=-90, singleAxis=False, **plot_kwargs):
    """
    Plot previously recorded traces

    - traceData
        dict(trace_name -> h.Vector()) containing recorded traces

    - traceXforms
        dict(trace_name -> function) containg transformation to apply
        to trace before plotting

    - timeRange ([start:stop])
        Time range of spikes shown; if None shows all (default: None)
    
    - oneFigPer ('cell'|'trace')
        Whether to plot one figure per cell (showing multiple traces) 
        or per trace (showing multiple cells) (default: 'cell')
    
    - showFig (True|False)
        Whether to show the figure or not (default: True)
    
    - labelTime (True|False)
        whether to show the time axis label (default:True)
    
    - includeTraces (['V_soma', ...])
        traces to include in this plot
    
    - excludeTraces (['V_soma', ...])
        traces to exclude in this plot
    
    - yRange
        y limit for range, e.g. (0, 60). Can be a tuple, list or dict with traces as keys
    
    - traceSharex
        if True, all x-axes will be shared (maintained during zooming/panning), else
        if an axes object is provided, share x axis with that axes

    - **plot_kwargs : **dict
        Extra parameters for calls to plt.plot() and plt.figure()
        
    Returns figure handles
    """
    trace_vecs = traceData
    tracesList = traceData.keys()
    if includeTraces is not None:
        tracesList = [trace for trace in tracesList if trace in includeTraces]
    if excludeTraces is not None:
        tracesList = [trace for trace in tracesList if trace not in excludeTraces]

    # Convert to numpy
    traceData = collections.OrderedDict()
    for trace_name in tracesList:
        traceData[trace_name] = to_numpy(trace_vecs[trace_name])

    # time range
    first_trace = traceData[tracesList[0]]
    # if recordStep is None:
    #   recordStep = first_trace[1] - first_trace[0]
    if timeRange is None:
        timeRange = [0, first_trace.size*recordStep]

    if colorList is None:
        colorList = [[0.42,0.67,0.84], [0.90,0.76,0.00], [0.42,0.83,0.59], [0.90,0.32,0.00],
                [0.34,0.67,0.67], [0.90,0.59,0.00], [0.42,0.82,0.83], [1.00,0.85,0.00],
                [0.33,0.67,0.47], [1.00,0.38,0.60], [0.57,0.67,0.33], [0.5,0.2,0.0],
                [0.71,0.82,0.41], [0.0,0.2,0.5]]

    # Sharing axes for comparing signals
    shared_ax = None
    if isinstance(traceSharex, matplotlib.axes.Axes):
        shared_ax = traceSharex

    figs = []
    traces_axes = {t: None for t in tracesList}
    label_fontsize = fontSize
    fig_size = plot_kwargs.pop('figsize', None)
    plot_kwargs.setdefault('linewidth', 1.0)
    for itrace, trace in enumerate(tracesList):

        found_ax = False # have not found axis for this trace
        
        if oneFigPer == 'cell':
            
            if itrace == 0:

                figs.append(plt.figure(figsize=fig_size))
                ax = plt.subplot(len(tracesList), 1, itrace+1, sharex=shared_ax)
                if traceSharex:
                    shared_ax = ax
                found_ax = True    

            elif singleAxis:
                
                if isinstance(singleAxis, (list, tuple)):
                    trace_grp = next((g for g in singleAxis if trace in g), None)
                    if trace_grp is None:
                        found_ax = False
                    else:
                        ax = next((traces_axes[t] for t in trace_grp if traces_axes[t] is not None), None)
                        found_ax = ax is not None

                else:
                    pass # ax = ax

            if not found_ax:
                ax = plt.subplot(len(tracesList), 1, itrace+1, sharex=shared_ax)
        
        else: # one separate figure per trace
            
            fig = plt.figure(figsize=fig_size)
            figs.append(fig)
            
            ax = plt.subplot(111, sharex=shared_ax)
            if itrace==0 and traceSharex:
                shared_ax = ax

        traces_axes[trace] = ax

        # Get data to plot
        if (traceXforms is not None) and (trace in traceXforms):
            xform = traceXforms[trace]
            tracevec = xform(traceData[trace])
        else:
            tracevec = traceData[trace]
        
        data = tracevec[int(timeRange[0]/recordStep):int(timeRange[1]/recordStep)]
        t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep)

        plt.plot(t[:len(data)], data, color=colorList[itrace%len(colorList)],
                 **plot_kwargs)

        # Axes ranges/labels
        # plt.ylabel(trace, fontsize=label_fontsize)
        ax.set_ylabel(trace, rotation=labelRotation, fontsize=label_fontsize)
        
        if labelTime:
            plt.xlabel('Time (ms)', fontsize=label_fontsize)
        
        if isinstance(yRange, dict):
            if trace in yRange:
                plt.ylim(yRange[trace])
        
        elif yRange is not None:
            plt.ylim(yRange)
        
        plt.xlim(timeRange)

        # Customize the grid
        ax.grid(showGrid)

    if title:
        plt.suptitle(title) # suptitle() is Fig title, title() is ax title

    if showFig:
        plt.show(block=False)
    
    return figs

# Define visually appealing colors
allcolors = [
    '#7742f4', # Dark purple
    [0.90,0.76,0.00], # Ochre
    [0.42,0.83,0.59], # soft pastel green
    [0.90,0.32,0.00], # pastel red brick
    [0.90,0.59,0.00], # OrangeBrown
    '#f442c5', # Pink
    '#c2f442', # Lime
    [1.00,0.85,0.00], # hard yellow
    [0.33,0.67,0.47], # dark pastel green
    [1.00,0.38,0.60], [0.57,0.67,0.33], [0.5,0.2,0.0],
    [0.71,0.82,0.41], [0.0,0.2,0.5],
]
greenish = [
    [0.42,0.83,0.59], # soft pastel green
    '#c2f442', # Lime
    [0.33,0.67,0.47], # dark pastel green
]
redish = [
    [0.90,0.32,0.00], # pastel red brick
    '#f442c5', # Bright pink
    [0.90,0.59,0.00], # OrangeBrown
]
blueish = [
    '#7742f4', # Dark purple
    '#c2f442', # Soft cyan blue
    '#0066FF', # pastel blue
]
solid_styles = ['-']
broken_styles = ['--', '-.', ':']


def pick_line_options(hue, solidity, index):
    """
    Pick a line color and solidity specifier based on description
    of the class and an index that wraps around.

    @return     tuple(color, linestyle) that can be passed to the keyword
                arguments 'color' and 'linestyl' of matplotlib.plot
    """
    if hue in ('green', 'greenish'):
        hue_clist = greenish
    elif hue in ('red', 'redish'):
        hue_clist = redish
    elif hue in ('blue', 'blueish'):
        hue_clist = blueish
    elif hue in ('any', 'all'):
        hue_clist = allcolors
    else:
        raise ValueError(hue)
    if solidity == 'solid':
        sol_list = solid_styles
    elif solidity == 'broken':
        sol_list = broken_styles
    else:
        raise ValueError(solidity)
    return hue_clist[index % len(hue_clist)], sol_list[index % len(sol_list)]


def pick_line(trace_name, trace_index, solid_only=False):
    """
    Pick a line style and color based on the trace name
    """
    style_map = {
        'I': (allcolors, solid_styles),
        'V': (allcolors, solid_styles),
        'C': (allcolors, solid_styles),
        'A': (greenish, broken_styles),
        'B': (redish, broken_styles),
        'O': (blueish, broken_styles),
    }
    if solid_only:
        style_map['A'] = (greenish, solid_styles)
        style_map['B'] = (redish, solid_styles)
        style_map['O'] = (blueish, solid_styles)
    default_style = (allcolors, broken_styles)

    # pick a line style
    match_prefix = re.search(r'^[a-zA-Z]', trace_name) # first letter
    if match_prefix:
        prefix = match_prefix.group()
        colors, styles = style_map.get(prefix, default_style)
    else:
        colors, styles = default_style
    return colors[trace_index%len(colors)], styles[trace_index%len(styles)]


def match_traces(recData, matchfun, orderfun=None, reverse=False, pop=False):
    """
    Get ordered dictionary with matching traces.

    @param  matchfun: lambda string -> bool
            Filter function that matches trace names

    @param  pop : bool
            Pop (remove) matched traces from recData
    """

    traces = ((name,data) for name,data in recData.iteritems() if matchfun(name))

    if orderfun is not None:
        traces = sorted(traces, key=orderfun, reverse=reverse) # if no orderfun: keep recData order

    if pop:
        for trace_name, _ in traces:
            recData.pop(trace_name)
    
    return collections.OrderedDict(traces)


def cumulPlotTraces(traceData, recordStep, timeRange=None, cumulate=False,
                    includeTraces=None, excludeTraces=None,
                    showFig=True, colorList=None, lineList=None, 
                    traceSharex=None, fig=None, showGrid=True,
                    yRange=None, twinFilter=None, includeFilter=None,
                    yRangeR=None, ax1=None, solid_only=False, **plot_kwargs):
    """
    Cumulative plot of traces

    @param      **plot_kwargs : **dict
                Extra parameters for calls to plt.plot() and plt.figure()

    @param      fig              if provided, add axs as subplot to this fig

    @param      traceSharex      a matplotlib.axes.Axes object to share the
                            x-axis with (x-axis locked while zooming/panning)

    @param      twinFilter       function that takes as argument a trace name, 
                            returns True if trace should be plotted on right
                            axis (using ax.twinx()), False if on left
    """

    tracesList = traceData.keys()
    if includeFilter is not None:
        tracesList = [trace for trace in tracesList if includeFilter(trace)]
    if includeTraces is not None:
        tracesList = [trace for trace in tracesList if trace in includeTraces]
    if excludeTraces is not None:
        tracesList = [trace for trace in tracesList if trace not in excludeTraces]
    if len(tracesList) == 0:
        print('WARNING: No traces left for plotting after applying filter! Empty figure will be returned.')
        return fig

    if timeRange is None:
        timeRange = [0, traceData[tracesList[0]].size()*recordStep]

    fontsiz=12

    # Get the axes to draw on
    if not fig and not ax1: # create fig and axis
        fig = plt.figure(figsize=plot_kwargs.pop('figsize', None))
        ax1 = fig.add_subplot(111, sharex=traceSharex)
    elif not fig and ax1: # get fig from given axis
        fig = ax1.figure
    elif fig and not ax1: # add new axis to figure
        nax = len(fig.axes)
        for i, ax in enumerate(fig.axes):
            ax.change_geometry(nax+1, 1, i+1) # change grid and position in grid
        ax1 = fig.add_subplot(nax+1, 1, nax+1, sharex=traceSharex)
    if twinFilter is not None:
        ax2 = ax1.twinx()
    else:
        ax2 = None
        twinFilter = lambda x: False

    t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep)
    cumulTrace = np.zeros(int((timeRange[1]-timeRange[0])/recordStep))

    # plot each trace
    lines = []
    for itrace, tracename in enumerate(tracesList):
        tracevec = traceData[tracename].as_numpy()
        data = tracevec[int(timeRange[0]/recordStep):int(timeRange[1]/recordStep)]
        if twinFilter(tracename):
            pax = ax2
        else:
            pax = ax1
        line_color, line_style = pick_line(tracename, itrace, solid_only=solid_only)
        li, = pax.plot(t[:len(data)], data+cumulTrace, label=tracename, 
                        color=line_color, linestyle=line_style)
        lines.append(li)
        if cumulate: cumulTrace += data

    # Axes ranges/labels
    plt.xlim(timeRange)
    if yRange is not None:
        ax1.set_ylim(yRange)
    if yRangeR is not None:
        ax2.set_ylim(yRangeR)

    # Set labels
    L_prefixes = []
    R_prefixes = []
    for tracename in tracesList:
        match = re.search(r'^._', tracename)
        if match and twinFilter(tracename):
            R_prefixes.append(match.group()[:-1])
        elif match:
            L_prefixes.append(match.group()[:-1])
    ax1.set_ylabel(', '.join(set(L_prefixes)), fontsize=fontsiz)
    if ax2:
        ax2.set_ylabel(', '.join(set(R_prefixes)), fontsize=fontsiz)

    # legend for lines plotted in all axes
    ax1.legend(lines, [l.get_label() for l in lines])

    # Customize the grid
    ax1.grid(showGrid)

    if showFig:
        plt.show(block=False)
    return fig


def plotRaster(spikeData, timeRange=None, showFig=True, 
                includeTraces=None, excludeTraces=None, 
                showGrid=True, title=None, **plot_kwargs):
    """
    Plot rastergram of spike times.

    ARGUMENTS:

    - spikeData
        dict<str, array> : {trace_name: spike_times}
        array can be Hoc.Vector(), numpy.array, or any iterable

    - timeRange ([start:stop])
        Time range of spikes shown; if None shows all (default: None)
    
    - showFig (True|False)
        Whether to show the figure or not (default: True)

    - marker
        valid matplotlib marker, e.g. ',' '|'
        see https://matplotlib.org/api/markers_api.html

    @return         tuple(fig, ax)
    """
    plot_kwargs.setdefault('marker', '|')   # clearest for raster plot
    plot_kwargs.setdefault('linestyle', '') # no line through markers
    plot_kwargs.setdefault('snap', True)    # crisper lines
    plot_kwargs.setdefault('color', 'b')

    # Select traces to be plotted
    traceNames = spikeData.keys()
    if includeTraces is not None:
        traceNames = [trace for trace in traceNames if trace in includeTraces]
    if excludeTraces is not None:
        traceNames = [trace for trace in traceNames if trace not in excludeTraces]
    traceNames = reversed(traceNames) # rasterplot top -> bottom

    # create X and Y data for scatter plot
    spike_vecs = [to_numpy(spikeData[trace]) for trace in traceNames]
    x_data = np.concatenate(spike_vecs) # X data is concatenated spike times
    y_data = np.concatenate([np.zeros_like(vec)+j for j, vec in enumerate(spike_vecs)]) # Y-data is trace IDs
    
    # Filter data within given time interval
    if timeRange is not None:
        mask = (x_data > timeRange[0]) & (x_data < timeRange[1])
        x_data = x_data[mask]
        y_data = y_data[mask]
    
    # Plot data as scatter plot
    fig, ax = plt.subplots(figsize=plot_kwargs.get('figsize', None))
    # ax.scatter(x_data, y_data, s=4, c=color, lw=0, marker=marker)
    ax.plot(x_data, y_data, **plot_kwargs) # crispest lines
    ax.set_xlabel('time (ms)')
    
    # X-Y limits
    if timeRange is not None:
        ax.set_xlim(timeRange)

    # X-ticks labels examples
    # https://matplotlib.org/examples/ticks_and_spines/ticklabels_demo_rotation.html
    # https://matplotlib.org/devdocs/gallery/ticks_and_spines/tick_labels_from_values.html
    ax.set_yticklabels(range(len(spike_vecs)), traceNames, rotation='horizontal')
    # plt.margins(0.2) # Pad margins so that markers don't get clipped by the axes
    fig.subplots_adjust(left=0.15) # Tweak spacing to prevent clipping of tick-labels

    # Cosmetics
    ax.grid(showGrid, axis='x')
    
    if title:
        fig.suptitle(title) # suptitle() is Fig title, title() is ax title

    if showFig:
        plt.show(block=False)

    return fig, ax


def test_plotRaster():
    """
    Test for plotRaster() function
    """
    trace_names = ['trace'+str(i) for i in range(12)]
    trace_dict = dict([(tname, h.Vector(np.random.uniform(0,2000,50))) for tname in trace_names])
    plotRaster(trace_dict, timeRange=(0,2000))


def plot_connectivity_matrix(
        W, px=10, py=10, 
        pop0='A', pop1='B', 
        cmap='Oranges', show=True, figsize=(8,6), title=None):
    """
    Plot connectivity matrix given as string.

    @param  W : np.array
            2D matrix containing connection weights

    @param  px, py : int
            Number of cells per (sub)population in 1st/2nd dimension of W.
            Used to draw gridlines.
    """
    from matplotlib import patches

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    fig.suptitle(title)
    # # fig.subplots_adjust(left=0.02) # Less space on left
    # fig.subplots_adjust(right=0.98) # Less space on right
    # fig.subplots_adjust(top=0.96) # Less space on bottom
    # # fig.subplots_adjust(bottom=0.02) # Less space on bottom
    # fig.subplots_adjust(wspace=0) # More space between
    # fig.subplots_adjust(hspace=0) # More space between

    # Plot matrix as image
    colormap = plt.get_cmap(cmap) # see https://matplotlib.org/examples/color/colormaps_reference.html
    plt.imshow(W, interpolation='none', cmap=colormap)

    # Plot grid lines
    y_popsize, x_popsize = W.shape
    y_nticks = y_popsize/py + 1
    x_nticks = x_popsize/px + 1
    yticks_pos = np.arange(y_nticks)*py
    xticks_pos = np.arange(x_nticks)*px
    for p in yticks_pos:
        # Plot gridlines (population boundaries)
        # plt.plot(np.array([0, x_popsize])-0.5, 
        #          np.array([p, p])-0.5,
        #          'k-', linewidth=1.0, snap=True)
        # plt.plot(np.array([p, p])-0.5,
        #          np.array([0, y_popsize])-0.5, 
        #          'k-', linewidth=1.0, snap=True)
        # Add rectangles on diagonal
        ax.add_patch(patches.Rectangle((p-0.5, p-0.5),
                                        px, # Width
                                        py, # Height
                                        facecolor="none",
                                        edgecolor='g',
                                        linewidth="1"))
    # Grid instead of manual gridlines
    plt.grid(True)

    # Configure the x and y axis
    ax.set_xticks(xticks_pos-0.5)
    ax.set_yticks(yticks_pos-0.5)
    ax.set_xticklabels(xticks_pos)
    ax.set_yticklabels(yticks_pos)
    ax.xaxis.set_ticks_position('top')

    ax.set_xlabel('{} cell index'.format(pop1))
    ax.set_ylabel('{} cell index'.format(pop0))
    
    plt.xlim(-0.5, x_popsize - 0.5)
    plt.ylim(y_popsize - 0.5 ,-0.5)

    # Add color bar to measure weights
    plt.clim(0, abs(W).max())
    plt.colorbar()
    if show:
        plt.show(block=False)

if __name__ == '__main__':
    test_plotRaster()
