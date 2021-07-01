"""
Utilities for recording NEURON traces using PyNN.

@author     Lucas Koelman

@date       20/03/2018
"""

import re
from datetime import datetime

from pyNN import recording
from pyNN.neuron.recording import Recorder # recordable_pattern
from pyNN.neuron import simulator

import numpy
import neo
import quantities as pq

import neuron
h = neuron.h

from bgcellmodels.common import nrnutil

logger = recording.logger


class TraceSpecRecorder(Recorder):
    """
    Extension of pyNN Recorder class for NEURON simulator that understands
    trace specifications in the format used by NetPyne

    @see        Based on NetPyne's Cell.RecordTraces in file:
                https://github.com/Neurosim-lab/netpyne/blob/master/netpyne/cell.py


    Example
    -----

        >>> from pyNN.neuron import Population
        >>> from bgcellmodels.extensions.pynn.recording import TraceSpecRecorder
        >>> Population._recorder_class = TraceSpecRecorder
        >>> ... (recording setup code)


    Attributes
    ----------

    @attr   sampling_interval : float
            Sampling interval for recording

    @attr   recorded : dict[str, set(int)]
            Cell IDs of recorded cells for each trace name

    @attr   trace_opts : dict[str, dict[str, object]]
            Options for each trace.


    Developer Notes
    ---------------

    See common Population and Recorder class for method chains involved
    in recording and writing traces:

    https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/common/populations.py
    https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/recording/__init__.py

        - sampling_interval is common to all traces/signals of the recorder

        - PyNN.common.populations.BasePopulation.write_data(io, variables, gather,
                                                            clear, annotations)
            `-> PyNN.recording.Recorder.write(...)
                `-> data = self.get(variables, gather, ...)
                    `-> segment = self._get_current_segment(...)
                    `-> self.clear()
                `-> neo_io.write_block(data)
    """

    def __init__(self, *args, **kwargs):
        self.trace_opts = {}
        super(TraceSpecRecorder, self).__init__(*args, **kwargs)


    def record(self, variables, ids,
               sampling_interval=None,
               write_interval=None):
        """
        Override the default record() method to accept trace specifications
        in the format used by NetPyne.

        @override   pyNN.recording.Recorder.record()


        @param      variables : iterable(tuple(str, <dict or str>))

                    Any iterable where the first element of each item
                    is the trace name, and the second element the trace
                    specification. The trace specification can be either a
                    string like the default PyNN variable names, or a dict
                    that specifies the trace in NetPyne format.

        @note       The original graph is as follows:
                    -> Recorder.record() @ /pyNN/recording/__init__.py
                    -> Recorder._record() @ /pyNN/neuron/recording.py
                    -> Recorder._record_state_variable() @ /pyNN/neuron/recording.py

        Example
        -----

        >>> trace_specs = {
        >>>     'GP_cai':{'sec':'soma[0]','loc':0.5,'var':'cai'},
        >>>     'GP_ainf':{'sec':'soma[0]','loc':0.5,'mech':'gpRT','var':'a_inf'},
        >>>     'GP_r':{'sec':'soma[0]','loc':0.5,'mech':'gpRT','var':'r'},
        >>> }
        >>> pop.record(trace_specs.items())

        """
        logger.debug('Recorder.record(<cell ids: {}>)'.format(ids))

        # if sampling_interval is not None:
        #     if ((sampling_interval != self.sampling_interval) and (len(self.recorded) > 0) and
        #             not (len(self.recorded) == 1 and 'spikes' in self.recorded)):
        #         raise ValueError("All neurons in a population must be recorded with the same sampling interval.")

        ids = set([id for id in ids if id.local])

        for variable in variables:
            if isinstance(variable, str):
                trace_name = variable
            else:
                trace_name, trace_spec = variable
            # if not self.population.can_record(trace_spec):
            #     raise errors.RecordingError(trace_spec, self.population.celltype)

            # Get cells that aren't recording this trace yet
            new_ids = ids.difference(self.recorded[trace_name])

            # Set the sampling interval for this trace
            sampling_interval = sampling_interval or self._simulator.state.dt
            if trace_name in self.trace_opts:
                if self.trace_opts[trace_name]['sampling_interval'] != sampling_interval:
                    raise ValueError("All neurons in a population must be recorded with the same sampling interval.")
            else:
                self.trace_opts[trace_name] = {}
            self.trace_opts[trace_name]['sampling_interval'] = sampling_interval

            self.recorded[trace_name] = self.recorded[trace_name].union(ids)
            self._record(variable, new_ids)


    def _record(self, variable, new_ids):
        """
        Set up NEURON recordings for givell cell ids.
        """
        if variable == 'spikes':
            for id in new_ids:
                if id._cell.rec is not None:
                    id._cell.rec.record(id._cell.spike_times)
        else:
            for id in new_ids:
                # Override to pass ID instead of cell
                self._record_state_variable(id, variable)


    def _get_trace_opt(self, trace_name, opt_name):
        """
        Get recording options for given trace name (variable).
        """
        # trace_name might be Synapse1 with Synapse{:d} present as key
        if trace_name not in self.trace_opts:
            match = re.match(r'(\w+?)(\d+)$', trace_name)
            if not match:
                raise KeyError(trace_name)
            trace_prefix = match.group(1)
            trace_template = next((k for k in self.trace_opts.keys() if
                                   re.match(trace_prefix + r'\{.+?\}', k)), None)
            if not trace_template:
                raise KeyError(trace_name)
            return self.trace_opts[trace_template][opt_name]
        return self.trace_opts[trace_name][opt_name]


    def _record_state_variable(self, id, variable):
        """
        Record the variable specified by the object 'variable'.

        @override   pyNN.neuron.recording.Recorder._record_state_variable()

        @param      cell : ID._cell
                    Instantiated cell model (CellType.model(**params))

        @param      variable : str OR tuple(str, <str or dict>)
                    A trace specifier consisting of the trace name as first
                    element and full trace specification as second element.
        """
        cell = id._cell
        if isinstance(variable, str):
            trace_name, trace_spec = variable, variable
        else:
            trace_name, trace_spec = variable

        recorded = False
        sampling_interval = self._get_trace_opt(trace_name, 'sampling_interval')

        # First try to interpret spec as PyNN format (string)
        if hasattr(cell, 'recordable') and trace_spec in cell.recordable:
            hoc_var = cell.recordable[trace_spec]
        elif trace_spec == 'v':
            hoc_var = cell.source_section(0.5)._ref_v  # or use "seg.v"?
        elif trace_spec == 'gsyn_exc':
            hoc_var = cell.esyn._ref_g
        elif trace_spec == 'gsyn_inh':
            hoc_var = cell.isyn._ref_g
        elif trace_spec == 'spikes':
            cell.rec.record(cell.spike_times)
            return # was implemented in _record() -> don't execute rest of body

        elif trace_spec == 'lfp':
            # assume cell has attribute 'lfp_tracker'
            vec = h.Vector()
            pp = cell.lfp_tracker.summator
            hoc_ref = cell.lfp_tracker.summator._ref_summed
            vec.record(pp, hoc_ref, sampling_interval)
            cell.traces[trace_name] = vec
            recorded = True

        elif isinstance(trace_spec, str):
            source, var_name = self._resolve_variable(cell, trace_spec)
            hoc_var = getattr(source, "_ref_%s" % var_name)

        elif 'sec' in trace_spec.keys():
            # spec is in NetPyne format
            hoc_obj = cell.resolve_section(trace_spec['sec'])
            vec, marker = self._record_trace(hoc_obj, trace_spec, sampling_interval,
                            threshold=cell.get_threshold())
            if marker is not None:
                _pp_markers = getattr(cell, '_pp_markers', [])
                _pp_markers.append(marker)
                cell._pp_markers = _pp_markers
            cell.traces[trace_name] = vec
            recorded = True

        elif 'secs' in trace_spec.keys():
            # spec is in NetPyne format
            hoc_objs = cell.resolve_section(trace_spec['secs'], multiple=True)
            for i, hoc_obj in enumerate(hoc_objs):
                vec, marker = self._record_trace(hoc_obj, trace_spec, sampling_interval,
                                threshold=cell.get_threshold())
                if marker is not None:
                    _pp_markers = getattr(cell, '_pp_markers', [])
                    _pp_markers.append(marker)
                    cell._pp_markers = _pp_markers

                numbered_trace_name = trace_name.format(i)
                cell.traces[numbered_trace_name] = vec
                self.recorded[numbered_trace_name] = \
                    self.recorded[numbered_trace_name].union((id,))

            self.recorded[trace_name] = set() # remove unformatted trace name
            recorded = True

        elif 'syn' in trace_spec.keys():
            # trace spec
            #   'Syn{:d}' : {'syn': 'Exp2Syn[:]', 'var':'g'}
            #   'Syn{:d}' : {'syn': 'Exp2Syn[0]', 'var':'g'}
            pp_list = cell.resolve_synapses(trace_spec['syn'])
            for i, pp in enumerate(pp_list):
                vec = h.Vector()
                hoc_ptr = getattr(pp, '_ref_'+trace_spec['var'])
                vec.record(pp, hoc_ptr, sampling_interval)

                numbered_trace_name = trace_name.format(i)
                cell.traces[numbered_trace_name] = vec
                self.recorded[numbered_trace_name] = \
                    self.recorded[numbered_trace_name].union((id,))

            # This is a hack: remove the unformatted trace name from recorded
            # signals for this ID
            self.recorded[trace_name] = set()
            recorded = True

        # Record global variable
        if not recorded:
            vec = h.Vector()
            if sampling_interval == self._simulator.state.dt:
                vec.record(hoc_var)
            else:
                vec.record(hoc_var, sampling_interval)
            cell.traces[trace_name] = vec

        # Record global time variable 't' if not recorded already
        if not cell.recording_time:
            cell.record_times = h.Vector()
            if sampling_interval == self._simulator.state.dt:
                cell.record_times.record(h._ref_t)
            else:
                cell.record_times.record(h._ref_t, sampling_interval)

            cell.recording_time += 1


    def _record_trace(self, hoc_obj, spec, rec_dt, duration=None, threshold=None):
        """
        Record the given traces from section

        @param  threshold : float
                Spike threshold if trace spec contains 'spikes' as its var

        WARNING: For multithreaded execution, section _must_ have POINT_PROCESS
                 associated with it to identify the tread. Hence, a PointProcessMark
                 (see ppmark.mod) will be inserted if no PP present.
        """
        pp_marker = None
        hoc_ptr = None  # pointer to Hoc variable that will be recorded
        pp = None       # Hoc POINT_PROCESS instance
        vec = None      # Hoc.Vector that will be recorded into
        recorded = False

        # Get Section and segment
        if isinstance(hoc_obj, neuron.nrn.Segment):
            seg = hoc_obj
            sec = seg.sec
            if 'loc' in spec:
                seg = sec(spec['loc'])
        else:
            sec = hoc_obj
            seg = sec(spec['loc'])


        if 'loc' in spec: # hoc_obj is Section

            # Get pointer/reference to variable to record
            if spec['var'] == 'spikes':
                pp_marker = h.NetCon(seg._ref_v, None, threshold, 0.0, 0.0, sec=sec)
                vec = h.Vector(0)
                pp_marker.record(vec)
                recorded = True

            # Mechanism RANGE variable, e.g. soma(0.5).hh._ref_gna
            elif 'mech' in spec:
                mech_instance = getattr(seg, spec['mech'])
                hoc_ptr = getattr(mech_instance, '_ref_'+spec['var'])

            # Section RANGE variable, e.g. soma(0.5)._ref_v
            else:
                hoc_ptr = getattr(seg, '_ref_'+spec['var'])

            # find a POINT_PROCESS in segment to improve efficiency
            # seg_pps = seg.point_processes()
            # if any(seg_pps):
            #     pp = seg_pps[0]
            # else:
            if pp_marker is None:
                pp = h.PointProcessMark(seg)
                pp_marker = pp


        elif 'pointp' in spec: # hoc_obj is POINT_PROCESS
            # Look for the point process in Section
            mech_name = spec['pointp']
            seg_pps = seg.point_processes()

            pp = next((hobj for hobj in seg_pps if
                        nrnutil.get_mod_name(hobj)==mech_name), None)

            if pp is None:
                raise ValueError("Could not find point process '{}' "
                    "in segment {}".format(spec['pointp'], seg))

            hoc_ptr = getattr(pp, '_ref_'+spec['var'])


        else: # global vars, e.g. h.t
            hoc_ptr = getattr(h, '_ref_'+spec['var'])


        # Record from pointer into Vector
        if duration is not None:
            vec = h.Vector(duration/rec_dt+1).resize(0)
        else:
            vec = h.Vector()

        if pp and not recorded:
            vec.record(pp, hoc_ptr, rec_dt)
        elif not recorded:
            vec.record(hoc_ptr, rec_dt)

        return vec, pp_marker


    def _get_current_segment(self, filter_ids=None, variables='all', clear=False):
        """
        Rewrite of same method in module pyNN.recording.Recorder to support
        distinct sampling period per trace.
        """
        segment = neo.Segment(name="segment%03d" % self._simulator.state.segment_counter,
                              description=self.population.describe(),
                              rec_datetime=datetime.now())  # would be nice to get the time at the start of the recording, not the end
        variables_to_include = set(self.recorded.keys())
        if variables != 'all':
            variables_to_include = variables_to_include.intersection(set(variables))
        for variable in variables_to_include:
            if variable.startswith('spikes'):
                t_stop = self._simulator.state.t * pq.ms  # must run on all MPI nodes
                sids = sorted(self.filter_recorded(variable, filter_ids))
                data = self._get_spiketimes(sids, variable)
                spiketrains = [
                    neo.SpikeTrain(data.get(int(id),[]),
                                   t_start=self._recording_start_time,
                                   t_stop=t_stop,
                                   units='ms',
                                   source_population=self.population.label,
                                   source_id=int(id),
                                   source_index=self.population.id_to_index(int(id)))
                    for id in sids]

                # FIXME: spike recording
                # if variable == 'spikes':
                #     segment.spiketrains = spiketrains
                # else:
                #     setattr(segment, variable, spiketrains)
                segment.spiketrains.extend(spiketrains)
            else:
                ids = sorted(self.filter_recorded(variable, filter_ids))
                signal_array = self._get_all_signals(variable, ids, clear=clear)
                t_start = self._recording_start_time
                sampling_period = self._get_trace_opt(variable, 'sampling_interval') * pq.ms
                current_time = self._simulator.state.t * pq.ms
                mpi_node = self._simulator.state.mpi_rank  # for debugging
                if signal_array.size > 0:  # may be empty if none of the recorded cells are on this MPI node
                    units = self.population.find_units(variable)
                    source_ids = numpy.fromiter(ids, dtype=int)
                    source_index = numpy.array([self.population.id_to_index(id) for id in ids])
                    signal = neo.AnalogSignal(
                                    signal_array,
                                    units=units,
                                    t_start=t_start,
                                    sampling_period=sampling_period,
                                    name=variable,
                                    source_population=self.population.label,
                                    source_ids=source_ids,
                                    source_indices=source_index)
                    signal.channel_index = neo.ChannelIndex(
                            index=numpy.arange(source_ids.size),
                            channel_ids=source_index)
                    segment.analogsignals.append(signal)
                    logger.debug("%d **** ids=%s, channels=%s", mpi_node, source_ids, signal.channel_index)
                    # assert abs(signal.t_stop - current_time) < 2 * sampling_period + 1e-10 * pq.ms
                    # need to add `Unit` and `RecordingChannelGroup` objects
        return segment


    def _get_spiketimes(self, id, trace_name='spikes'):
        """
        Copied from PyNN master (latest version) to solve bug
        """
        if hasattr(id, "__len__"):
            all_spiketimes = {}
            for cell_id in id:
                if trace_name == 'spikes':
                    spikes = numpy.array(cell_id._cell.spike_times)
                else:
                    spikes = numpy.array(cell_id._cell.traces[trace_name])
                # add OR mask: self._recording_start_time <= spikes
                mask = ((spikes >= self._recording_start_time) &
                        (spikes <= simulator.state.t + 1e-9))
                all_spiketimes[cell_id] = spikes[mask]
            return all_spiketimes
        else:
            spikes = numpy.array(id._cell.spike_times)
            return spikes[spikes <= simulator.state.t + 1e-9]
