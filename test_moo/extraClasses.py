#! /ems/elsc-labs/segev-i/yoni.leibner/anaconda2/bin/ipython
#############################################################
#
#maker: YONI LEIBNER
#
# parameter that compansate for spine area in CM or g_pas
#
#############################################################
from bluepyopt import ephys
from bluepyopt.ephys.parameters import NrnParameter,NrnRangeParameter
from bluepyopt.ephys.parameterscalers import *

import logging
logger = logging.getLogger(__name__)
##########################################################
neuron_start_time = 300

class NrnSegmentSomaDistanceScaler_(ParameterScaler, DictMixin):
    """Scaler based on distance from soma"""
    # SERIALIZED_FIELDS = ('name', 'comment', 'distribution',)

    def __init__(
            self,
            name=None,
            comment='',
            dist_param_names=None,
            dist_thresh_apical=60,
            dist_thresh_basal=60,
            F_factor = 1.9,
            shrinckage_factor = 1.0/0.7):

        """Constructor

        Args:
            name (str): name of this object
            distribution (str): distribution of parameter dependent on distance
                from soma. string can contain `distance` and/or `value` as
                placeholders for the distance to the soma and parameter value
                respectivily
            dist_params (list): list of names of parameters that parametrise
                the distribution. These names will become attributes of this
                object.
                The distribution string should contain these names, and they
                will be replaced by values of the corresponding attributes
        """

        super(NrnSegmentSomaDistanceScaler_, self).__init__(name, comment)

        self.dist_param_names = dist_param_names

        if self.dist_param_names is not None:
            for dist_param_name in self.dist_param_names:
                if dist_param_name not in self.distribution:
                    raise ValueError(
                        'NrnSegmentSomaDistanceScaler: "{%s}" '
                        'missing from distribution string "%s"' %
                        (dist_param_name, None))
                setattr(self, dist_param_name, None)

        self.dist_thresh_apical = dist_thresh_apical
        self.dist_thresh_basal = dist_thresh_basal
        self.F_factor = F_factor
        self.shrinckage_factor=shrinckage_factor


    def scale(self, value, sec, sim=None):
        """Scale a value based on a segment"""
        # TODO soma needs other addressing scheme
        # Initialise origin
        try:
            soma = sec.cell().soma[0]
            sim.neuron.h.distance(0, 0.5, sec=soma)
            L = sec.L / sec.nseg
            distance = sim.neuron.h.distance(1, sec(0.5).x, sec=sec)
            if sec.name().find("spine") > -1:
                spine_factor = 1
            elif sec in sec.cell().apical:
                if distance > self.dist_thresh_apical:
                    spine_factor = self.F_factor
                else:
                    spine_factor = 1  # no change
            elif sec in sec.cell().basal:
                if distance > self.dist_thresh_basal:
                    spine_factor = self.F_factor
                else:
                    spine_factor = 1  # no change
            else:
                spine_factor = 1
            if spine_factor < 1:
                print(sec, "spine factor < 1")
            sim.neuron.h.pop_section()
            # print(sec.name(), distance, value, value*spine_factor)

            return value * spine_factor# * self.shrinckage_factor
        except:
            return value #* self.shrinckage_factor
        # if self.name.find('cm') > -1: return value * spine_factor
        # if self.name.find('g_pas') > -1: return value * spine_factor
        #
        # raise Exception('this should be use only for CM or g_pas values!!!!!!!!!')

    def __str__(self):
        """String representation"""
        return "apic spine start:"+str(self.dist_thresh_apical)+\
               "\nbasal spine start:"+str(self.dist_thresh_basal)+\
               "\n spine factor is:"+str(self.F_factor)


from bluepyopt.ephys.efeatures import EFeature


class EFeatureRDSM(EFeature):
    def __init__(self, T, V, extra_rec=""):
        """

        :param T: the time objective
        :param V: the voltage objective
        """
        super(EFeatureRDSM, self).__init__('RDSM', '')
        # if len(fit_times) < 1:
        #     raise Exception('need time for RDSM fit')
        # self.fit_times = fit_times
        self.T_objective = T
        self.V_objective = V
        self.extra_rec = extra_rec

    def calculate_feature(self, responses, raise_warnings=False):
        """Calculate feature value"""
        import numpy
        if self.extra_rec +'soma.v' not in responses: return 0
        T = numpy.array(responses[self.extra_rec +'soma.v']['time'])
        V = numpy.array(responses[self.extra_rec +'soma.v']['voltage'])
        seatell = numpy.where(T > neuron_start_time)[0][0]
        seatell2 = int(neuron_start_time / (T[1] - T[0]))
        T = T - T[seatell]
        T = T[seatell:]
        V = V[seatell:]
        temp =numpy.sqrt(numpy.sum((self.V_objective - V[:len(self.V_objective)]) ** 2))
        # print(temp)
        if temp < 30:
            a=1
        return numpy.sqrt(numpy.sum((self.V_objective - V[:len(self.V_objective)]) ** 2))
               # + abs(numpy.max(self.V_objective)-numpy.max(V))*100 + abs(peak_pos1-peak_pos2) * dt * 20


    def calculate_score(self, responses, trace_check=False):
        return self.calculate_feature(responses, raise_warnings=trace_check)


class EFeatureHalfWidth(EFeature):
    def __init__(self, T, V):
        """

        :param T: the time objective
        :param V: the voltage objective
        """
        super(EFeatureHalfWidth, self).__init__('HalfWidth', '')
        # if len(fit_times) < 1:
        #     raise Exception('need time for HalfWidth fit')
        # self.fit_times = fit_times
        self.T_objective = T
        self.V_objective = V
        import numpy
        dt = T[1] - T[0]
        pas = V[0]
        peak_V = numpy.max(V)
        mid = pas + abs(peak_V - pas)/2
        # V = numpy.round(V,1)
        up_ = numpy.where(V > mid)[0][0]
        down_ = numpy.where(V[up_] > mid)[0][-1]
        self.halfWidth =(down_ - up_)*dt


    def calculate_feature(self, responses, raise_warnings=False):
        """Calculate feature value"""
        import numpy
        T = numpy.array(responses['soma.v']['time'])
        V = numpy.array(responses['soma.v']['voltage'])
        seatell = numpy.where(T > neuron_start_time)[0][0]
        T = T - T[seatell]
        T = T[seatell:]
        V = V[seatell:]
        dt = T[1]-T[0]
        pas = V[0]
        peak_V = numpy.max(V)
        mid = pas + abs(peak_V - pas)/2
        up_ = numpy.where(V > mid)[0][0]
        down_ = numpy.where(V[up_] > mid)[0][-1]
        return abs(dt*(down_ - up_) - self.halfWidth)


    def calculate_score(self, responses, trace_check=False):
        return self.calculate_feature(responses, raise_warnings=trace_check)



class EFeaturePeakTime(EFeature):
    def __init__(self, T, V):
        """

        :param T: the time objective
        :param V: the voltage objective
        """
        super(EFeaturePeakTime, self).__init__('peak_time', '')
        # if len(fit_times) < 1:
        #     raise Exception('need time for RDSM fit')
        # self.fit_times = fit_times
        self.T_objective = T
        self.V_objective = V

    def calculate_feature(self, responses, raise_warnings=False):
        """Calculate feature value"""
        import numpy
        T = numpy.array(responses['soma.v']['time'])
        V = numpy.array(responses['soma.v']['voltage'])
        seatell = numpy.where(T > neuron_start_time)[0][0]
        T = T - T[seatell]
        T = T[seatell:]
        V = V[seatell:]
        dt = T[1]-T[0]

        peak_pos1=numpy.where(V==numpy.max(V))[0][0]
        peak_pos2=numpy.where(self.V_objective==numpy.max(self.V_objective))[0][0]

        return  abs(peak_pos1-peak_pos2) * dt


    def calculate_score(self, responses, trace_check=False):
        return self.calculate_feature(responses, raise_warnings=trace_check)


class EFeaturePeak(EFeature):
    def __init__(self, T, V, exp_std=1):
        """

        :param T: the time objective
        :param V: the voltage objective
        """
        super(EFeaturePeak, self).__init__('RDSM', '')
        # if len(fit_times) < 1:
        #     raise Exception('need time for RDSM fit')
        # self.fit_times = fit_times
        self.T_objective = T
        self.V_objective = V
        self.exp_std=exp_std

    def calculate_feature(self, responses, raise_warnings=False):
        """Calculate feature value"""
        import numpy
        T = numpy.array(responses['soma.v']['time'])
        V = numpy.array(responses['soma.v']['voltage'])
        seatell = numpy.where(T > neuron_start_time)[0][0]
        T = T - T[seatell]
        T = T[seatell:]
        V = V[seatell:]
        return abs(numpy.max(self.V_objective)-numpy.max(V))


    def calculate_score(self, responses, trace_check=False):
        return self.calculate_feature(responses, raise_warnings=trace_check)/self.exp_std

class EFeatureImpadance(EFeature):
    def __init__(self, Rin, injEnd, injCur, extra_rec="Rin1"):
        """

        :param T: the time objective
        :param V: the voltage objective
        """
        super(EFeatureImpadance, self).__init__('impadance', '')
        # if len(fit_times) < 1:
        #     raise Exception('need time for RDSM fit')
        # self.fit_times = fit_times
        self.Rin = Rin # in M ohm
        self.inj_end = injEnd
        self.inj = injCur
        self.extra_rec=extra_rec+"."

    def calculate_feature(self, responses, raise_warnings=False):
        """Calculate feature value"""
        import numpy
        if self.extra_rec +'soma.v' not in responses: return 0

        T = numpy.array(responses[self.extra_rec+'soma.v']['time'])
        V = numpy.array(responses[self.extra_rec+'soma.v']['voltage'])
        # plt.plot(T,V)
        # plt.show()
        # plt.close()
        end_pos = numpy.where(T < self.inj_end)[0][-1]
        Rin_cur = abs((V[end_pos]-V[-1])/self.inj)
        print(self.Rin, Rin_cur)
        return abs(self.Rin - Rin_cur)


    def calculate_score(self, responses, trace_check=False):
        return self.calculate_feature(responses, raise_warnings=trace_check)

from bluepyopt.ephys import parameterscalers

class NrnSectionParameterPas(NrnParameter, DictMixin):

    """Parameter of a section"""
    SERIALIZED_FIELDS = ('name', 'value', 'frozen', 'bounds', 'param_name',
                         'value_scaler', 'locations', )

    def __init__(
            self,
            name,
            tau,
            Rin,
            value=None,
            frozen=False,
            bounds=None,
            param_name=None,
            value_scaler=None,
            locations=None):
        """Contructor
        Args:
            name (str): name of the Parameter
            value (float): Value for the parameter, required if Frozen=True
            frozen (bool): Whether the parameter can be varied, or its values
            is permently set
            bounds (indexable): two elements; the lower and upper bounds
                (Optional)
            param_name (str): name used within NEURON
            value_scaler (float): value used to scale the parameter value
            locations (list of ephys.locations.Location): locations on which
                to instantiate the parameter
        """

        super(NrnSectionParameterPas, self).__init__(
            name,
            value=value,
            frozen=frozen,
            bounds=bounds)

        self.locations = locations
        self.param_name = param_name
        # TODO value_scaler has to be made more general
        self.value_scaler = value_scaler
        # TODO add a default value for a scaler that is picklable
        if self.value_scaler is None:
            self.value_scaler = parameterscalers.NrnSegmentLinearScaler()
        self.value_scale_func = self.value_scaler.scale
        self.tau=tau
        self.Rin=Rin

    def instantiate(self, sim=None, icell=None):
        """Instantiate"""
        if self.value is None:
            raise Exception(
                'NrnSectionParameter: impossible to instantiate parameter "%s"'
                ' without value' %
                self.name)

        cm = self.value
        Rm = 1000.0 * (self.tau / float(cm))
        Ra = (((self.Rin/1000.0)**2)*4)/Rm

        for location in self.locations:
            iseclist = location.instantiate(sim=sim, icell=icell)
            for section in iseclist:
                setattr(section, "cm",
                        self.value_scale_func(cm, section, sim=sim))
                setattr(section, "g_pas",
                        self.value_scale_func(1.0/Rm, section, sim=sim))
                # setattr(section, "Ra", Ra)

            logger.debug(
                'Set cm in %s to %s\n'
                'Set Rm in %s to %s\n'
                'Set Ra in %s to %s'
                ,location,cm,location,Rm,location,Ra)

    def __str__(self):
        """String representation"""
        return '%s: %s cm = %s' % (self.name,
                                   [str(location)
                                    for location in self.locations],
                                   self.value if self.frozen else self.bounds)


class NrnNetstimWeightParameter(NrnParameter, DictMixin):

    """Parameter of a section"""
    SERIALIZED_FIELDS = ('name', 'value', 'frozen', 'bounds', 'param_name',
                         'value_scaler', 'locations', )

    def __init__(
            self,
            name,
            value=None,
            frozen=False,
            bounds=None,
            locations=None,
            param_name=None,
            reletive_strength = None):
        """Constructor

        Args:
            name (str): name of the Parameter
            value (float): Value for the parameter, required if Frozen=True
            frozen (bool): Whether the parameter can be varied, or its values
            is permently set
            bounds (indexable): two elements; the lower and upper bounds
                (Optional)
            locations: an iterator of the point process locations you want to
                       set the parameters of
            param_name (str): name of parameter used within the point process
        """

        super(NrnNetstimWeightParameter, self).__init__(
            name,
            value=value,
            frozen=frozen,
            bounds=bounds)

        self.locations = locations
        self.param_name = param_name
        self.total_duration=0 # this is a dummy variable
        self.reletive_strength = reletive_strength


    def instantiate(self, sim=None, icell=None):
        """Instantiate"""
        if self.value is None:
            raise Exception(
                'NrnSectionParameter: impossible to instantiate parameter "%s"'
                ' without value' %
                self.name)
        i=0
        for location in self.locations:
            all = location.instantiate(sim=sim, icell=icell)
            for pprocess in location.connections:
                if self.reletive_strength is None:
                    location.connections[pprocess][0][0].weight[0] = self.value
                else:
                    location.connections[pprocess][0][0].weight[0] = self.value * self.reletive_strength[i]
                i+=1
                logger.debug(
                        'Set %s to %s for point process',
                        self.param_name,
                        self.value)
        if not i == len(self.reletive_strength):
            print('error in reletive_strength, the length of reletive_strengths is too long')

    def __str__(self):
        """String representation"""
        return '%s: %s = %s' % (self.name,
                                self.param_name,
                                self.value if self.frozen else self.bounds)


# from bluepyopt.ephys.protocols import Protocol
import collections
from bluepyopt.ephys import simulators, locations

class SweepProtocolRin2(ephys.protocols.Protocol):

    """Sweep protocol"""

    def __init__(
            self,
            name=None,
            stimuli=None,
            recordings=None,
            cvode_active=None,
            dt=None,
            e_pas = None):
        """Constructor

        Args:
            name (str): name of this object
            stimuli (list of Stimuli): Stimulus objects used in the protocol
            recordings (list of Recordings): Recording objects used in the
                protocol
            cvode_active (bool): whether to use variable time step
        """

        super(SweepProtocolRin2, self).__init__(name)
        self.stimuli = stimuli
        self.recordings = recordings
        self.cvode_active = cvode_active

        self.dt=dt
        self.e_pas=e_pas

    @property
    def total_duration(self):
        """Total duration"""

        return max([stimulus.total_duration for stimulus in self.stimuli])

    def subprotocols(self):
        """Return subprotocols"""

        return collections.OrderedDict({self.name: self})

    def _run_func(self, cell_model, param_values, sim=None):
        """Run protocols"""
        print("in _run_func")
        try:
            cell_model.freeze(param_values)
            cell_model.instantiate(sim=sim)

            self.instantiate(sim=sim, icell=cell_model.icell)

            try:
                for sec in cell_model.icell.all:
                    sec.e_pas = self.e_pas
                sim.run(self.total_duration, cvode_active=self.cvode_active, dt=self.dt)
            except (RuntimeError, simulators.NrnSimulatorException):
                logger.debug(
                    'SweepProtocol: Running of parameter set {%s} generated '
                    'an exception, returning None in responses',
                    str(param_values))
                responses = {recording.name:
                             None for recording in self.recordings}
            else:
                responses = {
                    recording.name: recording.response
                    for recording in self.recordings}

            self.destroy(sim=sim)

            cell_model.destroy(sim=sim)

            cell_model.unfreeze(param_values.keys())

            return responses
        except BaseException:
            import sys
            import traceback
            raise Exception(
                "".join(
                    traceback.format_exception(*sys.exc_info())))

    def run(self, cell_model, param_values, sim=None, isolate=None):
        """Instantiate protocol"""

        if isolate is None:
            isolate = True

        if isolate:
            def _reduce_method(meth):
                """Overwrite reduce"""
                return (getattr, (meth.__self__, meth.__func__.__name__))

            import copyreg
            import types
            copyreg.pickle(types.MethodType, _reduce_method)

            import multiprocessing
            print("in SweepProtocolRin2 tring to run function named: ", self._run_func.__name__)
            pool = multiprocessing.Pool(1, maxtasksperchild=1)
            responses = pool.apply(
                self._run_func,
                kwds={
                    'cell_model': cell_model,
                    'param_values': param_values,
                    'sim': sim})

            pool.terminate()
            pool.join()
            del pool
        else:
            responses = self._run_func(
                cell_model=cell_model,
                param_values=param_values,
                sim=sim)

        return responses

    def instantiate(self, sim=None, icell=None):
        """Instantiate"""

        for stimulus in self.stimuli:
            stimulus.instantiate(sim=sim, icell=icell)

        for recording in self.recordings:
            try:
                recording.instantiate(sim=sim, icell=icell)
            except locations.EPhysLocInstantiateException:
                logger.debug(
                    'SweepProtocol: Instantiating recording generated '
                    'location exception, will return empty response for '
                    'this recording')

    def destroy(self, sim=None):
        """Destroy protocol"""

        for stimulus in self.stimuli:
            stimulus.destroy(sim=sim)

        for recording in self.recordings:
            recording.destroy(sim=sim)

    def __str__(self):
        """String representation"""

        content = '%s:\n' % self.name

        content += '  stimuli:\n'
        for stimulus in self.stimuli:
            content += '    %s\n' % str(stimulus)

        content += '  recordings:\n'
        for recording in self.recordings:
            content += '    %s\n' % str(recording)

        return content
