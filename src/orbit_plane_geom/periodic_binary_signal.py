import numpy
import bisect
import pylab

class PBSignal(object):

    def __init__(self, initial_val, transitions, period):
        assert initial_val == 0 or initial_val == 1
        self.initial_val = initial_val
        assert len(transitions)%2 == 0, "the number of transitions for a periodic signal must be even"
        assert len(set(transitions)) == len(transitions), "transitions must be unique"
        self.transitions = numpy.array(sorted(transitions))
        assert all(self.transitions < period)
        self.period = period


    @classmethod
    def from_random_transitions(cls, num_transitions, period):
        transitions = period*numpy.random.random(num_transitions)
        initial_val = round(numpy.random.random())
        return cls(initial_val, transitions, period)


    @classmethod
    def from_or_of_signals(cls, signal1, signal2):
        # TODO: handle integer multiple periods
        period = signal1.period
        assert signal2.period == period
        #
        transitions1 = list(signal1.transitions)
        transitions2 = list(signal2.transitions)
        #
        initial_val = signal1.initial_val or signal2.initial_val
        cur_val = initial_val
        transitions = []
        possible_transitions = numpy.array(sorted(set(transitions1) | set(transitions2)))
        for t in possible_transitions:
            # v1
            if t in transitions1:
                idx = transitions1.index(t)
                idx_is_even = (idx%2 == 0)
                if idx_is_even:
                    v1 = (signal1.initial_val + 1)%2
                else:
                    v1 = signal1.initial_val
            else:
                v1 = signal1.get_val(t)
            # v2
            if t in transitions2:
                idx = transitions2.index(t)
                idx_is_even = (idx%2 == 0)
                if idx_is_even:
                    v2 = (signal2.initial_val + 1)%2
                else:
                    v2 = signal2.initial_val
            else:
                v2 = signal2.get_val(t)
            #
            if cur_val == 0 and (v1 == 1 or v2 == 1):
                transitions.append(t)
                cur_val = 1
            elif cur_val == 1 and v1 == 0 and v2 == 0:
                transitions.append(t)
                cur_val = 0
        return cls(initial_val, transitions, period)


    @classmethod
    def from_shifted_signal(cls, signal, shift):
        period = signal.period
        initial_val = signal.get_val(shift)
        transitions = (signal.transitions - shift) % period
        return cls(initial_val, transitions, period)


    def get_val(self, t):
        # early out
        if len(self.transitions) == 0:
            return self.initial_val
        #
        t = t%self.period
        idx = bisect.bisect_left(self.transitions, t)
        idx_is_even = (idx%2 == 0)
        if idx_is_even:
            val = self.initial_val
        else:
            val = (self.initial_val + 1)%2
        return val


    def get_time_till_on(self, t):
        # early outs
        if self.get_val(t) == 1:
            return 0
        if len(self.transitions) == 0:
            assert self.initial_val == 0
            return None
        #
        t = t%self.period
        # find the next_transition
        idx = bisect.bisect_left(self.transitions, t)%len(self.transitions)
        # this transition has to be an on transition
        # since we are currently off
        dt = (self.transitions[idx] - t)%self.period
        return dt


    def get_on_durations(self):
        if self.initial_val == 1:
            off_idxs = numpy.arange(0, len(self.transitions), 2)
        else:
            off_idxs = numpy.arange(1, len(self.transitions), 2)
        durations = []
        for off_idx in off_idxs:
            on_idx = (off_idx - 1) % len(self.transitions)
            duration = (self.transitions[off_idx] - self.transitions[on_idx])%self.period
            durations.append(duration)
        return numpy.array(durations)


    def get_off_durations(self):
        if self.initial_val == 0:
            on_idxs = numpy.arange(0, len(self.transitions), 2)
        else:
            on_idxs = numpy.arange(1, len(self.transitions), 2)
        durations = []
        for on_idx in on_idxs:
            off_idx = (on_idx - 1) % len(self.transitions)
            duration = (self.transitions[on_idx] - self.transitions[off_idx])%self.period
            durations.append(duration)
        return numpy.array(durations)


    def plot(self, ax=None, offset=0):
        if ax is None:
            ax = pylab.subplot(1,1,1)
        xs = [0]
        ys = [self.initial_val]
        for t in self.transitions:
            xs = xs + [t,t]
            ys = ys + [ys[-1], (ys[-1] + 1)%2]
        xs.append(1.0)
        ys.append(ys[-1])
        ax.plot(xs, numpy.array(ys)+offset)
        return ax


######
# tests
######

def test_basic():
    my_signal = PBSignal(initial_val=0, transitions=[0.1, 0.35, 0.6, 0.85], period=1)
    #
    assert my_signal.get_val(0.05) == 0
    assert my_signal.get_val(0.25) == 1
    assert my_signal.get_val(0.4) == 0
    assert my_signal.get_val(0.7) == 1
    assert my_signal.get_val(0.9) == 0
    assert my_signal.get_val(1.05) == 0
    #
    assert numpy.isclose(my_signal.get_time_till_on(0.05), 0.05)
    assert numpy.isclose(my_signal.get_time_till_on(0.25), 0)
    assert numpy.isclose(my_signal.get_time_till_on(0.4), 0.2)
    assert numpy.isclose(my_signal.get_time_till_on(0.7), 0)
    assert numpy.isclose(my_signal.get_time_till_on(0.9), 0.2)
    assert numpy.isclose(my_signal.get_time_till_on(1.05), 0.05)
    #
    assert numpy.allclose(my_signal.get_on_durations(), [0.25, 0.25])
    assert numpy.allclose(my_signal.get_off_durations(), [0.25, 0.25])


def test_on_off_sum():
    my_signal = PBSignal.from_random_transitions(num_transitions=10, period=1)
    on_durations = my_signal.get_on_durations()
    off_durations = my_signal.get_off_durations()
    assert numpy.isclose(sum(on_durations) + sum(off_durations), my_signal.period)


def test_from_or_of_signals(doplot=True):
    signal1 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    signal2 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    signal3 = PBSignal.from_or_of_signals(signal1, signal2)
    #
    if doplot:
        import pylab
        ax = pylab.subplot(1,1,1)
        signal1.plot(ax, offset=0)
        signal2.plot(ax, offset=-2)
        signal3.plot(ax, offset=-4)
        pylab.show
    #
    for t in numpy.random.random(1000):
        if signal1.get_val(t) == 0 and signal2.get_val(t) == 0:
            assert signal3.get_val(t) == 0
        else:
            assert signal3.get_val(t) == 1


def test_from_shifted_signal():
    signal1 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    signal2 = PBSignal.from_shifted_signal(signal=signal1, shift=numpy.random.random())
    on_durations1 = sorted(signal1.get_on_durations())
    on_durations2 = sorted(signal2.get_on_durations())
    assert numpy.allclose(on_durations1, on_durations2)
    off_durations1 = sorted(signal1.get_off_durations())
    off_durations2 = sorted(signal2.get_off_durations())
    assert numpy.allclose(off_durations1, off_durations2)


def test_signal_mean():
    signal1 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    # theory
    mu_theory = sum(signal1.get_on_durations())/signal1.period
    # experiment
    nsamps = 10000
    samps = []
    for i in range(nsamps):
        t = i*signal1.period/nsamps
        samps.append(signal1.get_val(t))
    mu_exp = numpy.mean(samps)
    print(mu_theory, mu_exp)
    assert numpy.isclose(mu_theory, mu_exp, rtol=0.01)


def test_signal_var():
    signal1 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    # theory
    mu_theory = sum(signal1.get_on_durations())/signal1.period
    var_theory = mu_theory - mu_theory*mu_theory
    # experiment
    nsamps = 10000
    samps = []
    samps2 = []
    for i in range(nsamps):
        t = i*signal1.period/nsamps
        val = signal1.get_val(t)
        samps.append(val)
        samps2.append(val*val)
    mu_exp = numpy.mean(samps)
    var_exp = numpy.mean(samps2) - mu_exp**2
    print(var_theory, var_exp)
    assert numpy.isclose(var_theory, var_exp, rtol=0.01)


def test_expected_time_till_on1():
    signal1 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    # theory: sum of the on times divided by the period
    mu_theory = 0.5*sum(signal1.get_off_durations()**2)/signal1.period
    # experiment
    nsamps = 10000
    samps = []
    for i in range(nsamps):
        t = i*signal1.period/nsamps
        tto = signal1.get_time_till_on(t)
        samps.append(tto)
    mu_exp = numpy.mean(samps)
    print(mu_theory, mu_exp)
    assert numpy.isclose(mu_theory, mu_exp, rtol=0.01)


def test_total_time_off_stats():
    # expected total off time for two signals with random phase
    signal1 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    off_duration1s = signal1.get_off_durations()
    total_off_time1 = sum(off_duration1s)
    #
    signal2 = PBSignal.from_random_transitions(num_transitions=10, period=1)
    off_duration2s = signal2.get_off_durations()
    total_off_time2 = sum(off_duration2s)
    #
    total_off_time3s = []
    nsamp = 1000
    for i in range(nsamp):
        offset = i*signal2.period/nsamp
        signal2o = PBSignal.from_shifted_signal(signal=signal2, shift=offset)
        signal3 = PBSignal.from_or_of_signals(signal1, signal2o)
        total_off_time3 = sum(signal3.get_off_durations())
        total_off_time3s.append(total_off_time3)
    mean3_experiment = numpy.mean(total_off_time3s)
    mean3_theory = total_off_time1*total_off_time2
    pct_error = (mean3_theory - mean3_experiment)/mean3_experiment
    assert abs(pct_error) < 1e-5


def do_test_expected_time_till_on2():
    nsamp = 100000
    da = numpy.random.random()
    db = da*numpy.random.random()
    while da+db > 1:
       da = numpy.random.random()
       db = da*numpy.random.random()
    p = 1
    ttos = []
    count1 = 0
    count2 = 0
    count3 = 0
    tto1s = []
    tto2s = []
    tto3s = []
    for i in range(nsamp):
        print(i)
        t1 = numpy.random.random()
        t2 = t1+da
        if t2 > 1:
            signal1 = PBSignal(initial_val=1, transitions=[t2%p, t1], period=p)
        else:
            signal1 = PBSignal(initial_val=0, transitions=[t1, t2], period=p)
        t1 = numpy.random.random()
        t2 = t1+db
        if t2 > 1:
            signal2 = PBSignal(initial_val=1, transitions=[t2%p, t1], period=p)
        else:
            signal2 = PBSignal(initial_val=0, transitions=[t1, t2], period=p)
        signal3 = PBSignal.from_or_of_signals(signal1, signal2)
        t = numpy.random.random()
        tto = signal3.get_time_till_on(t)
        ttos.append(tto)
        #
        on_durations = signal3.get_on_durations()
        if len(on_durations) == 2:
           count1 += 1
           tto1s.append(tto)
        else:
           assert len(on_durations) == 1
           on_duration = on_durations[0]
           if numpy.isclose(on_duration, da):
               count2 += 1
               tto2s.append(tto)
           else:
               assert on_duration > da
               count3 += 1
               tto3s.append(tto)
    p1 = (p-da-db)/p
    tto1 = (p-da-db)**2/(3*p)
    p2 = (da-db)/p
    tto2 = (p-da)**2/(2*p)
    p3 = 2*db/p
    tto3 = (p-da-db/2)**2/(2*p)
    #
    print(numpy.mean(ttos))
    print(p1*tto1 + p2*tto2 + p3*tto3)



'''
n = 1000000
c = 0
for i in range(n):
    r = 0.5*numpy.random.random()
    s = (1-r)*numpy.random.random()
    c = c+s
print(c/n)


from sympy import *

p, da, db = symbols("p da db")

p1 = (p-da-db)/p
tto1 = (p-da-db)**2/(3*p)
p2 = (da-db)/p
tto2 = (p-da)**2/(2*p)
p3 = 2*db/p
tto3 = (p-da-db/2)**2/(2*p)

result = simplify(p1*tto1 + p2*tto2 + p3*tto3)
'''
