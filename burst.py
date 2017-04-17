"""burst.py
module to create single pulses for FRB type searches
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
import PSS_utils as utils

class Burst(object):
    def __init__(self, Signal_in, burst_width = 2048, amplitude=1, DM_broadening=False): #period in milliseconds
        """burst_width given in microseconds"""
          self.Signal_in = Signal_in
          self.signal = self.Signal_in.signal
          self.f0 = self.Signal_in.f0
          self.bw = self.Signal_in.bw
          self.Nf = self.Signal_in.Nf
          self.Nt = self.Signal_in.Nt
          self.TotTime = self.Signal_in.TotTime
          self.TimeBinSize = self.TotTime/self.Nt
          self.Time_location = self.Nt//2
          self.burst_width_time = burst_width
          self.burst_width_bins = self.burst_width_time // self.TimeBinSize
          self.phase = np.linspace(0., 1., self.burst_width_bins)
          self.profile = 1./np.sqrt(2.*np.pi)/0.05 * np.exp(-0.5 * ((self.phase-0.25)/0.05)**2)
          self.BurstDict = dict(Profile="gaussian", peak=0.25, width=0.05, amplitude=1.)


    def draw_intensity_pulse(self, reps):
        """draw_intensity_pulse(pulse)
        draw a single pulse as bin by bin random process (gamma distr) from input template
        """
        pr = 100*np.tile(self.profile, reps)
        pulse = np.random.gamma(4., pr/4.) #pull from gamma distribution


        return pulse

    def draw_voltage_pulse(self):
        """draw_voltage_pulse(pulse)
        draw a single pulse as bin by bin random process (normal distr) from input template
        """
        pr = self.profile
        L = len(pr)
        pulse = pr * np.random.normal(0, 1 , L) #pull from gaussian distribution

        return pulse


    #First profiles
    def gauss_template(self, peak=0.25, width=0.05, amp=1.):

        #TODO: error checking for array length consistency?
        #TODO: if any param is a not array, then broadcast to all entries of other arrays?

        try: # is this an array
            peak = np.array(peak)
            width = np.array(width)
            amp = np.array(amp)
            self.PulsarDict["amplitude"] = amp
            self.PulsarDict["Profile"] = "multiple gaussians"
            amp = amp/amp.sum()  # normalize sum
            profile = np.zeros(self.burst_width_bins)
            for ii in range(amp.size):
                norm = amp[ii]/np.sqrt(2.*np.pi)/width[ii]
                self.profile += norm * np.exp(-0.5 * ((self.phase-peak[ii])/width[ii])**2)
        except: # one gaussian
            norm = 1./np.sqrt(2.*np.pi)/width
            self.profile = norm * np.exp(-0.5 * ((self.phase-peak)/width)**2)
            self.PulsarDict["amplitude"] = amp
            self.PulsarDict["Profile"] = "gaussian"

        self.PulsarDict["peak"] = peak
        self.PulsarDict["width"] = width

    def user_template(self,template):
        # Function to make any given 1-dimensional numpy array into the profile
        #TODO Allow other input files
        #TODO Adds error messages if not the correct type of file.
        self.PulsarDict["Profile"] = "user_defined"
        #TODO Add I/O for adding attributes for user defined templates.
        self.PulsarDict["peak"] = "None"
        self.PulsarDict["width"] = "None"
        self.PulsarDict["amplitude"] = "None"
        #self.nBinsPeriod = len(template)
        #self.profile = template
        self.nBinsTemplate = len(template)

        if self.nBinsTemplate==self.burst_width_bins:
            self.profile = template

        elif self.nBinsTemplate > self.burst_width_bins:
            self.profile = utils.rebin(template, self.burst_width_bins)
            print("User supplied template has been downsampled.")
            print("Input array length= ", self.nBinsTemplate,". Pulse template length= ",self.profile.size,".")

        else:
            TempPhase = np.linspace(0,1,len(template))
            ProfileFcn = sp.interpolate.interp1d(TempPhase, template, kind='cubic', bounds_error=True)
            self.profile = ProfileFcn(self.phase)
            print("User supplied template has been interpolated using a cubic spline.")
            print("Input array length was ", self.nBinsTemplate," bins. New pulse template length is ",self.profile.size,".")

        self.MinCheck = np.amin(self.profile)
        if self.MinCheck < 0 :
            self.profile = np.where(self.profile > 0, self.profile, self.profile-self.MinCheck)


    def burst(self, SignalType = "intensity"):
        #Function that makes a pulse using the defined profile template
        pulseType = {"intensity":"draw_intensity_pulse", "voltage":"draw_voltage_pulse"}
        pulseTypeMethod = getattr(self, pulseType[SignalType])

        if self.DM_broadening==True:
            profile_table = np.zeros((self.Nf, self.burst_width_bins))
            #for ii in range(self.NF):
                #t_h_w = utils.top_hat_width(sub_band_width, sub_bandwidth_center, DM)#Need to figure out when to include DM
            self.signal[ii,self.Time_location:self.Time_location+len(self.profile)] = np.tile(pulseTypeMethod(1),(self.Nf,1))
        else:
            self.signal[:,self.Time_location:self.Time_location+len(self.profile)] = np.tile(pulseTypeMethod(1),(self.Nf,1))

        self.BurstDict["SignalType"] = SignalType

        self.Signal_in.MetaData.AddInfo(self.BurstDict)