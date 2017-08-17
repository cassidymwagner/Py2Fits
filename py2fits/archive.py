'''
Jeffrey Hazboun 2017

based on archive class in PyPulse by Michael Lam

PSRFITS specification: www.atnf.csiro.au/people/pulsar/index.html?n=Main.Psrfits

EXAMPLE USAGE:
ar = Archive(filename)
ar.tscrunch()
ar.pavplot()

TODO:

Check POL_TYPE in the above to figure out how to pscrunch

Add emulate_psrchive mode?

Allow for chopping of subints, frequencies, etc?

Flip order of arguments in scrunching?

Check if time() still works.
'''

import numpy as np
import numpy.ma as ma
import gc as g
import matplotlib.pyplot as plt
import time
import py2fits.utils as u
#import pypulse.singlepulse as SP
import py2fits.par as par
Par = par.Par
#import pypulse.calibrator as calib
#Calibrator = calib.Calibrator
import decimal as d
Decimal = d.Decimal
from importlib import import_module
import inspect
import tempfile
import os
try:
    import astropy.io.fits as pyfits
except:
    import pyfits
import astropy.coordinates as coordinates
import astropy.units as units
import sys
if sys.version_info.major == 2:
    fmap = map
elif sys.version_info.major == 3:
    fmap = lambda x,*args: list(map(x,*args))
    xrange = range

PSR = "PSR"
CAL = "CAL"
FON = "FON"
FOF = "FOF"
PCM = "PCM"
SEARCH = "SEARCH"





class Archive:
<<<<<<< HEAD
    def __init__(self,filename,prepare=False,lowmem=False,verbose=True,weight=False,center_pulse=False,baseline_removal=False,wcfreq=False,thread=False,cuda=False):
        ## Parse filename here
        self.pypulse_history = []
        #self.record(inspect.currentframe())
=======
    def __init__(self,filename,prepare=True,lowmem=False,verbose=True,weight=True,center_pulse=True,baseline_removal=True,wcfreq=True,thread=False,cuda=False):
        ## Parse filename here
        self.pypulse_history = []
        self.record(inspect.currentframe())
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        self.filename = str(filename) #fix unicode issue
        self.prepare = prepare
        self.lowmem = lowmem
        self.verbose = verbose
        self.center_pulse = center_pulse
        self.baseline_removal = baseline_removal
        self.wcfreq = wcfreq
        self.thread = thread
        self.cuda = cuda

        if verbose:
            print("Loading: %s" % self.filename)
            t0=time.time()

        self.load(self.filename,prepare=prepare,center_pulse=center_pulse,baseline_removal=baseline_removal,weight=weight,wcfreq=wcfreq)
        if not self.lowmem:
            self.data_orig = np.copy(self.data)
            self.weights_orig = np.copy(self.weights)
        if verbose:
            t1=time.time()
            print("Load time: %0.2f s" % (t1-t0))

        #self.reset(False) #put prepare into here?, copying to arch is done here

        #if prepare:
        #    self.pscrunch()

        #if verbose and prepare:
        #    t2=time.time()
        #    print("Prep time: %0.2f s" % (t2-t1))

    def __repr__(self):
        return "Archive(%r,prepare=%r,lowmem=%r,verbose=%r)" % (self.filename,self.prepare,self.lowmem,self.verbose)
    def __str__(self):
        return self.filename



<<<<<<< HEAD
    def load(self,filename,prepare=False,center_pulse=False,baseline_removal=False,weight=False,wcfreq=False):
        '''
        Loads a PSRFITS file
        http://www.atnf.csiro.au/people/pulsar/index.html?n=PsrfitsDocumentation.Txt

        Parameters:
        ==========

        filename: pathname of PSRFITS file

        prepare: flag
            see reset function: 'Replace the arch 
            with the original clone'


        center_pulse: flag
            not sure what this does
            possibly used elsewhere in py2fits

        baseline_removal: flag
            not sure what this does
            possibly used elsewhere in py2fits


        weight: possibly a flag?
            related to the weights of data 
            possibly read in from the PSRFITS template file
            
        wcfreq: flag
            not sure what this does
            possibly used elsewhere in py2fits

        '''
=======
    def load(self,filename,prepare=True,center_pulse=True,baseline_removal=True,weight=True,wcfreq=False):
        """
        Loads a PSRFITS file and processes
        http://www.atnf.csiro.au/people/pulsar/index.html?n=PsrfitsDocumentation.Txt
        """
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        if filename is None: #Needed?
            filename = self.filename
        try:
            if self.lowmem:
<<<<<<< HEAD
                hdulist = pyfits.open(filename,ignore_missing_end=True,memmap=True) #Probably don't need...
=======
                hdulist = pyfits.open(filename,ignore_missing_end=True,memmap=True)
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
            else:
                hdulist = pyfits.open(filename,ignore_missing_end=True)
        except IOError:
            print("Filename not found")
            raise SystemExit
        self.header = hdulist[0].header
        self.keys = fmap(lambda x: x.name,hdulist)
        tablenames = self.keys[:] #temporary list for checking other tables

<<<<<<< HEAD

        if 'HISTORY' in self.keys:
            tablenames.remove('HISTORY')
            self.history = History(hdulist['HISTORY'])
            self.nsubint = self.history.getLatest("NSUB")
            self.npol = self.history.getLatest("NPOL")
            self.nchan = self.history.getLatest("NCHAN")
            self.nbin = self.history.getLatest("NBIN")
            self.nsblk = self.history.getLatest("NSBLK")
        else:
            self.history = None
            self.nsubint = hdulist['SUBINT'].header['NAXIS2']
            self.nbin,self.nchan,self.npol,self.nsblk = fmap(int,hdulist['SUBINT'].columns[-1].dim[1:-1].split(","))

        nsubint = self.nsubint
        npol = self.npol
        nchan = self.nchan
        nbin = self.nbin
        nsblk = self.nsblk
        observer = self.header[9]

=======
        if 'HISTORY' in self.keys:
            tablenames.remove('HISTORY')
            self.history = History(hdulist['HISTORY'])
            nsubint = self.history.getLatest("NSUB")
            npol = self.history.getLatest("NPOL")
            nchan = self.history.getLatest("NCHAN")
            nbin = self.history.getLatest("NBIN")
        else:
            self.history = None
            nsubint = hdulist['SUBINT'].header['NAXIS2']
            nbin,nchan,npol,nsblk = fmap(int,hdulist['SUBINT'].columns[-1].dim[1:-1].split(","))
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9

        if 'PSRPARAM' in self.keys:
            tablenames.remove('PSRPARAM')
            self.paramheaderlist = hdulist['PSRPARAM'].header.keys()
            self.paramheader = dict()
            for key in self.paramheaderlist:
                self.paramheader[key] = hdulist['PSRPARAM'].header[key]
            self.params = Par(fmap(lambda x: x[0],hdulist['PSRPARAM'].data),numwrap=float)
        elif 'PSREPHEM' in self.keys:
            tablenames.remove('PSREPHEM')
            paramkeys = fmap(lambda x: x.name,hdulist['PSREPHEM'].columns)
            paramvals = hdulist['PSREPHEM'].data[0]
            paramstrs = fmap(lambda x,y: "%s %s"%(x,y),paramkeys,paramvals)
            #self.paramheaderlist = hdulist['PSREPHEM'].header.keys()
            #self.paramheaderdict = dict()
            #for key in self.paramheaderlist:
            #    self.paramheaderdict[key] = hdulist['PSREPHEM'].header[key]
            self.params = Par(paramstrs,numwrap=float)
        else:
            self.params = None



        if 'POLYCO' in self.keys:
            tablenames.remove('POLYCO')
            self.polyco = Polyco(hdulist['POLYCO'],MJD=self.getMJD(full=True))
        else:
            self.polyco = None


        tablenames.remove('PRIMARY')
        if 'SUBINT' in tablenames:
            tablenames.remove('SUBINT')

        isFluxcal = False
        if 'FLUX_CAL' in tablenames: #need to account for is this appropriate?
            tablenames.remove('FLUX_CAL')
            isFluxcal = True
        self.tables = list()
        for tablename in tablenames: #remaining table names to store
            self.tables.append(hdulist[tablename].copy())



        #if self.header['OBS_MODE'] == 'PCM':
        if isFluxcal:

            raise SystemExit




<<<<<<< HEAD
=======

>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        self.subintinfo = dict()
        self.subintinfolist = fmap(lambda x: x.name, hdulist['SUBINT'].columns[:-5])
        for i,column in enumerate(hdulist['SUBINT'].columns[:-5]):
            self.subintinfo[column.name] = (column.format,column.unit,hdulist['SUBINT'].data[column.name])
<<<<<<< HEAD
        
        #TODO possibly make self.primaryinfo
        #not sure if needed

        #Creating primary header dictionary
        #in the same way as the secondary header
        #dictionary is created
        self.primaryheader = dict()
        self.primaryheaderlist = hdulist['PRIMARY'].header.keys()
        for i,key in enumerate(hdulist['PRIMARY'].header):
            self.primaryheader[key] = hdulist['PRIMARY'].header[key]

        #Create secondary header dictionary        
=======
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        self.subintheader = dict()
        self.subintheaderlist = hdulist['SUBINT'].header.keys()#for ordering
        for i,key in enumerate(hdulist['SUBINT'].header):
            self.subintheader[key] = hdulist['SUBINT'].header[key]

        DATA = hdulist['SUBINT'].data['DATA']
<<<<<<< HEAD
        #if np.ndim(DATA)==5:
        #    DATA = DATA[:,0,:,:,:] #remove the nsblk column
=======
        if np.ndim(DATA)==5:
            DATA = DATA[:,0,:,:,:] #remove the nsblk column
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        DATA = np.ascontiguousarray(DATA)

        #Definitions in Base/Formats/PSRFITS/ProfileColumn.C
        DAT_FREQ = hdulist['SUBINT'].data['DAT_FREQ']
        DAT_WTS = np.ascontiguousarray(hdulist['SUBINT'].data['DAT_WTS'])
        if not weight:
            DAT_WTS = np.ones(np.shape(DAT_WTS))
        DAT_SCL = np.ascontiguousarray(hdulist['SUBINT'].data['DAT_SCL'])
        DAT_OFFS = np.ascontiguousarray(hdulist['SUBINT'].data['DAT_OFFS'])# + 0.5 #testing
        self.DAT_SCL = DAT_SCL #testing



        #print DAT_WTS,np.max(DAT_WTS),np.min(DAT_WTS)
        if np.size(DAT_WTS) == 1:
            #DAT_WTS[0] = 1.0
            #DAT_WTS[0,0] = 1.0
            self.weights = np.ones((1,1))
        else:
            #DAT_WTS /= np.max(DAT_WTS) #close???
            #DAT_WTS *= 50
            self.weights = DAT_WTS


        if self.lowmem: #Replace the data arrays with memmaps to reduce memory load



            SHAPE = np.shape(DATA)
            tfDATA = tempfile.NamedTemporaryFile()
            fp = np.memmap(tfDATA.name,dtype=np.int16,mode='w+',shape=SHAPE)
            fp[:] = DATA[:]
            del fp
            DATA = np.memmap(tfDATA.name,dtype=np.int16,mode='r',shape=SHAPE)

            '''
            SHAPE = np.shape(DAT_WTS)
            tfDAT_WTS = tempfile.NamedTemporaryFile()
            fp = np.memmap(tfDAT_WTS.name,dtype=np.float32,mode='w+',shape=SHAPE)
            fp[:] = DAT_WTS[:]
            del fp
            DAT_WTS = np.memmap(tfDAT_WTS.name,dtype=np.float32,mode='r',shape=SHAPE)

            SHAPE = np.shape(DAT_SCL)
            tfDAT_SCL = tempfile.NamedTemporaryFile()
            fp = np.memmap(tfDAT_SCL.name,dtype=np.float32,mode='w+',shape=SHAPE)
            fp[:] = DAT_SCL[:]
            del fp
            DAT_SCL = np.memmap(tfDAT_SCL.name,dtype=np.float32,mode='r',shape=SHAPE)

            SHAPE = np.shape(DAT_OFFS)
            tfDAT_OFFS = tempfile.NamedTemporaryFile()
            fp = np.memmap(tfDAT_OFFS.name,dtype=np.float32,mode='w+',shape=SHAPE)
            fp[:] = DAT_OFFS[:]
            del fp
            DAT_OFFS = np.memmap(tfDAT_OFFS.name,dtype=np.float32,mode='r',shape=SHAPE)
            '''

<<<<<<< HEAD
            #tf = tempfile.NamedTemporaryFile()
            self.data = DATA#np.memmap(tf.name,dtype=np.float32,mode='w+',shape=(nsubint,nsblk,npol,nchan,nbin))

        else:
            self.data = DATA#np.zeros((nsubint,nsblk,npol,nchan,nbin))
=======
            tf = tempfile.NamedTemporaryFile()
            self.data = np.memmap(tf.name,dtype=np.float32,mode='w+',shape=(nsubint,npol,nchan,nbin))

        else:
            self.data = np.zeros((nsubint,npol,nchan,nbin))
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9

        #self.data = np.zeros((nsubint,npol,nchan,nbin))
        #data = np.zeros((nsubint,npol,nchan,nbin))

        I = range(nsubint)
        J = range(npol)
        K = range(nchan)


        self.freq = DAT_FREQ

<<<<<<< HEAD
=======



        if nsubint == 1 and npol == 1 and nchan == 1:
            self.data = (DAT_SCL*DATA+DAT_OFFS)#*DAT_WTS
        elif nsubint == 1 and npol == 1:
            for k in K:
                self.data[0,0,k,:] = (DAT_SCL[0,k]*DATA[0,0,k,:]+DAT_OFFS[0,k])#*DAT_WTS[0,k] #dat WTS[0]?
        elif nsubint == 1 and nchan == 1:
            for j in J:
                self.data[0,j,0,:] = (DAT_SCL[0,j]*DATA[0,j,0,:]+DAT_OFFS[0,j])#*DAT_WTS[0]
        elif npol == 1 and nchan == 1:
            for i in I:
                self.data[i,0,0,:] = (DAT_SCL[i,0]*DATA[i,0,0,:]+DAT_OFFS[i,0])#*DAT_WTS[0]
        else: #if nsubint == 1 or npol == 1 or nchan == 1 this works, or all three are not 1, might want to split this up
            t0 = time.time()
            cudasuccess = False
            if self.cuda:
                try:
                    gpuarray = import_module('pycuda.gpuarray')
                    compiler = import_module('pycuda.compiler')
                    driver = import_module('pycuda.driver')
                    autoinit = import_module('pycuda.autoinit')
                    cudasuccess = True
                except ImportError:
                    print("PyCUDA not imported")
                    #                __global__ void combine(float *retval, float *DAT_SCL, float *DATA, float *DAT_OFFS, int nbin)
                mod = compiler.SourceModule("""
                __global__ void combine(int16_t *retval, float *DAT_SCL, int16_t *DATA, float *DAT_OFFS, int nbin, int size)
                {
                    //uint idx = threadIdx.x + threadIdx.y*4;
                    uint xidx = blockDim.x * blockIdx.x + threadIdx.x;
                    uint yidx = blockDim.y * blockIdx.y + threadIdx.y;
                    uint zidx = blockDim.z * blockIdx.z + threadIdx.z;
                    uint idx = xidx + 5*yidx;
                    if (idx < size)
                        retval[idx] = 1;//DATA[idx];

                    //int subidx = idx/nbin;
                    //retval[idx] = ((DATA[idx] & 0x00FF) <<8) | ((DATA[idx]>>8) & 0x00ff); //swap-endian
                    //retval[idx] = (double)(DAT_SCL[subidx]*DATA[idx]+DAT_OFFS[subidx]);
                }
                """)
                combine = mod.get_function("combine")

                maxX = autoinit.device.get_attributes()[driver.device_attribute.MAX_BLOCK_DIM_X]
                maxY = autoinit.device.get_attributes()[driver.device_attribute.MAX_BLOCK_DIM_Y]
                maxZ = autoinit.device.get_attributes()[driver.device_attribute.MAX_BLOCK_DIM_Z]






                #combine(driver.Out(self.data),driver.In(np.ascontiguousarray(DAT_SCL)),driver.In(np.ascontiguousarray(DATA)),driver.In(np.ascontiguousarray(DAT_OFFS)),nbin,block=(4,4,1))
                data = np.zeros(np.shape(self.data),dtype='>i2')

                if nsubint <= maxX:
                    X = int(nsubint)
                    gridX = 1
                else:
                    pass
                # Assume maxZ >= 4 for now
                Z = int(npol)
                gridZ = 1
                if nchan <= maxY:
                    Y = int(nchan)
                    gridY = 1
                else:
                    pass

                print(np.shape(data),np.size(data))
                retval = np.zeros(np.shape(self.data),dtype='>i2')
                print(X,Y,Z)
                combine(driver.Out(retval),driver.In(DAT_SCL),driver.In(data),driver.In(DAT_OFFS),np.int32(nbin),np.int32(np.size(data)),block=(X,1,1),grid=(1,Z,Y))

                raise SystemExit

                #combine(driver.Out(DATA),driver.In(DAT_SCL),driver.In(DATA),driver.In(DAT_OFFS),nbin,block=(4,4,1))

                combine(driver.Out(retval),driver.In(DAT_SCL),driver.In(data),driver.In(DAT_OFFS),nbin,block=(4,4,4))
                #combine(data_gpu,driver.In(DAT_SCL),data_gpu,driver.In(DAT_OFF),nbin,block=(4,4,1))

                #driver.memcpy_dtoh(retval, data_gpu)
                #combine(driver.Out(data),driver.In(DAT_SCL),driver.In(DATA),driver.In(DAT_OFFS),nbin,block=(4,4,1))
                #print "Break1"
                #print retval,np.all(retval==256)
                #print np.shape(retval)
                #print len(np.where(retval==256)[0])
                #for i in retval:
                #    print i,
                #print "break"
                #print data.#dtype,DATA.dtype
                raise SystemExit

            if self.thread and cudasuccess == False:
                def loop_func(i):
                    for j in J:
                        jnchan = j*nchan
                        for k in K:
                            self.data[i,j,k,:] = (DAT_SCL[i,jnchan+k]*DATA[i,j,k,:]+DAT_OFFS[i,jnchan+k])#*DAT_WTS[i,k]
                u.parmap(loop_func,I)
            elif cudasuccess == False:
                for i in I:
                    for j in J:
                        jnchan = j*nchan
                        for k in K:
                            self.data[i,j,k,:] = (DAT_SCL[i,jnchan+k]*DATA[i,j,k,:]+DAT_OFFS[i,jnchan+k])#*DAT_WTS[i,k]
            t1 = time.time()
            #print t1-t0



>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        bw = self.getBandwidth()


        # All time-tagging info
        self.durations = self.getSubintinfo('TSUBINT')
        self.subint_starts = np.array(fmap(Decimal,self.getSubintinfo('OFFS_SUB')),dtype=np.dtype(Decimal))#+self.getTbin(numwrap=Decimal)*Decimal(nbin/2.0)#+self.getMJD(full=False,numwrap=Decimal) #converts center-of-bin times to start-of-bin times, in seconds, does not include the integer MJD part. This means that a template sitting in the center will have zero extra time
        self.channel_delays = np.zeros(nchan,dtype=np.dtype(Decimal)) #used to keep track of frequency-dependent channel delays, in time units.

<<<<<<< HEAD
=======
        if prepare and not self.isCalibrator():
            self.pscrunch()
            self.dedisperse(wcfreq=wcfreq)

        self.calculateAverageProfile()


        if center_pulse and not self.isCalibrator() and prepare: #calibrator is not a pulse, prepare must be run so that dedisperse is run?
            self.center()

        if baseline_removal and not self.isCalibrator():
            self.removeBaseline()

>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        hdulist.close()
        return


    def save(self,filename):
<<<<<<< HEAD
        '''
        Save the file to a new FITS file. 
        Initializes a new primary header and then
        appends the secondary header at the end 
        of the function.

        Parameters:
        ==========
        filename: name for new FITS file

        '''

        #Creating the new primary header to be able to change pieces
        #of the header when needed.
        primaryhdr = pyfits.Header()
        for key in self.primaryheaderlist:
            primaryhdr[key] = self.primaryheader[key]
        primaryhdu = pyfits.PrimaryHDU(header=primaryhdr) #need to make alterations to header
        
        hdulist = pyfits.HDUList(primaryhdu)    #Don't need to append anything to hdulist 
                                                #for the primary header since this initializes 
                                                #the primary header
=======
        """Save the file to a new FITS file"""

        primaryhdu = pyfits.PrimaryHDU(header=self.header) #need to make alterations to header
        hdulist = pyfits.HDUList(primaryhdu)
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9

        if self.history is not None:
            cols = []
            for name in self.history.namelist:
                fmt,unit,array = self.history.dictionary[name]
                #print name,fmt,unit,array
                col = pyfits.Column(name=name,format=fmt,unit=unit,array=array)
                cols.append(col)
            historyhdr = pyfits.Header()
            for key in self.history.headerlist:
                historyhdr[key] = self.history.header[key]
            historyhdu = pyfits.BinTableHDU.from_columns(cols,name='HISTORY',header=historyhdr)
            hdulist.append(historyhdu)
            # Need to add in PyPulse changes into a HISTORY
        #else: #else start a HISTORY table


        if self.params is not None:
            cols = [pyfits.Column(name='PSRPARAM',format='128A',array=self.params.filename)] #PARAM and not PSRPARAM?
            paramhdr = pyfits.Header()
            for key in self.paramheaderlist:
                paramhdr[key] = self.paramheader[key]
            paramhdu = pyfits.BinTableHDU.from_columns(cols,name='PSRPARAM')
            hdulist.append(paramhdu)
            # Need to include mode for PSREPHEM


<<<<<<< HEAD
=======


>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        if self.polyco is not None:
            cols = []
            for name in self.polyco.namelist:
                fmt,unit,array = self.polyco.dictionary[name]
                #print name,fmt,unit,array
                col = pyfits.Column(name=name,format=fmt,unit=unit,array=array)
                cols.append(col)
            polycohdr = pyfits.Header()
            for key in self.polyco.headerlist:
                polycohdr[key] = self.polyco.header[key]
            polycohdu = pyfits.BinTableHDU.from_columns(cols,name='POLYCO',header=polycohdr)
            hdulist.append(polycohdu)


<<<<<<< HEAD
=======


>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        if len(self.tables) > 0:
            for table in self.tables:
                hdulist.append(table)

        cols = []
        for name in self.subintinfolist:
            fmt,unit,array = self.subintinfo[name]
            col = pyfits.Column(name=name,format=fmt,unit=unit,array=array)
            cols.append(col)
            # finish writing out SUBINT!

        cols.append(pyfits.Column(name='DAT_FREQ',format='%iE'%np.shape(self.freq)[1],unit='MHz',array=self.freq)) #correct size? check units?
        cols.append(pyfits.Column(name='DAT_WTS',format='%iE'%np.shape(self.weights)[1],array=self.weights)) #call getWeights()

<<<<<<< HEAD
        nsubint,nsblk,npol,nchan,nbin = self.data.shape

        DAT_OFFS = np.zeros((nsubint,npol*nchan),dtype=np.float32)
        DAT_SCL = np.zeros((nsubint,npol*nchan),dtype=np.float32)
        DATA = self.data #getData(squeeze=False,weight=False)
        saveDATA = self.data #np.zeros(self.shape(squeeze=False),dtype=np.int16)
        # Following Base/Formats/PSRFITS/unload_DigitiserCounts.C
        #for i in xrange(nsubint):
        #    for j in xrange(npol):
        #        jnchan = j*nchan
        #        for k in xrange(nchan):
        #            MIN = np.min(DATA[i,j,k,:])
        #            MAX = np.max(DATA[i,j,k,:])
        #            RANGE = MAX - MIN
        #            if MAX == 0 and MIN == 0:
        #                DAT_SCL[i,jnchan+k] = 1.0
        #            else:
        #                DAT_OFFS[i,jnchan+k] = 0.5*(MIN+MAX)
        #                DAT_SCL[i,jnchan+k] = (MAX-MIN)/32766.0 #this is slightly off the original value? Results in slight change of data
        #
        #            saveDATA[i,j,k,:] = np.floor((DATA[i,j,k,:] - DAT_OFFS[i,jnchan+k])/DAT_SCL[i,jnchan+k] + 0.5) #why +0.5?

        cols.append(pyfits.Column(name='DAT_OFFS',format='%iE'%np.size(DAT_OFFS[0]),array=DAT_OFFS))
        cols.append(pyfits.Column(name='DAT_SCL',format='%iE'%np.size(DAT_SCL[0]),array=DAT_SCL))
        cols.append(pyfits.Column(name='DATA',format='%iI'%np.size(saveDATA[0]),array=saveDATA,unit='Jy',dim='(%s,%s,%s,%s,%s)'%(nbin,nsblk,npol,nchan,nbin))) #replace the unit here

        #Basically does the same thing as the primary
        #header initialization at the beginning of the function
=======
        nsubint = self.getNsubint()
        npol = self.getNpol()
        nchan = self.getNchan()
        nbin = self.getNbin()
        DAT_OFFS = np.zeros((nsubint,npol*nchan),dtype=np.float32)
        DAT_SCL = np.zeros((nsubint,npol*nchan),dtype=np.float32)
        DATA = self.getData(squeeze=False,weight=False)
        saveDATA = np.zeros(self.shape(squeeze=False),dtype=np.int16)
        # Following Base/Formats/PSRFITS/unload_DigitiserCounts.C
        for i in xrange(nsubint):
            for j in xrange(npol):
                jnchan = j*nchan
                for k in xrange(nchan):
                    MIN = np.min(DATA[i,j,k,:])
                    MAX = np.max(DATA[i,j,k,:])
                    RANGE = MAX - MIN
                    if MAX == 0 and MIN == 0:
                        DAT_SCL[i,jnchan+k] = 1.0
                    else:
                        DAT_OFFS[i,jnchan+k] = 0.5*(MIN+MAX)
                        DAT_SCL[i,jnchan+k] = (MAX-MIN)/32766.0 #this is slightly off the original value? Results in slight change of data

                    saveDATA[i,j,k,:] = np.floor((DATA[i,j,k,:] - DAT_OFFS[i,jnchan+k])/DAT_SCL[i,jnchan+k] + 0.5) #why +0.5?

        cols.append(pyfits.Column(name='DAT_OFFS',format='%iE'%np.size(DAT_OFFS[0]),array=DAT_OFFS))
        cols.append(pyfits.Column(name='DAT_SCL',format='%iE'%np.size(DAT_SCL[0]),array=DAT_SCL))
        cols.append(pyfits.Column(name='DATA',format='%iI'%np.size(saveDATA[0]),array=saveDATA,unit='Jy',dim='(%s,%s,%s)'%(nbin,nchan,npol))) #replace the unit here

>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        subinthdr = pyfits.Header()
        for key in self.subintheaderlist:
            subinthdr[key] = self.subintheader[key]
        subinthdu = pyfits.BinTableHDU.from_columns(cols,name='SUBINT',header=subinthdr)
<<<<<<< HEAD

        #Appends the secondary header
        hdulist.append(subinthdu)


        #Writes the data
        hdulist.writeto(filename,overwrite=True)#clobber=True?
=======
        hdulist.append(subinthdu)



        hdulist.writeto(filename,clobber=True)#clobber=True?
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9



    def unload(self,filename):
        return self.save(filename)

    def gc(self):
        """Manually clear the data cube for python garbage collection"""
        if self.verbose:
            t0=time.time()
        self.data = None
        self.data_orig = None
        self.weights = None
        self.weights_orig = None
        if self.verbose:
            t1=time.time()
            print("Unload time: %0.2f s" % (t1-t0))
        g.collect()


    def shape(self,squeeze=True):
        """Return the current shape of the data array"""
        return np.shape(self.getData(squeeze=squeeze))
    def reset(self,prepare=True):
        """Replace the arch with the original clone"""
<<<<<<< HEAD
=======
        self.record(inspect.currentframe()) #temporary, actually change the history!
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        if self.lowmem:
            self.load(self.filename,prepare=prepare)
        else:
            self.data = np.copy(self.data_orig)
            self.weights = np.copy(self.weights_orig)
        self.durations = self.getSubintinfo('TSUBINT')
        #if prepare:
        #    self.scrunch()

<<<<<<< HEAD

=======
    def getLevels(self,differences=False):
        """ Returns calibration levels if this is a calibrator"""
        if not self.isCalibrator():
            print("Not a calibration file")
            return
        # Ensure data are scrunched in time, or ignore this and simply calculate the weighted time average?
        self.tscrunch()

        #Pre-define
        data = self.getData()
        npol = self.getNpol()
        nchan = self.getNchan()
        nbin = self.getNbin()

        # Check header info CAL_DCYC, CAL_NPHS, etc, to determine on-diode
        # or take an absolute value?
        first = np.mean(data[0,:,:nbin/2])
        second = np.mean(data[0,:,nbin/2:])
        if first > second:
            highinds = np.arange(0,nbin/2)
            lowinds = np.arange(nbin/2,nbin)
        else:
            lowinds = np.arange(0,nbin/2)
            highinds = np.arange(nbin/2,nbin)

        # Calculate calibrations
        freqs = self.getAxis('F')
        if differences:
            caldata = np.zeros((npol,nchan))
            calerrs = np.zeros((npol,nchan))
            for i in xrange(npol):
                for j in xrange(nchan):
                    caldata[i,j] = np.mean(data[i,j,highinds]) - np.mean(data[i,j,lowinds])
                    calerrs[i,j] = np.sqrt(np.std(data[i,j,highinds])**2 / len(highinds) + np.std(data[i,j,lowinds])**2 / len(lowinds))
        else:
            caldatalow = np.zeros((npol,nchan))
            caldatahigh = np.zeros((npol,nchan))
            calerrslow = np.zeros((npol,nchan))
            calerrshigh = np.zeros((npol,nchan))
            for i in xrange(npol):
                for j in xrange(nchan):
                    caldatalow[i,j] = np.mean(data[i,j,lowinds])
                    caldatahigh[i,j] = np.mean(data[i,j,highinds])
                    calerrslow[i,j] = np.std(data[i,j,lowinds]) / np.sqrt(len(lowinds))
                    calerrshigh[i,j] = np.std(data[i,j,highinds]) / np.sqrt(len(highinds))



        if differences:
            return freqs,caldata,calerrs
        return freqs,caldatalow,caldatahigh,calerrslow,calerrshigh

    def calibrate(self,psrcalar,fluxcalonar=None,fluxcaloffar=None):
        """Calibrates using another archive"""
        self.record(inspect.currentframe())
        if not (isinstance(psrcalar,Calibrator) or (isinstance(psrcalar,Archive) and psrcalar.isCalibrator())):
            raise ValueError("Require calibration archive")
        # Check if cals are appropriate?


        if isinstance(psrcalar,Calibrator):
            cal = psrcalar
        else:
            # Calculate calibration levels
            psrcalfreqs,psrcaldata,psrcalerrs = psrcalar.getLevels(differences=True)

            # Check if cal has the correct dimensions, if not perform interpolation
            freqs = self.getAxis('F')
            if len(freqs) != len(psrcalfreqs):
                pass

            cal = Calibrator(psrcalfreqs,psrcaldata,psrcalerrs)
        #if fluxcalon is not None:
        #    fluxcaloncal = Calibrator(fluxcalonfreqs,fluxcalondata,fluxcalonerrs)
        #    fluxcaloffcal = Calibrator(fluxcalofffreqs,fluxcaloffdata,fluxcalofferrs)



        if fluxcalonar is not None and fluxcaloffar is not None:
            fluxcaldata = np.zeros((npol,nchan))
            for i in xrange(npol):
                for j in xrange(nchan):
                    fluxcaldata[i,j] = np.mean(fdata[i,j,highinds]) - np.mean(fdata[i,j,lowinds])

        cal.applyCalibration(self)
        return cal
        # Apply calibrations









    def getData(self,squeeze=True,setnan=None,weight=True):
        """Returns the data array, fully squeezed"""
        if weight:
            data = np.zeros_like(self.data)
            I,J,K,L = np.shape(self.data)
            I = range(I)
            J = range(J)
            K = range(K)
            for i in I:
                #print np.shape(self.weights),np.shape(data),np.shape(self.data),i#,j,k
                for j in J:
                    for k in K:
                        data[i,j,k,:] = self.data[i,j,k,:]*self.weights[i,k]
        else:
            data = self.data

        if squeeze:
            data = data.squeeze()

        if setnan is not None:
            data = np.where(data==setnan,np.nan,data)

        return np.copy(data) #removes pointer to data
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
    def setData(self,newdata):
        """Sets the data, very dangerous!"""
        self.record(inspect.currentframe())
        if np.shape(newdata) == np.shape(self.data):
            self.data = np.copy(newdata)
<<<<<<< HEAD
    def getWeights(self,squeeze=False):
=======
    def getWeights(self,squeeze=True):
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
        """ Return copy of weights array """
        weights = self.weights
        if squeeze:
            weights = weights.squeeze()
        return np.copy(weights)
    def setWeights(self,val,t=None,f=None):
        """
        Set weights to a certain value
        Can be used for RFI routines
        """
        self.record(inspect.currentframe())
        if t is None and f is None:
            self.weights[:,:] = val
        elif t is None:
            self.weights[:,f] = val
        elif f is None:
            self.weights[t,f] = val
        else:
            self.weights[t,f] = val



    def saveData(self,filename=None,ext='npy',ascii=False):
        """Save the data array to a different format"""
        if filename is None:
            filename = self.filename
            filename = ".".join(filename.split(".")[:-1])+"."+ext
        if self.verbose:
            print("Saving: %s" % filename)
        if ascii:
            shape = self.shape(squeeze=False)
            nsubint = self.getNsubint()
            npol = self.getNpol()
            nchan = self.getNchan()
            nbin = self.getNbin()
            output = ""
            if shape[0] == 1 and shape[1] == 1 and shape[2] == 1:
                np.savetxt(filename,self.getData())
                return
            elif ((shape[0] == 1 and shape[1] == 1) or
                  (shape[0] == 1 and shape[2] == 1) or
                  (shape[1] == 1 and shape[2] == 1)):
                np.savetxt(filename,self.getData())
                return
            for i in xrange(nsubint):
                for j in xrange(npol):
                    for k in xrange(nchan):
                        for l in xrange(nbin):
                            output += "%i %i %i %i %.18e\n" % (i,j,k,l,self.data[i,j,k,l])

            FILE = open(filename,'w')
            FILE.write(output)
            FILE.close()
        else:
            np.save(filename,self.getData())
        return



    def outputPulses(self,filename):
        """ Write out a .npy file"""
        np.save(filename,self.getData())
        return



    ### ==============================
    ### Simple get functions
    ### ==============================


    def getNsubint(self):
        """Returns number of subintegrations"""
        return self.shape(squeeze=False)[0]
<<<<<<< HEAD
    def getNsblk(self):
        """Returns number of NSBLK"""
        return self.data.shape[1]
=======
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
    def getNpol(self):
        """Returns number of polarizations"""
        return self.shape(squeeze=False)[1]
    def getNchan(self):
        """Returns number of channels"""
        return self.shape(squeeze=False)[2]
    def getNbin(self):
        """Returns number of phase bins"""
        return self.shape(squeeze=False)[3]
<<<<<<< HEAD
=======
    def getPeriod(self,header=False):
        """Returns period of the pulsar"""
        if self.isCalibrator():
            return 1.0/self.header['CAL_FREQ']
        if self.params is None:
            return None
        if header or self.polyco is None:
            return self.params.getPeriod()
        else:
            P0 = self.polyco.calculatePeriod()
            #print P0,self.params.getPeriod()
            if np.abs(P0)<1e-5: #Problem with large DT POLYCO values?
                return self.params.getPeriod()
            else:
                ratio = (P0-self.params.getPeriod())/self.params.getPeriod()
                if ratio < 0.5 or ratio > 2:
                    return self.params.getPeriod()
            return P0
            return self.polyco.calculatePeriod()
    # Best replacement for without PSRCHIVE
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
    def getValue(self,value):
        """Looks for a key in one of the headers and returns"""
        if value in self.header.keys():
            return self.header[value]
        if value in self.subintinfo.keys():
            return self.subintinfo[value][-1]
        if self.params is None:
            return None
        return self.params.get(value) #will return None if non-existent
    def getSubintinfo(self,value):
        """Returns value from subintinfo"""
        if value in self.subintinfo.keys():
            return self.subintinfo[value][-1]
        return None
    def getName(self):
        """Returns pulsar name"""
        return self.header['SRC_NAME']
    def getMJD(self,full=False,numwrap=float):
        """Returns MJD of observation"""
        if full:
            return numwrap(self.header['STT_IMJD'])+(numwrap(self.header['STT_SMJD'])+numwrap(self.header['STT_OFFS']))/numwrap(86400)
        return numwrap(self.header['STT_IMJD'])+numwrap(self.header['STT_OFFS'])
    def getTbin(self,numwrap=float):
        """Returns the time per bin"""
        return numwrap(self.getPeriod()) / numwrap(self.getNbin())
    def getDM(self):
        """Returns the data header DM"""
        return self.subintheader['DM']
        #if self.params is None:
        #    return
        #return self.params.getDM()
    def getRM(self):
        """Returns the data header RM"""
        return self.subintheader['RM']
    def getCoords(self,string=False,parse=False):
        """Returns the coordinate info in the header"""
        if string:
            RA = self.header['RA']
            dec = self.header['DEC']
            return RA,dec
        elif parse:
            RA = tuple(map(float,self.header['RA'].split(":")))
            dec = tuple(map(float,self.header['DEC'].split(":")))
            return RA,dec
        return coordinates.SkyCoord("%s %s"%(self.header['RA'],self.header['DEC']),unit=(units.hourangle,units.degree))
    getPulsarCoords = getCoords

    def getTelescopeCoords(self):
        """Returns the telescope coordinates"""
        return self.header['ANT_X'],self.header['ANT_Y'],self.header['ANT_Z']
    def getBandwidth(self,header=False):
        """Returns the observation bandwidth"""
        if header:
            return self.header['OBSBW']
        else:
            return self.subintheader['CHAN_BW']*self.subintheader['NCHAN']
    def getDuration(self):
        """Returns the observation duratrion"""
        #return np.sum(self.subintinfo['TSUBINT']) #This is constant.
        return np.sum(self.getSubintinfo('TSUBINT')) #This is constant.
    def getDurations(self):
        """Returns the subintegration durations"""
        return self.durations
    def getCenterFrequency(self,weighted=False):
        """Returns the center frequency"""
        if weighted:
            return np.sum(self.freq*self.weights)/np.sum(self.weights)
        if "HISTORY" in self.keys:
            return self.history.getLatest('CTR_FREQ')
        else:
            return self.header['OBSFREQ'] #perhaps do an unweighted version from DAT_FREQ?
    def getTelescope(self):
        """Returns the telescope name"""
        return self.header['TELESCOP']
    def getFrontend(self):
        """Returns the frontend name"""
        return self.header['FRONTEND']
    def getBackend(self):
        """Returns the backend name"""
        return self.header['BACKEND']
    def getSN(self):
        """Returns the average pulse S/N"""
        return self.spavg.getSN()

    def isCalibrator(self):
        if self.header['OBS_MODE'] == CAL or self.header['OBS_MODE'] == FON or self.header['OBS_MODE'] == FOF:
            return True
        return False


    #def varargtest(self,*args):
    #    self.record(inspect.currentframe())
    def record(self,frame):
        args, varargs, keywords, values = inspect.getargvalues(frame)
        funcname = frame.f_code.co_name
        string = "%s("%funcname
        for arg in args[1:]:
            if type(values[arg]) == str:
                string += "%s=\"%s\","%(arg,values[arg])
            else:
                string += "%s=%s,"%(arg,values[arg])
        if varargs is not None:  # Var args typically not implemented in PyPulse
            argdict = values[varargs]
            string += "%s,"%(str(argdict)[1:-1].replace(", ",","))
        if keywords is not None:
            kwargdict = values[keywords]
            for kwarg in kwargdict:
                string += "%s=%s,"%(kwarg,kwargdict[kwarg])
        if string[-1] == "(":
            string += ")"
        else:
            string = string[:-1] + ")"
        self.pypulse_history.append(string)

    def print_pypulse_history(self):
        for elem in self.pypulse_history:
            print(elem)






# Takes hdulist['HISTORY']
class History:
    def __init__(self,history):
        """Intializer"""
        self.header = dict()
        self.headerlist = history.header.keys()
        for key in self.headerlist:
            self.header[key] = history.header[key]
        self.dictionary = dict()
        self.namelist = list()
        for col in history.columns:
            self.namelist.append(col.name)
            self.dictionary[col.name] = (col.format,col.unit,list(col.array)) #make a np.array?
    def getValue(self,field,num=None):
        """Returns a dictionary array value for a given numeric entry"""
        if num is None:
            return self.dictionary[field][-1]
        else:
            try:
                return self.dictionary[field][-1][num]
            except IndexError:
                print("Entry out of range")
                return None
    def getLatest(self,field):
        """Returns the latest key value"""
        return self.getValue(field,-1)
    def printEntry(self,i):
        """Prints the i-th history entry"""
        for name in self.namelist:
            value = self.getValue(name,i)
            if value is None:
                return
            print(name,self.getValue(name,i))

# Takes hdulist['POLYCO']
# Similar to History class
class Polyco:
    def __init__(self,polyco,MJD=None):
        """Initializer"""
        self.MJD = MJD
        self.header = dict()
        self.headerlist = polyco.header.keys()
        for key in self.headerlist:
            self.header[key] = polyco.header[key]
        self.dictionary = dict()
        self.namelist = list()
        for col in polyco.columns:
            self.namelist.append(col.name)
            self.dictionary[col.name] = (col.format,col.unit,list(col.array)) #make a np.array?
    def getValue(self,field,num=None):
        """Returns a dictionary array value for a given numeric entry"""
        if num is None:
            return self.dictionary[field][-1]
        else:
            return self.dictionary[field][-1][num]
    def getLatest(self,field):
        """Returns the latest key value"""
        return self.getValue(field,-1)
    def calculate(self,MJD=None):
        if self.MJD is None and MJD is None:
            pass
        elif MJD is None:
            MJD = self.MJD

        #NSITE = self.getValue('NSITE',num=0)
        REF_FREQ = self.getValue('REF_FREQ',num=0)
        #PRED_PHS = self.getValue('PRED_PHS',num=0)
        REF_MJD = self.getValue('REF_MJD',num=0)
        REF_PHS = self.getValue('REF_PHS',num=0)
        REF_F0 = self.getValue('REF_F0',num=0)
        COEFF = self.getValue('COEFF',num=0)
        #print "POLYCO",REF_FREQ,REF_MJD,REF_PHS,REF_F0,COEFF
        #http://tempo.sourceforge.net/ref_man_sections/tz-polyco.txt
        DT = (MJD-REF_MJD)*1440.0
        #print "DT",DT,MJD,REF_MJD
        PHASE = REF_PHS + DT*60*REF_F0
        #print "PHASE",PHASE
        FREQ = 0.0
        for i,c in enumerate(COEFF):
            PHASE += c*np.power(DT,i)
            if i==0:
                continue
            FREQ += c*i*np.power(DT,i-1)
        FREQ = REF_F0 + FREQ/60.0
        return PHASE,FREQ
    def calculatePeriod(self,MJD=None):
        PHASE,FREQ = self.calculate(MJD=MJD)
<<<<<<< HEAD
        return 1.0/FREQ
=======
        return 1.0/FREQ
>>>>>>> 24e4811192ebccb5168b7ecf9be88ce476b769c9
