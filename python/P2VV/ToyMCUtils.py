import optparse
import sys
import os

def __check_req_kw__( name, kwargs ) :
    if not name in kwargs : raise KeyError( 'Must provide kw argument %s' % name )

class Toy(object):
    def __init__(self):
        self._parser = optparse.OptionParser(usage = '%prog')

        self._parser.add_option("-o", "--output", dest = "output", default = 'toy.root',
                                type = 'string', help = "set output filename")
        self._parser.add_option("-n", "--ntoys", dest = "ntoys", default = 100,
                                type = 'int', help = 'number of toys to run')
        self._parser.add_option("--ncpu", dest = "ncpu", default = 1,
                                type = 'int', help = 'number of CPUs to use')
        self._parser.add_option("-e", "--nevents", dest = "nevents", default = 10000,
                                type = 'int', help = 'number of events to generate')
        self._parser.add_option("-s", "--snapshot", dest = "snapshot", default = '',
                                action = 'store', help = 'Extract a snapshot to current directory.')
        self._parser.add_option("--protodata", dest = "protodata", default = '',
                                action = 'store', help = 'use protodata')

        self.__fit_opts = dict(Save = True, Optimize = 1, Verbose = True, Minos = (),
                               Minimizer = 'Minuit2')
        self.__transform = None
        
    def parser(self):
        return self._parser

    def configure(self):
        (self._options, self._args) = self._parser.parse_args()

        if not self._options.snapshot:
            return (self._options, self._args)

        
        import tarfile
        archive = None
        try:
            archive = tarfile.open(self._options.snapshot, 'r:bz2')
            python_dirs = []
            for member in archive.getmembers():
                if member.isfile() and os.path.exists(member.path):
                    print "File %s already exists, skipping" % member.path
                else:
                    archive.extract(member)
                if member.isdir() and member.path.endswith('python'):
                    python_dirs.append(member.path)
            sys.path.extend([os.path.join(os.path.realpath('.'), d) for d in python_dirs])
        except OSError, e:
            print e
            sys.exit(-2)
        finally:
            if archive: archive.close()
        try:
            print 'Running with snapshot %(comment)s' % archive.pax_headers
        except KeyError:
            pass

        from ROOT import gROOT
        gROOT.SetBatch(True)

        return (self._options, self._args)

    def run(self, **kwargs):
        pass

    def data(self):
        return self._data
        
    def write_output(self):
        # Write the results to a file
        from ROOT import TFile
        output_file = TFile.Open(self.options().output, 'recreate')
        output_file.WriteTObject(self._data, self._data.GetName())
        gp = self.gen_params()
        if gp:
            output_file.WriteTObject(gp, 'gen_params')
        output_file.Close()
        
    def set_fit_opts(self, **opts):
        self.__fit_opts.update(opts)

    def set_transform(self, t):
        self.__transform = t

    def transfrom(self):
        return self.__transform

    def gen_params(self):
        return None

    def fit_opts(self):
        return self.__fit_opts

    def options(self):
        return self._options
    
def FitToy(Toy):
    def __init__(self, *args, **kwargs):
        Toy.__init__(self, *args, **kwargs)

    def run(self, **kwargs):
        from ROOT import RooArgSet

        __check_req_kw__('Observables', kwargs)
        __check_req_kw__('Pdf', kwargs)

        observables = kwargs.pop('Observables')
        obs_set = RooArgSet(*observables)
        
        pdf = kwargs.pop('Pdf')
        genPdf = kwargs.pop('GenPdf', pdf)

        gen_obs_set = RooArgSet()
        for o in list(observables) + list(genPdf.ConditionalObservables()):
            gen_obs_set.add(o._target_())
        gen_pdf_params = genPdf.getParameters(gen_obs_set).snapshot(True)

        genPdf = genPdf.clone(genPdf.GetName() + "_toy_clone")
        genPdf.recursiveRedirectServers(gen_pdf_params)

        fit_obs_set = RooArgSet()
        for o in list(observables) + list(pdf.ConditionalObservables()):
            fit_obs_set.add(o._target_())
        params = pdf.getParameters(fit_obs_set)

        pdf_params = RooArgSet()
        for p in params:
            if p.isConstant(): continue
            pdf_params.add(p)
        ## for param in pdf_params:
        ##     if param.GetName() not in ['Gamma', 'dGamma']:
        ##         param.setConstant()
        self._gen_params = pdf_params.snapshot(True)

        # Make another ArgSet to put the fit results in
        result_params = RooArgSet(pdf_params, "result_params")

        transfrom = self.transform():
        if transfrom:
            trans_params = transform.gen_params(gen_obs_set)
            result_params.add(trans_params)
            
        # Some extra numbers of interest
        from ROOT import RooRealVar
        NLL = RooRealVar('NLL', '-log(Likelihood)', 1.)
        ngen = RooRealVar('ngen', 'number of generated events', self.options().nevents)
        seed = RooRealVar('seed', 'random seed', 0.)
        from ROOT import RooCategory
        status = RooCategory('status', 'fit status')
        status.defineType('success', 0)
        status.defineType('one', 1)
        status.defineType('two', 2)
        status.defineType('three', 3)
        result_params.add(status)
        result_params.add(NLL)
        result_params.add(ngen)
        result_params.add(seed)

        # The dataset to store the results
        from ROOT import RooDataSet
        self._data = RooDataSet('result_data', 'result_data', result_params)
        data_params = self._data.get()

        from ROOT import RooRandom
        import struct, os

        i = 0
        while i < self.options().ntoys:
            # Get a good random seed, set it and store it
            s = struct.unpack('I', os.urandom(4))[0]    
            RooRandom.randomGenerator().SetSeed(s)
            seed.setVal(s)

            # Reset pdf parameters to initial values. Note: this does not reset the estimated errors...
            pdf_params.assignValueOnly(self.gen_params()) 
            args = dict(NumEvents = self.options().nevents)
            if 'ProtoData' in kwargs:
                args['ProtoData'] = kwargs.pop('ProtoData')
            
            genPdf.getParameters(obs_set).assignValueOnly(gen_pdf_params)
            data = genPdf.generate(obs_set, **args)
            if transform:
                data = transform(data)
                if not data:
                    # Transform has failed
                    transform.set_params(data_params)
                    self._data.fill()
                    continue
            
            if data.isWeighted() and 'SumW2Error' not in self.fit_opts():
                self.fit_opts()['SumW2Error'] = False

            while j < 4: 
                fit_result = pdf.fitTo(data, NumCPU = self.options().ncpu, **(self.fit_opts()))
                if fit_result.status() == 0:
                    break
                j += 1
            if fit_result.status() != 0:
                print 'Fit result status = %s' % fit_result.status()
            NLL.setVal(fit_result.minNll())
            status.setIndex(fit_result.status())
            for result_param in result_params:
                data_param = data_params.find(result_param.GetName())
                if isinstance(result_param, RooCategory):
                    data_param.setIndex(result_param.getIndex())
                else:
                    data_param.setVal(result_param.getVal())
                    # This sets a symmetric error, but since we don't run Minos, that's ok
                    data_param.setError(result_param.getError())
            if transform:
                transform.set_params(self, data_params)
                
            self._data.fill()

        return self.data()

    def gen_params(self):
        return self._gen_params

def DilutionToy(Toy):
    def __init__(self, *args, **kwargs):
        Toy.__init__(self, *args, **kwargs)

    def run(self, **kwargs):
        from ROOT import RooArgSet

        __check_req_kw__('Observables', kwargs)
        __check_req_kw__('Pdf', kwargs)
        __check_req_kw__('Sigmat', kwargs)
        __check_req_kw__('SigmatCat', kwargs)

        observables = kwargs.pop('Observables')
        obs_set = RooArgSet(*observables)
        
        pdf = kwargs.pop('Pdf')
        sigmat_cat = kwargs.pop('SigmatCat')
        sigmat = kwargs.pop('Sigmat')
        
        gen_obs_set = RooArgSet(*observables).snapshot(True)

        # Make another ArgSet to put the fit results in
        result_params = RooArgSet("result_params")

        da = RealVar('da', Observable = True, MinMax = (0.01, 1.1))
        dft = RealVar('dft', Observable = True, MinMax = (0.01, 1.1))
        result_params.add(da._target_())
        result_params.add(dft._target_())
        
        transfrom = self.transform():
        if transfrom:
            trans_params = transform.gen_params(gen_obs_set)
            result_params.add(trans_params)
            
        # Some extra numbers of interest
        from ROOT import RooRealVar
        seed = RooRealVar('seed', 'random seed', 0.)
        result_params.add(seed)
        
        # The dataset to store the results
        from ROOT import RooDataSet
        self._data = RooDataSet('result_data', 'result_data', result_params)
        data_params = self._data.get()

        from ROOT import RooRandom
        import struct, os

        while self._data.numEntries() < self.options().ntoys:
            # Get a good random seed, set it and store it
            s = struct.unpack('I', os.urandom(4))[0]    
            RooRandom.randomGenerator().SetSeed(s)
            seed.setVal(s)

            # Reset pdf parameters to initial values. Note: this does not reset the estimated errors...
            args = dict(NumEvents = self.options().nevents)
            if 'ProtoData' in kwargs:
                args['ProtoData'] = kwargs.pop('ProtoData')
            
            data = pdf.generate(obs_set, **args)
            if self.transform():
                data = self.transform()(data)
                if not data:
                    transform.set_params(data_params)
                    self._data.fill()
                    continue

            st_cat = data.addColumn(sigmat_cat._target_())
            from P2VV import Dilution
            d_ft = Dilution.dilution_bins(sdata, self.__t, sigmat, st_cat, t_range = 2)
            d_a = Dilution.signal_dilution_dg(sdata, sigmat, 1.2, 0.2, 2)
            da.setVal(d_a[0])
            da.setError(d_a[1])
            dft.setVal(d_ft[0])
            dft.setError(d_ft[1])
        
            if transform:
                transform.set_params(data_params)
            
            self._data.fill()

        return self.data()

    def gen_params(self):
        return self._gen_params

class SWeightTransform(object):
    def __init__(self, pdf, component, fit_opts):
        self.__comp = component
        self.__pdf = pdf
        self.__fit_opts = fit_opts
        self.__result = None

        self.__status = RooCategory('sweight_status', 'sweight fit status')
        self.__status.defineType('success', 0)
        self.__status.defineType('one', 1)
        self.__status.defineType('two', 2)
        self.__status.defineType('three', 3)
        
    def __call__(self, data):
        pdf_pars = self.__pdf.getParameters(data.get())
        if not hasattr(self, '__parameters'):
            self.__parameters = pdf_pars.snapshot(True)
            self.__parameters.add(self.__status)
        else:
            for p in self.__parameters:
                pdf_par = pdf_pars.find(p.GetName())
                if not pdf_par:
                    continue
                pdf_par.setVal(p.getVal())
                pdf_par.setError(p.getError())

        success = False
        for i in range(3):
            self.__result = self.__pdf.fitTo(data, **self.__fit_opts)
            if self.__result.status() == 0:
                success = True
                break

        self.__status.setIndex(self.__result.status())
        if success:
            from P2VV.Utilities.SWeights import SData
            sData = SData(Pdf = self.__pdf, Data = data, Name = 'MassSPlot')
            return sData.data(self.__comp)
        else:
            return None

    def gen_params(self, observables = None):
        from ROOT import RooArgSet
        if hasattr(self, '__parameters'):
            return self.__parameters
        else:
            if observables and not isinstance(observables, RooArgSet):
                obs = RooArgSet()
                for o in observables:
                    obs.add(o._target_() if hasattr(o, '_target_') else o)
                observables = obs
            params = self.__pdf.getParameters(observables)
            params.add(self.__status)
            return params

    def result_params(self):
        if not self.__result:
            return []
        else:
            return [p for p in self.__result.floatParsFinal()] + [self.__status]

    def set_params(self, data_params):
        from ROOT import RooCategory
        for trans_param in self.result_params():
            data_param = data_params.find(trans_param.GetName())
            if isinstance(trans_param, RooCategory):
                data_param.setIndex(trans_param.getIndex())
            else:
                data_param.setVal(trans_param.getVal())
                # This sets a symmetric error, but since we don't run Minos, that's ok
                data_param.setError(result_param.getError())
