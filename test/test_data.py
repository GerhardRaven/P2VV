class TestData(object):

    def setup_class(cls):
        # setup (singleton) workspace
        from ROOT import RooWorkspace
        from RooFitWrappers import RooObject
        ws = RooObject()
        ws.setWorkspace(RooWorkspace("test_workspace"))

    def test_mapping(self):
        # now (consistently!) create/declare observables
        import os
        import inspect
        from ROOT import TFile
        from RooFitWrappers import RealVar
        from Helpers import Mapping

        m = RealVar('m', Observable = True, Unit = 'MeV/c^2', MinMax=(5000, 6000))
        t = RealVar('t', Observable = True, MinMax = (-1, 14), Unit = 'ps', Value = 0)

        import inspect
        thisFile = inspect.getfile(inspect.currentframe())
        thisDir = os.path.dirname(thisFile)
        rootFile = TFile.Open(os.path.join(thisDir, 'data.root'))
        dataset = rootFile.Get('data')
        
        m = Mapping({m : 'data_m', t :'data_t'}, dataset)

        m = RealVar('m', Observable = True, Unit = 'MeV/c^2', MinMax=(5000, 6000))
        assert m['Observable'] == True
        assert m['Unit'] == 'MeV/c^2'
        assert m['MinMax'] == (5000, 6000)
        assert m['Name'] == 'data_m'

        t = RealVar('t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)
        assert t['Name'] == 'data_t'
        assert t['Observable'] == True
        assert t['Unit'] == 'ps'
        assert t['MinMax'] == (-1,14)
