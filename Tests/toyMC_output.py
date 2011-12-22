from optparse import OptionParser
import sys
parser = OptionParser(usage = '%prog jobid nbins')

(options, args) = parser.parse_args()
if len(args) != 2:
    print parser.usage
    sys.exit(-1)
jobid = int(args[0])
nbins = int(args[1])

from ROOT import (RooGenFitStudy, RooStudyManager,
                  RooFit, TFile)

root_file = TFile.Open('workspace2D_%d.root' % nbins)
w = root_file.Get('w')

gfs = RooGenFitStudy();
gfs.setGenConfig("acc_pdf", "m,t", RooFit.NumEvents(10000))
gfs.setFitConfig("acc_pdf", "m,t", RooFit.Minimizer("Minuit2"))

mgr = RooStudyManager(w, gfs)

## mgr.run(1)
import os
basedir = '/opt/project/bfys/raaij/gangadir/workspace/raaij/LocalXML/%d' % jobid
suffix = 'output/study_result_data2D_%d_*.root' % nbins
for d in os.listdir(basedir):
    try:
        int(d)
        f = os.path.join(basedir, d, suffix)
        mgr.processBatchOutput(f)
    except ValueError:
        pass
    
data = gfs.summaryData()
output_file = TFile.Open('data2D_%d.root' % nbins, 'recreate')
output_file.WriteTObject(data, 'data')
output_file.Close()
