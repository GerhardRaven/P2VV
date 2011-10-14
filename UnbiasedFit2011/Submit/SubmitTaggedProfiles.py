from math import pi

g = GaudiPython()
g.script = ('TaggedProfiles.py')
#g.project = 'Erasmus'
#g.version = 'v3r1'

#g.user_release_area = '/afs/cern.ch/user/d/dvaneijk/cmtuser'

###############################
# Make loop for arguments
###############################
npoints = 40

list = []

for j in range(0,npoints):
    #list.append(['-p',str(Phi_s_in),'-s',str(j+1)])
    list.append(['-j',str(j),'-n',str(npoints)])

print 'list =',list
J = Job(name = 'TaggedProfile_%s_steps_floatAs2'%(npoints), application = g, backend=PBS())
s = ArgSplitter(args=list)
J.splitter = s
#J.outputdata = ['testdata.root']
J.inputsandbox = ['libp2vv.so','RooFitDecorators.py','ModelBuilders.py','TaggedWS.root']
J.outputsandbox = ['profile.root']
J.submit()

    
