from math import pi

g = GaudiPython()
g.script = ('UnbiasedProfiles_Coeff.py')

###############################
# Make loop for arguments
###############################
npoints = 40

list = []

for i in range(0,npoints):
    list.append(['-x',str(i),'-n',str(npoints)])

print 'list =',list
J = Job(name = 'UnbiasedProfile_Zoom_Coeff_%s_steps'%(npoints), application = g, backend=PBS())
s = ArgSplitter(args=list)
J.splitter = s
J.inputsandbox = ['libp2vv.so','RooFitDecorators.py','ModelBuilders.py','UnbiasedWS_Coeff.root']
J.outputsandbox = ['profile.root']
J.submit()
