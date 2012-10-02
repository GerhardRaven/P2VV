import os

acceptance = 0

if acceptance:
    fitname = 'cfit_acc_floatdeltams'
else:
    fitname = 'cfit_noacc_floatdeltams'

scriptname = 'Profiles_floatdeltams.py'
nsigma = 17
npoints = 40
paramlist = ['dM']
    
for param in paramlist:
    # Create the job
    j = Job()
    j.application = Executable()
    j.application.exe = 'python'

    # Create the right environment
    env = {}
    root_location = '/project/bfys/jleerdam/local/root/5.32.00-patches_P2VV'
    env_vars = {'PATH' : 'bin',
                'PYTHONPATH' : 'lib',
                'LD_LIBRARY_PATH' : 'lib',
                'ROOTSYS' : ''}
    for var, d in env_vars.iteritems():
        val = os.environ[var]
        s = val.split(':')
        found = False
        for k in s:
            if k.startswith(root_location):
                found = True
                break
        if not d:
            env[var] = root_location
        elif found:
            env[var] = val
        else:
            env[var] = os.path.join(root_location, d) + ':' + val
    j.application.env = env

    # Add the inputsandbox
    j.inputsandbox = ['/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/%s'%(scriptname)
                      ,'/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/%s_result'%(fitname)
                      ,'/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/%sparams.txt'%(fitname)
                      ]

    # Add the outputsandbox
    j.outputsandbox = ['*.txt']


    # Add the splitter
    list = []
    for i in range(0,npoints):
        list.append([scriptname, '-x',str(i),'-n',str(npoints),'-v',param,'-s',nsigma,'-a',acceptance])

    s = ArgSplitter(args=list)
    j.splitter = s

    j.name = '%s_Profiles_%s'%(param,fitname)

    # backend
    #j.backend = PBS(queue = 'stbcq')
    #j.backend = PBS(queue = 'qlong')
    j.backend = PBS()

    j.submit()
