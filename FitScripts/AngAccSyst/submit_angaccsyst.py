import os

acceptance = 1

namelist = ['000plus','000minus','020plus','020minus','022plus','022minus','200plus','200minus']

if acceptance:
    fitname = 'cfit_acc'
else:
    fitname = 'cfit_noacc'

scriptname = 'Thesis_angaccsyst.py'

for name in namelist:
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
    j.inputsandbox = ['/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/angaccsyst/%s'%(scriptname)
                      ,'/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/angaccsyst/%s_result'%(fitname)
                      ,'/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/angaccsyst/%sparams.txt'%(fitname)
                      ]

    # Add the outputsandbox
    j.outputsandbox = ['*.txt','*_result']

    # Add the splitter
    list = []
    list.append([scriptname,'-n',name])
    print list
    s = ArgSplitter(args=list)
    j.splitter = s

    j.name = 'angaccsyst_%s'%(name)
        
    # backend
    #j.backend = PBS(queue = 'qlong')
    j.backend = PBS()
    
    j.submit()
