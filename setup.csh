setenv P2VVROOT $PWD

if ( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH $P2VVROOT/lib:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH $P2VVROOT/lib
endif

if ( $?PYTHONPATH ) then
  setenv PYTHONPATH $P2VVROOT/python:${PYTHONPATH}
else
  setenv PYTHONPATH $P2VVROOT/python
endif
