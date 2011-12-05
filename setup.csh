if ( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH $PWD/lib:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH $PWD/lib
endif

if ( $?PYTHONPATH ) then
  setenv PYTHONPATH $PWD/python:${PYTHONPATH}
else
  setenv PYTHONPATH $PWD/python
endif
