# prof.pyx
cdef extern from "google/profiler.h":
    void ProfilerStart( char* fname )
    void ProfilerStop()
 
def profiler_start(fname):
    ProfilerStart(<char *> fname)
 
def profiler_stop():
    ProfilerStop()

cdef extern from "google/heap-profiler.h":
    void HeapProfilerStart( char* fname )
    void HeapProfilerStop()

def heap_profiler_start(fname):
    HeapProfilerStart(<char *> fname)

def heap_profiler_stop():
    HeapProfilerStop()
