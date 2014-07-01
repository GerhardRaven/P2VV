import os
import re

class CacheFiles(object):
    def __init__(self, directory, filename):
        self._filename = filename
        self._directory = directory
        
        if self._filename.endswith('.root'):
            self._filename = self._filename[:-5]

        if not os.path.exists(directory):
            raise OSError("directory %s does not exist" % directory)
        self._rootFiles = []

        self._updateFiles()

    def _updateFiles(self):
        self._files = []
        self._regex = re.compile('%s_(?P<index>[\w_0-9]+)\.root' % self.filename())
        for f in os.listdir(self.directory()):
            if self._regex.match(f):
                self._files.append(f)
        self._files = sorted(self._files)

        from ROOT import TFile
        rf = set(f.GetName() for f in self._rootFiles)
        for f in self._files:
            f = os.path.join(self.directory(), f)
            if f not in rf:
                tf = TFile.Open(os.path.join(self.directory(), f))
                self._rootFiles.append(tf)
        
    def getCacheFiles(self):
        self._updateFiles()
        return self._rootFiles

    def filename(self):
        return self._filename

    def directory(self):
        return self._directory

    def getFromCache(self, path):
        for cf in self.getCacheFiles():
            o = cf.Get(path)
            if o: return o, cf
        else:
            return (None, None)

class WritableCacheFile(CacheFiles):
    def __init__(self, arg1, arg2):
        if type(arg1) != str:
            self.__initCacheFiles(arg1, arg2)
        else:
            self.__initDirFile(arg1, arg2)

    def __initCacheFiles(self, cf, directory):
        self._file_dir = directory
        self._filename = cf.filename()
        self._directory = cf.directory()
        
        if not os.path.exists(self._directory):
            raise OSError("directory %s does not exist" % self._directory)
        
        self._rootFiles = []
        self._updateFiles()
        self._rootFile = None
        self._rootFile = self._getWritableCacheFile()

    def __initDirFile(self, directory, filename):
        CacheFiles.__init__(self, directory, filename)
        self._file_dir = None
        self._rootFile = None
        print 'Get writable file'
        self._rootFile = self._getWritableCacheFile()

    def _getFirstCacheFile(self):
        if self._rootFiles:
            cf = self._rootFiles[0]
        else:
            from ROOT import TFile
            cf = TFile.Open(os.path.join(self._directory, '%s_cache.root' % self._filename), 'new')
            self._rootFiles.append(cf)
        return cf
    
    def _getWritableCacheFile(self):
        if self._file_dir:
            if not self._rootFile:
                co, cf = self.getFromCache(self._file_dir)
                if not cf:
                    cf = self._getFirstCacheFile()
            else:
                return self._rootFile
        else:
            if self._rootFile:
                if not self._rootFile.Get(self._file_dir):
                    self._rootFile.mkdir(self._file_dir)
                return self._rootFile
            else:
                cf = self._getFirstCacheFile()
        
        from time import sleep
        printed = False
        while self.locked(cf):
            sleep(5)
            if not printed:
                print 'Waiting for lock file'
                printed = True
        self._rootFile = self._openForWriting(cf)
        if self._file_dir and not self._rootFile.Get(self._file_dir):
            self._rootFile.mkdir(self._file_dir)
        return self._rootFile

    def _openForWriting(self, rf):
        self._rootFiles.remove(rf)
        rf.Close()
        tmp = rf.Open(rf.GetName(), "update")
        self.lock(tmp)
        self._rootFiles.append(tmp)
        return tmp
        
    def getCacheFiles(self):
        self._updateFiles()
        if self._rootFile:
            return [self._rootFile] + self._rootFiles
        else:
            return self._rootFiles

    def getCacheFile(self):
        return self._rootFile

    def _lockName(self, root_file):
        filename = root_file.GetName().rsplit('/', 1)[-1]
        lock_name = ('.lock_%d' % hash(filename)).replace('-', 'm')
        return os.path.join(self.directory(), lock_name)

    def locked(self, root_file):
        return os.path.exists(self._lockName(root_file))

    def lock(self, root_file):
        from time import sleep
        lock_name = self._lockName(root_file)
        while os.path.exists(lock_name):
            sleep(5)
        with open(lock_name, 'w') as f:
            pass

    def unlock(self, root_file):
        lock_name = self._lockName(root_file)
        if os.path.exists(lock_name):
            os.remove(lock_name)

    def __enter__(self):
        return self._rootFile
    
    def __exit__(self, type, value, traceback):
        ## Files open for reading are not closed, ROOT will deal with those in
        ## the end.
        if self._rootFile:
            self._rootFile.Write()
            self._rootFile.Purge()
            self._rootFile.Flush()
            self._rootFile.Close()
            self.unlock(self._rootFile)

class Cache(object):
    def __init__(self, cache_def, directory):
        from P2VV.CacheUtils import CacheFiles
        self.__cache_files = CacheFiles(*cache_def)
        self.__directory = directory
        self.__cache_dir, self.__cache_file = self.__cache_files.getFromCache(directory)

    def directory(self):
        return self.__directory

    def cache_files(self):
        return self.__cache_files

    def cache_dir(self):
        return self.__cache_dir
        
    def _get_make_dir(self, d, top = None):
        if top == None:
            top = self.cache_dir()
        tmp = top.Get(d)
        if not tmp:
            top.mkdir(d)
            tmp = top.Get(d)
        return tmp

    def read_data(self):
        ## Read sdata
        cache_dir = self.cache_dir()
        if not cache_dir:
            return None, {}
        data_dir = cache_dir.Get('data')
        if not data_dir:
            return None, {}
        keys = data_dir.GetListOfKeys()
        if keys.GetEntries() != 1:
            return None, {}
        else:
            data = data_dir.Get(keys.At(0).GetName())

        sdata_dir = cache_dir.Get('sdata')
        sdatas = {}
        if not sdata_dir:
            return data, sdatas
        dss = []
        for e in sdata_dir.GetListOfKeys():
            if e.GetClassName() == 'RooDataSet':
                dss.append(e.GetName())
        for e in dss:
            sdata = sdata_dir.Get(e)
            if not sdata:
                return data, {}
            else:
                sdatas[e] = sdata
        return data, sdatas

    def read_results(self):
        # Read results
        results = {}
        cd = self.cache_dir()
        if not cd:
            return results
        rd = cd.Get('results')
        if not rd:
            return {}
        for key in rd.GetListOfKeys():
            if key.GetClassName() != 'RooFitResult':
                continue
            results[key.GetName()] = rd.Get(key.GetName())
        return results

    def write_cut(self, cut):
        with WritableCacheFile(self.cache_files(), self.directory()) as cache_file:
            cache_dir = cache_file.Get(self.directory())
            from ROOT import TObjString
            cut_string = TObjString(cut)
            cache_dir.WriteTObject(cut_string, 'cut')

    def write_data(self, data, sdatas):
        with WritableCacheFile(self.cache_files(), self.directory()) as cache_file:
            cache_dir = cache_file.Get(self.directory())

            from ROOT import TObject    
            sdata_dir = self._get_make_dir('sdata', cache_dir)
            sdata_sub_dir = None
            for name, ds in sdatas.iteritems():
                split_name = name.split('/')
                if len(split_name) > 1:
                    sub_dir, name = split_name
                    sdata_sub_dir = self._get_make_dir(sub_dir + '/sdata', cache_dir)
                    sdata_sub_dir.WriteTObject(ds, name, 'Overwrite')
                else:
                    sdata_dir.WriteTObject(ds, name, "Overwrite")
            if sdata_sub_dir:
                sdata_sub_dir.Write(sdata_sub_dir.GetName(), TObject.kOverwrite)
            sdata_dir.Write(sdata_dir.GetName(), TObject.kOverwrite)

            data_dir = self._get_make_dir('data', cache_dir)
            data_dir.WriteTObject(data, data.GetName(), "Overwrite")
            data_dir.Write(data_dir.GetName(), TObject.kOverwrite)

            # Delete the input TTree which was automatically attached.
            cache_file.Delete('%s;*' % 'DecayTree')

    def write_results(self, results):
        from ROOT import TObject
        with WritableCacheFile(self.cache_files(), self.directory()) as cache_file:
            cache_dir = cache_file.Get(self.directory())
            ## Write mass fit results
            result_dir = self._get_make_dir('results', cache_dir)
            result_sub_dir = None
            for name, result in results.iteritems():
                split_name = name.split('/')
                if len(split_name) > 1:
                    sub_dir, name = split_name
                    result_sub_dir = self._get_make_dir(sub_dir + '/results', cache_dir)
                    result_sub_dir.WriteTObject(result, name, 'Overwrite')
                else:
                    result_dir.WriteTObject(result, name, "Overwrite")
            if result_sub_dir:
                result_sub_dir.Write(result_sub_dir.GetName(), TObject.kOverwrite)
            result_dir.Write(result_dir.GetName(), TObject.kOverwrite)

    def write_plots(self, plots):
        from ROOT import TObject
        with WritableCacheFile(self.cache_files(), self.directory()) as cache_file:
            cache_dir = cache_file.Get(self.directory())
            ## Write plots
            dirs = {}
            for d in plots.iterkeys():
                if not d:
                    dirs[d] = self._get_make_dir('plots', cache_dir)
                else:
                    sub_dir = self._get_make_dir(d, cache_dir)
                    dirs[d] = self._get_make_dir('plots', sub_dir)
            for d, ps in plots.iteritems():
                for p in ps:
                    dirs[d].WriteTObject(p[0], p[0].GetName(), 'overwrite')
            for d in dirs.itervalues():
                d.Write(d.GetName(), TObject.kOverwrite)
                
class SimCache(Cache):
    def __init__(self, cache_def, directory, sub_dir):
        self.__sub_dir = sub_dir
        Cache.__init__(self, cache_def, directory)
        
    def read_results(self):
        cache_dir = self.cache_dir()
        results = {}
        if not cache_dir:
            return results
        for d in ('', self.__sub_dir + '/'):
            rd = cache_dir.Get(d + 'results')
            if not rd:
                return {}
            for key in rd.GetListOfKeys():
                if key.GetClassName() != 'RooFitResult':
                    continue
                results[d + key.GetName()] = rd.Get(key.GetName())
        return results

    def read_data(self):
        cache_dir = self.cache_dir()
        if not cache_dir:
            return None, {}
        data_dir = cache_dir.Get('data')
        if not data_dir:
            return None, {}
        keys = data_dir.GetListOfKeys()
        if keys.GetEntries() != 1:
            return None, {}
        else:
            data = data_dir.Get(keys.At(0).GetName())
        
        ## Read sdata
        cache_dir = self.cache_dir()
        sdatas = {}
        for d in ('', self.__sub_dir + '/'):
            sdata_dir = cache_dir.Get(d + 'sdata')
            if not sdata_dir:
                return data, {}
            for key in sdata_dir.GetListOfKeys():
                if key.GetClassName() != 'RooDataSet':
                    continue
                sdatas[d + key.GetName()] = sdata_dir.Get(key.GetName())
        return data, sdatas
