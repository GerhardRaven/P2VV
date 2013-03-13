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
        while self.locked(cf):
            sleep(5)
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
