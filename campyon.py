#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
import getopt
import os
import math 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot



def usage():
    print >>sys.stderr,"Usage: campyon -f filename"
    print >>sys.stderr,"Options:"
    print >>sys.stderr," -k [columns]     columns to keep, see column specification below. All others will be deleted"
    print >>sys.stderr," -d [columns]     columns to delete, see column specification below. All others will be kept"
    print >>sys.stderr," -e [encoding]    Encoding of the input file, defaults to utf-8"
    print >>sys.stderr," -D [delimiter]   Field delimiter (space by default)"
    print >>sys.stderr," -T               Set tab as delimiter"
    print >>sys.stderr," -o [outputfile]  Output to file instead of stdout"
    print >>sys.stderr," -i               Outputfile equals inputfile"
    print >>sys.stderr," -s [expression]  Select rows, expression is may use variables #1...#n for the columns, and operators and,or,not,>,<,!=,== (python syntax)."
    print >>sys.stderr," -S               Compute statistics"
    print >>sys.stderr," -H [columns]     Compute histogram on the specified columns"
    print >>sys.stderr," -C [char]        Ignore comments, line starting with the specified character. Example: -C #"
    print >>sys.stderr," -n               Number lines"
    print >>sys.stderr," -N               Number fields"
    print >>sys.stderr," -M [columns]     Mark/highlight columns"
    print >>sys.stderr," -1               First line is a header line, containing names of the columns"
    print >>sys.stderr," -y [columns]     Plot columns (use with -x)"
    print >>sys.stderr," -x [column]      Plot as a function of the specified column (use with -y)"  
    print >>sys.stderr," -X [samplesizes] Draw one or more random samples (non overlapping, comma separated list of sample sizes)"
    print >>sys.stderr," -A [columns]     Sort by columns, in ascending order"
    print >>sys.stderr," -Z [columns]     Sort by columns, in descending order"
    print >>sys.stderr," -a [column]=[columname]=[expression]   Adds a new column after the specified column"  
    print >>sys.stderr,"Column specification:"
    print >>sys.stderr," A comma separated list of column index numbers or column names (if -1 option is used). Column index numbers start with 1. Negative numbers may be used for end-aligned-indices, where -1 is the last column. Ranges may be specified using a colon, for instance: 3:6 equals 3,4,5,6. A selection like 3:-1 select the third up to the last column. A specification like ID,NAME selects the columns names as such."
    
        
def bold(s):
    CSI="\x1B["
    return CSI+"1m" + s + CSI + "0m"

def red(s):
    CSI="\x1B["
    return CSI+"31m" + s + CSI + "0m"  

def calcentropy(d, base = 2):
      """Compute the entropy of the distribution"""
      entropy = 0
      for type in d:
         if not base:
             entropy += d[type] * -math.log(d[type])     
         else:
             entropy += d[type] * -math.log(d[type], base)     
      return entropy


class CampyonError(Exception):
    pass

class Campyon(object):
    def _parsekwargs(self, key, default, kwargs):
        if key in kwargs:
            return kwargs[key]
        else:
            return default
    
    def parsecolumnindex(self, s):
        if s[0] == '-' and s[1:].isdigit():
            return self.fieldcount + int(s) + 1
        elif s.isdigit():
            return int(s)
        else:
            return self.indexbyname(s)
        
    def parsecolumns(self, settings):
        assert self.fieldcount > 0
        l = []
        for x in settings.split(','):
            if ':' in x:
                low,high = [ y for y in x.split(':') ]
                low = self.parsecolumnindex(low)
                high = self.parsecolumnindex(high) 
                for i in range(low, high + 1):
                    if i > self.fieldcount:  
                        print >>sys.stderr, "ERROR: Specified column " + str(i) + " is out of range"
                        sys.exit(4)
                    l.append(i)
            else:
                x = self.parsecolumnindex(x)
                if x > self.fieldcount:  
                    print >>sys.stderr, "ERROR: Specified column " + str(x) + " is out of range"
                    sys.exit(4)
                l.append(x)
        return l    
    
    def __init__(self, *args, **kwargs):
        try:
	        opts, args = getopt.getopt(args, "f:k:d:e:D:o:is:SH:TC:nNM:1x:y:A:Z:a:",["bar","plotgrid","plotxlog","plotylog","plotconf=","plotfile="])
        except getopt.GetoptError, err:
	        # print help information and exit:
	        print str(err)
	        usage()
	        sys.exit(2)           
        
        self.filename = self._parsekwargs('filename',"",kwargs)
        self.encoding = self._parsekwargs('encoding',"utf-8",kwargs)
        self.delete = self._parsekwargs('delete',[],kwargs)
        self.keep = self._parsekwargs('keep',[],kwargs)
        self.delimiter = self._parsekwargs('delimiter'," ",kwargs)
        self.overwriteinput = self._parsekwargs('overwriteinput',False,kwargs)
        self.outputfile = self._parsekwargs('outputfile',None,kwargs)
        self.DOSTATS = self._parsekwargs('DOSTATS',False,kwargs)
        self.hist = self._parsekwargs('hist',[],kwargs)
        self.sort = self._parsekwargs('sort',[],kwargs)
        self.select = self._parsekwargs('select',None,kwargs)
        self.commentchar = self._parsekwargs('commentchar',None,kwargs)
        self.numberfields = self._parsekwargs('numberfields',False,kwargs)
        self.numberlines =  self._parsekwargs('numberlines',False,kwargs)
        self.highlight =  self._parsekwargs('highlight',[],kwargs)
        self.DOHEADER = self._parsekwargs('DOHEADER',False,kwargs)
        self.x = self._parsekwargs('x',None,kwargs)
        self.y = self._parsekwargs('y',[],kwargs)

        
        self.plotgrid = self._parsekwargs('plotgrid',False,kwargs)
        self.plotxlog = self._parsekwargs('plotxlog',False,kwargs)
        self.plotylog = self._parsekwargs('plotylog',False,kwargs)
        self.plotconf = self._parsekwargs('plotconf',['r','g','b','y','m','c','b'],kwargs)
        self.plotfile = self._parsekwargs('plotfile',"",kwargs)
        
        self.fieldcount = 0
        self.header =  {}
        self.sortreverse = False
        self.inmemory = False
        self.xs = []
        self.ys = {}        
        
        keepsettings = ""
        deletesettings = ""
        histsettings = ""
        highlightsettings = ""
        sortsettings = ""
        plotxsettings = ""
        plotysettings = ""
    
        for o, a in opts:
            if o == "-e":	
                self.encoding = a
            elif o == "-f":	
                self.filename = a              
            elif o == "-k":	
                keepsettings = a
            elif o == "-d":	
                deletesettings = a      
            elif o == '-D':
                self.delimiter = a
            elif o == '-o':    
                self.outputfile = a
            elif o == '-s':
                self.select = a
            elif o == '-i':
                self.outputfile = self.filename + ".tmp"
                self.overwriteinput = True
            elif o == '-S':
                self.DOSTATS = True
            elif o == '-H':
                histsettings = a 
            elif o == '-T':
                self.delimiter = "\t"
            elif o == '-C':
                self.commentchar = a
            elif o == '-n':
                self.numberlines = True
            elif o == '-N':
                self.numberfields = True        
            elif o == '-M':    
                highlightsettings = a
            elif o == '-1':
                self.DOHEADER = True
            elif o == '-A':
                sortsettings = a
            elif o == '-Z':
                sortsettings = a
                self.sortreverse = True     
            elif o == '-x':           
                plotxsettings = a                
            elif o == '-y':           
                plotysettings = a           
            elif o == '--plotgrid':     
                self.plotgrid = True
            elif o == '--plotxlog':     
                self.plotxlog = True                
            elif o == '--plotylog':     
                self.plotylog = True
            elif o == '--plotconf':     
                self.plotconf = a.split(',')
            elif o == '--plotfile':     
                self.plotfile = a
            elif o == '-a':
                raise NotImplementedError
            else:
                raise Exception("invalid option: " + o)
                        
        if not self.filename:    
            usage()
            sys.exit(2)

        if not os.path.exists(self.filename):
            print >>sys.stderr,"No such file: " + self.filename
            sys.exit(2)
        
        f = codecs.open(self.filename,'r',self.encoding)
        for line in f:
            if line.strip() and (not self.commentchar or line[:len(self.commentchar)] != self.commentchar):                    
                fields = line.strip().split(self.delimiter)
                self.fieldcount = len(fields)
                print >>sys.stderr,"Number of fields: ", self.fieldcount
                if self.DOHEADER:
                    self.header = dict([ (x+1,y) for x,y in enumerate(fields) ])
                break            
        f.close()
        
        
        
        if keepsettings: self.keep = self.parsecolumns(keepsettings)
        if deletesettings: self.delete = self.parsecolumns(deletesettings)
        if histsettings: self.hist = self.parsecolumns(histsettings)
        if highlightsettings: self.highlight = self.parsecolumns(highlightsettings)
        if sortsettings: self.sort = self.parsecolumns(sortsettings)
        if plotxsettings: self.x = self.parsecolumnindex(plotxsettings)
        if plotysettings: self.y = self.parsecolumns(plotysettings)           
           
        if self.sort:
            self.inmemory = True
        

        
        
            
        self.memory = []    
        self.sumdata = {}
        self.nostats = set()
        self.freq = {}
        
        
        if self.keep:
            print >>sys.stderr, "Fields to keep: ",  " ".join([ str(x) for x in self.keep])
        if self.delete:
            print >>sys.stderr, "Fields to delete: ",  " ".join([ str(x) for x in self.delete])
        
        self.rowcount_in = 0
        self.rowcount_out = 0
        
    def __call__(self):        
        if self.outputfile:
            f_out = codecs.open(self.outputfile, 'w',self.encoding)
            
        for line, fields, linenum in self:
            if self.outputfile:                                          
                if self.numberlines: f_out.write("@" + str(linenum) + self.delimiter)
                f_out.write(line + "\n")
            else:
                if self.numberlines: print "@" + str(linenum) + self.delimiter,
                print line.encode(self.encoding)
              
        if self.outputfile:
            f_out.close()
        
        if self.overwriteinput:
            os.rename(self.outputfile,self.filename)
            
        if self.DOSTATS:            
            self.printstats()
            
        if self.hist:
            for fieldnum in sorted(self.freq):
                print >>sys.stderr, "Histogram for column #" + str(fieldnum) + "\ttypes=" + str(self.types(fieldnum)) + "\ttokens=" + str(self.tokens(fieldnum)) + "\tttr=" +  str(self.ttr(fieldnum)) + "\tentropy=" + str(self.entropy(fieldnum))
                print >>sys.stderr,"------------------------------------------------------------------------"
                self.printhist(fieldnum)
                
        if self.x and self.y:
            self.plot()
                 
    def __len__(self):        
        return self.rowcount_out
        
    def __iter__(self):            
        self.sumdata = {}
        self.nostats = set()
        self.freq = {}                
        self.rowcount_in = 0
        self.rowcount_out = 0 
        
        if self.keep: 
            default = 'delete'
        else:        
            default = 'keep'        
            
        headerfound = False
                           
        f = codecs.open(self.filename,'r',self.encoding)
        for line in f:
            self.rowcount_in += 1
            
            
            if not line.strip() or (self.commentchar and line[:len(self.commentchar)] == self.commentchar):
                self.rowcount_out += 1
                if not self.inmemory:
                    yield line.strip(), [], self.rowcount_out
                continue
            
            
            

                
            fields = line.strip().split(self.delimiter)
            if len(fields) != self.fieldcount:
                raise CampyonError("Number of columns in line " + str(self.rowcount_in) + " deviates, expected " + str(self.fieldcount) + ", got " + str(len(fields))) 
                
            if not headerfound:
                headerfound = True
                self.rowcount_out += 1
                if not self.inmemory:
                    yield line.strip(), fields, self.rowcount_out
                continue
    
        
            
            if self.select:
                currentselect = self.select
                for i in reversed(range(1,len(fields)+1)):
                    isdigit = True
                    try:
                        x = float(fields[i-1])
                    except:                    
                        isdigit = False
                    if isdigit:
                        currentselect = currentselect.replace('#' + str(i), fields[i-1])
                    else:
                        currentselect = currentselect.replace('#' + str(i), '"' + fields[i-1].replace('"',"\\\"") + '"')
                if not eval(currentselect):
                    continue
        
            self.rowcount_out += 1
            

            if self.hist:
                for fieldnum in self.hist:
                    if not fieldnum in self.freq:                
                        self.freq[fieldnum] = {}
                    if not fields[fieldnum] in self.freq:
                        self.freq[fieldnum][fields[fieldnum]] = 0
                    self.freq[fieldnum][fields[fieldnum]] += 1
                
            if self.DOSTATS:                    
                for i,field in enumerate(fields):
                    if not i in self.nostats: 
                        if '.' in field:
                            try:
                               x = float(field)
                            except:
                               self.nostats.add(i+1)
                               if i in self.sumdata: del self.sumdata[i+16]
                               continue
                        else:
                            try:
                               x = int(field)
                            except:
                               self.nostats.add(i+1)
                               if i in self.sumdata: del self.sumdata[i+1]
                               continue                        
                     
                        if not i in self.sumdata:
                            self.sumdata[i+1] = 0
                        self.sumdata[i+1] += x
                         
                
            
            newfields = []
            #k = [ x - 1 if x >= 0 else len(fields) + x for x in keep ]
            #d = [ x - 1 if x >= 0 else len(fields) + x for x in delete ]
            for i, field in enumerate(fields):
                fieldnum = i+1
                action = default 
                if fieldnum in self.keep:
                    action = 'keep'
                elif fieldnum in self.delete:
                    action = 'delete'
                if fieldnum in self.highlight:
                    field = bold(red(field))                
                if self.numberfields:
                    field = str(i) + '=' + field
                else:                
                    if field.isdigit() or field[0] == '-' and field[1:].isdigit():
                        field = int(field)
                    else:
                        isfloat = True
                        try:
                            f = float(field)
                        except:
                            isfloat = False
                        if isfloat:
                            field = f

                if self.x == fieldnum:
                    self.xs.append(field)
                
                if fieldnum in self.y:
                    if not isinstance(field, float) and not isinstance(field,int):
                        raise CampyonError("Can not plot non-numeric values: " + field) 
                        
                    if not fieldnum in self.ys:
                        self.ys[fieldnum] = []
                    self.ys[fieldnum].append(field)

                if action == 'keep':
                    newfields.append(field)
                    
            s = self.delimiter.join([ str(x) for x in newfields ])            
            if self.inmemory:
                self.memory.append( (newfields, self.rowcount_out) )
            else:                
                yield s, newfields, self.rowcount_out

        if self.inmemory:
            if self.sort:
               self.memory = sorted(self.memory, key=lambda x: tuple([ x[0][i-1] for i in self.sort ]), reverse=self.sortreverse)
                              
            if self.header:    
                s = self.delimiter.join( self.headerfields()  )
                yield s, self.headerfields(), 0
    
            for fields, linenum in self.memory:
                s = self.delimiter.join([ str(x) for x in fields])
                yield s, fields, linenum
            
        print >>sys.stderr,"Read " + str(self.rowcount_in) + " lines, outputted " + str(self.rowcount_out)
        
    def plot(self):        
        fig = matplotlib.pyplot.figure()
        if self.plotgrid:
            matplotlib.pyplot.grid(True)
        else:
            matplotlib.pyplot.grid(False)            
        if self.plotylog:            
            fig.set_yscale('log')
        if all([ isinstance(x,float) or isinstance(x,int) for x in self.xs ]):
            if self.plotxlog:
                fig.set_xscale('log')

            l = []
            for i, field in enumerate(self.y):
                l.append(self.xs)
                l.append(self.ys[field])
                l.append(self.plotconf[i])
            fig.plot(*l)
        else:            
           #TODO: implement bar chart
           raise NotImplementedError
            
        if self.plotfile:
            fig.savefig(self.plotfile, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        else:
            fig.show()
        
    def headerfields(self):          
        return [x[1] for x in sorted(self.header.items()) ]
        
    def printstats(self, out=sys.stderr):            
        for colnum, colname, s, average in self.stats():
            if colname == str(colnum): 
                colname == ""
            else:
                colname = '-' + colname
            out.write("column #" + str(colnum) + colname.encode(self.encoding) + " sum="+ str(s) + "\taverage=" + str(average) + "\n")
            
    def stats(self):
        for i in sorted(self.sumdata):
            if self.header:
                colname = self.header[i]
            else:
                colname = str(i)   
            
            yield i, colname, self.sumdata[i], self.sumdata[i] / float(self.rowcount_out)
                         
     
    def printhist(self, columnindex, out=sys.stderr):         
        for i, (word, count, f) in self.histdata(columnindex):
            print >>sys.stderr, str(i) + ")\t" + word.encode(self.encoding) + "\t" + str(count) + "\t" + str(f * 100) + '%'
        
    def entropy(self, columnindex):        
        return calcentropy(self.freq[columnindex])
    
    def tokens(self, columnindex):
        return sum(self.freq[columnindex].values())
        
    def types(self, columnindex):
        return len(self.freq[columnindex])
    
    def ttr(self, columnindex):
        return self.types(columnindex) / float(self.tokens(columnindex))
        
    def histdata(self, columnindex):   
        s = float(self.tokens(columnindex))
        for word, count in sorted(self.freq[columnindex].items(), key=lambda x: x[1] * -1):
            yield word, count, count / s
    
    def indexbyname(self, colname):
        for key, i in self.header.items():
            if key == colname:
                return i
        raise KeyError("Column " + colname + " not found")
        
if __name__ == "__main__":         
    campyon = Campyon(*sys.argv[1:])
    campyon()
            
