#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
import getopt
import os
import math 




def usage():
    print >>sys.stderr,"Usage: campyon -f filename"
    print >>sys.stderr,"Options:"
    print >>sys.stderr," -k [columns]     columns to keep (1-indexed), this is a comma separated list of column numbers, negative numbers are allowed for end-aligned-indices. All others will be deleted"
    print >>sys.stderr," -d [columns]     columns to delete (1-indexed), this is a comma separated list of column numbers, negative numbers are allowed for end-aligned-indices. All others will be kept"
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
    print >>sys.stderr," -p [columns]     Plot columns"

        
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


def parsecolumns(settings, fieldcount):
    assert fieldcount > 0
    print >>sys.stderr, "DEBUG=",fieldcount
    l = []
    for x in settings.split(','):
        if ':' in x:
            low,high = [ int(y) for y in x.split(':') ]                
            if low < 0: low = fieldcount + low + 1
            if high < 0: high = fieldcount + high + 1
            for i in range(low, high + 1):
                if i > fieldcount:  
                    print >>sys.stderr, "ERROR: Specified column " + str(i) + " is out of range"
                    sys.exit(4)
                l.append(i)
        else:
            if int(x) < 0: x = fieldcount + int(x) + 1
            if int(x) > fieldcount:  
                print >>sys.stderr, "ERROR: Specified column " + str(x) + " is out of range"
                sys.exit(4)
            l.append(int(x))
    return l

class CampyonError(Exception):
    pass

class Campyon(object):
    def __init__(self, *args, **kwargs):
        try:
	        opts, args = getopt.getopt(args, "f:k:d:e:D:o:is:SH:TC:nNM:1")
        except getopt.GetoptError, err:
	        # print help information and exit:
	        print str(err)
	        usage()
	        sys.exit(2)           
        
        self.filename = ""
        self.encoding = "utf-8"
        keepsettings = ""
        deletesettings = ""
        histsettings = ""
        self.delete = []
        self.keep = []
        self.delimiter = " "
        self.overwriteinput = False
        self.outputfile = None
        self.DOSTATS = False
        self.hist = []
        self.select = None
        self.fieldcount = 0
        self.commentchar = None
        self.numberfields = False
        self.numberlines = False
        highlightsettings = ""
        self.highlight = []
        self.DOHEADER = False
        self.header = []
    
        for o, a in opts:
            if o == "-e":	
                self.encoding = a
            elif o == "-f":	
                self.filename = a              
            elif o == "-k":	
                self.keepsettings = a
            elif o == "-d":	
                self.deletesettings = a      
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
                self.histsettings = a 
            elif o == '-T':
                self.delimiter = "\t"
            elif o == '-C':
                self.commentchar = a
            elif o == '-n':
                self.numberlines = True
            elif o == '-N':
                self.numberfields = True        
            elif o == '-M':    
                self.highlightsettings = a
            elif o == '-1':
                self.DOHEADER = True
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
                self.fieldcount = len(line.strip().split(self.delimiter))
                print >>sys.stderr,"Number of fields: ", self.fieldcount
                break
        f.close()
        
        
        
        if keepsettings: self.keep = parsecolumns(keepsettings, self.fieldcount)
        if deletesettings: self.delete = parsecolumns(deletesettings, self.fieldcount)
        if histsettings: self.hist = parsecolumns(histsettings, self.fieldcount)
        if highlightsettings: self.highlight = parsecolumns(highlightsettings, self.fieldcount)
           
        

        
        
            
            
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
            
        for line, fields in self:
            if self.outputfile:                                          
                if self.numberlines: f_out.write("@" + str(self.rowcount_in) + self.delimiter)
                f_out.write(line + "\n")
            else:
                if self.numberlines: print "@" + str(self.rowcount_in) + self.delimiter,
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
                           
        f = codecs.open(self.filename,'r',self.encoding)
        for line in f:
            self.rowcount_in += 1
            
            
            if not line.strip() or (self.commentchar and line[:len(self.commentchar)] == self.commentchar):
                self.rowcount_out += 1
                yield line.strip()
                continue
            
            
            

                
            fields = line.strip().split(self.delimiter)
            if len(fields) != self.fieldcount:
                raise CampyonError("Number of columns in line " + str(self.rowcount_in) + " deviates, expected " + str(self.fieldcount) + ", got " + str(len(fields))) 
                
            
            if self.DOHEADER and not self.header:
                self.header = fields
            
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
                action = default 
                if i in self.keep:
                    action = 'keep'
                elif i in self.delete:
                    action = 'delete'
                if i in self.highlight:
                    field = bold(red(field))                
                if self.numberfields:
                    field = str(i) + '=' + field                
                if action == 'keep':
                    newfields.append(field)
            s = self.delimiter.join(newfields)
            
            yield s, newfields

            
        print >>sys.stderr,"Read " + str(self.rowcount_in) + " lines, outputted " + str(self.rowcount_out)
        
        
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
        return calcentropy(self.freq[columnindex-1])
    
    def tokens(self, columnindex):
        return sum(self.freq[columnindex-1].values())
        
    def types(self, columnindex):
        return len(self.freq[columnindex-1])
    
    def ttr(self, columnindex):
        return self.types(columnindex) / float(self.tokens(columnindex))
        
    def histdata(self, columnindex):   
        s = float(self.tokens(columnindex-1))
        for word, count in sorted(self.freq[columnindex].items(), key=lambda x: x[1] * -1):
            yield word, count, count / s
    
if __name__ == "__main__":         
    campyon = Campyon(*sys.argv[1:])
    campyon()
            
