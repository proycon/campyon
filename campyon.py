#! /usr/bin/env python
# -*- coding: utf8 -*-



###############################################################
# CAMPYON
#       by Maarten van Gompel (proycon)
#       http://github.com/proycon/campyon
#
#       Centre for Language Studies
#       Radboud University Nijmegen
#       
#       Licensed under GPLv3
#
# Campyon is a command-line tool and Python library for viewing and manipulating columned data files.
# It supports various filters, statistics, and plotting.
#
###############################################################


import sys
import codecs
import getopt
import os
import math 
import re


if '-x' in sys.argv[1:]: #don't import if not used, to save time
    import numpy
    import matplotlib
    matplotlib.use('GTKAgg')
    import matplotlib.pyplot
if '-V' in sys.argv[1:]: #don't import if not used, to save time
    import gtk


def usage():
    print >>sys.stderr,"Campyon - by Maarten van Gompel - http://github.com/proycon/campyon"
    print >>sys.stderr," Campyon is a command-line tool and Python library for viewing and manipulating columned data files."
    print >>sys.stderr," It supports various filters, statistics, and plotting."
    print >>sys.stderr,"Usage: campyon [options] [filename] [[filename2]...]"
    print >>sys.stderr,"Options:"
    print >>sys.stderr," -k [columns]     columns to keep, see column specification below. All others will be deleted"
    print >>sys.stderr," -d [columns]     columns to delete, see column specification below. All others will be kept"
    print >>sys.stderr," -e [encoding]    Encoding of the input file, defaults to utf-8"
    print >>sys.stderr," -D [delimiter]   Field delimiter (space by default)"
    print >>sys.stderr," -T               Set tab as delimiter"
    print >>sys.stderr," -o [outputfile]  Output to file instead of stdout (will aggregate in case multiple input files are specified)"
    print >>sys.stderr," -i               Outputfile equals inputfile"
    print >>sys.stderr," -s [expression]  Select rows, see section on selector specification below for syntax"
    print >>sys.stderr," -S               Compute statistics"
    print >>sys.stderr," -H [columns]     Compute histogram on the specified columns"
    print >>sys.stderr," -C [char]        Ignore comments, line starting with the specified character. Example: -C #"
    print >>sys.stderr," -n               Number lines"

    print >>sys.stderr," -N               Number fields"
    print >>sys.stderr," -M [columns]     Mark/highlight columns"
    print >>sys.stderr," -1               First line is a header line, containing names of the columns"
    print >>sys.stderr," -y [columns]     Plot columns (use with -x)"
    print >>sys.stderr," -x [column]      Plot as a function of the specified column (use with -y)"  
    
    print >>sys.stderr," -A [columns]     Sort by columns, in ascending order"
    print >>sys.stderr," -Z [columns]     Sort by columns, in descending order"
    print >>sys.stderr," -v               Pretty view output, replaces tabs with spaces to nicely align columns. You may want to combine this with -n and --nl, and perhaps -N"
    print >>sys.stderr," -V               Pretty view output in a GUI"    
    print >>sys.stderr," --copysuffix=[suffix]       Output an output file with specified suffix for each inputfile (use instead of -o or -i)"
    print >>sys.stderr," --nl             Insert an extra empty newline after each line"
    print >>sys.stderr," --html           Output HTML table"
    print >>sys.stderr," --latex          Output LaTeX tabular"
    print >>sys.stderr,"Selection shortcuts:"
    print >>sys.stderr," -g [key]         Does a grep. Shortcut for: -s 'A() == \"key\"'"
    print >>sys.stderr," -G [key]         Does an inverse grep. Shortcut for: -s 'not (A() == \"key\"')"        
    print >>sys.stderr,"Options to be implemented still:"
    print >>sys.stderr," -X [samplesizes] Draw one or more random samples (non overlapping, comma separated list of sample sizes)"
    print >>sys.stderr," -a [column]=[columname]=[expression]   Adds a new column after the specified column"
    print >>sys.stderr," -J [sourcekey]:[filename]:[targetkey]:[selecttargetcolumns]:[insertafter]   Joins another data set with this one, on a specified column"  
    print >>sys.stderr,"Column specification:"
    print >>sys.stderr," A comma separated list of column index numbers or column names (if -1 option is used). Column index numbers start with 1. Negative numbers may be used for end-aligned-indices, where -1 is the last column. Ranges may be specified using a colon, for instance: 3:6 equals 3,4,5,6. A selection like 3:-1 select the third up to the last column. A specification like ID,NAME selects the columns names as such."
    print >>sys.stderr,"Selector specification:"
    print >>sys.stderr," The selection specification (-s) is normal python code and thus allows for a great deal of flexibility. You can use the normal boolean operators and, or, not to combine expressions. The following campyon-specific functions are available in selector context:"
    print >>sys.stderr,"    c(n)           Returns the value in column with index n"
    print >>sys.stderr,"    c('NAME')      Return the value in the column with the specified name"
    print >>sys.stderr,"    C((n,n....))   Match conjunction of multiple columns. An expression like   C((1,2)) > 4  is the same as: c(1) > 4 and c(2) > 4 . Names instead of numbers are also allowed. Note the double parentheses."
    print >>sys.stderr,"    D((n,n....))   Match disjunction of multiple columns. An expression like   D((1,2)) > 4  is the same as: c(1) > 4 or c(2) > 4 . Names instead of numbers are also allowed. Note the double parentheses."
    print >>sys.stderr,"    A()            Any field, applies a disjunction over all columns."
    print >>sys.stderr,"    r('REGEXP',n)  Apply regular expression on column n, uses standard python re.search(). Also:"
    print >>sys.stderr,"    r('REGEXP','NAME')"     
    print >>sys.stderr," The selection specification (-s) is normal python code and thus allows for a great deal of flexibility. The following functions are available in this context:"    
    print >>sys.stderr,"Plot options:"
    print >>sys.stderr," --plotgrid       Draw grid"
    print >>sys.stderr," --plotxlog       X scale is logarithmic"
    print >>sys.stderr," --plotylog       Y scale is logarithmic"
    print >>sys.stderr," --plotfile=[filename]      Save plot to PNG file"
    print >>sys.stderr," --lineplot       Sets lineplot defaults"
    print >>sys.stderr," --scatterplot    Sets scatterplot defaults"
            
def bold(s):
    CSI="\x1B["
    return CSI+"1m" + s + CSI + "0m"

def red(s):
    CSI="\x1B["
    return CSI+"31m" + s + CSI + "0m"  

def green(s):
    CSI="\x1B["
    return CSI+"32m" + s + CSI + "0m"   


def magenta(s):
    CSI="\x1B["
    return CSI+"35m" + s + CSI + "0m"   

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

class CampyonViewer(object):

    # close the window and quit
    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False

    def __init__(self, c, filename):
        # Create a new window
        self.window = gtk.Dialog() #gtk.WINDOW_TOPLEVEL)

        self.window.set_title(filename)
        self.window.set_border_width(0)
        self.window.set_size_request(640, 480)
        
        #self.window.set_size_request(640, 480)

        self.window.connect("delete_event", self.delete_event)
        
        self.scrollwindow = gtk.ScrolledWindow()
        self.scrollwindow.set_border_width(0)
        self.scrollwindow.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        self.scrollwindow.show()

        types = []
        first = True
        for line, fields, linenum in c.processmemory():  
            if first:
                for field in fields:      
                    if isinstance(field, float):
                        types.append(float)
                    elif isinstance(field, int):
                        types.append(int)
                    else:                        
                        types.append(str)                        
                first = False
            else:
                for i, field in enumerate(fields):
                    if types[i] == int and isinstance(field, float):
                        types[i] = float
                    elif isinstance(field, str) or isinstance(field, unicode):
                        types[i] = str
        if c.numberlines:
            types.insert(0,int)
                        
                
        print repr(types)
                
        # create a liststore with one string column to use as the model
        self.liststore = gtk.ListStore(*types) #(str, str, str, 'gboolean')

        # create the TreeView using liststore
        self.treeview = gtk.TreeView(self.liststore)

        # create the TreeViewColumns to display the data
        self.columns = []
        self.cellrenderers = []
        if c.numberlines:
                self.columns.append( gtk.TreeViewColumn('#') )
                self.cellrenderers.append( gtk.CellRendererText() )
                
        if c.DOHEADER:
            for colname in c.headerfields():
                self.columns.append( gtk.TreeViewColumn(colname) )
                self.cellrenderers.append( gtk.CellRendererText() ) 
        else:
            for num in range(1,len(types)+1):
                self.columns.append( gtk.TreeViewColumn(str(num)) )
                self.cellrenderers.append( gtk.CellRendererText() ) 

        #add data
        first = True
        for line, fields, linenum in c.processmemory(): 
            if c.numberlines: 
                fields = [linenum] + fields
            if (c.DOHEADER and not first) or not c.DOHEADER:             
                self.liststore.append(fields)
            first = False

        # add columns to treeview
        for col in self.columns:
            self.treeview.append_column(col)

        # create a CellRenderers to render the data
        
     
        

        # set background color property
        #self.cellpb.set_property('cell-background', 'yellow')
        for i,cr in enumerate(self.cellrenderers):
            if c.numberlines:
                fieldnum = i
            else:
                fieldnum = i + 1
            if fieldnum in c.highlight:
                cr.set_property('cell-background', 'yellow')
            else:
                cr.set_property('cell-background', 'white')
        #self.cell1.set_property('cell-background', 'pink')


        

            

        # set the cell attributes to the appropriate liststore column
        for i, (col, cr) in enumerate(zip(self.columns, self.cellrenderers)):
            col.pack_start(cr, True)         # add the cells to the columns 
            col.set_attributes(cr,text=i)
            col.set_sort_column_id(i)

        # make treeview searchable
        self.treeview.set_search_column(0)

        # Allow sorting on the column
        #self.tvcolumn.set_sort_column_id(0)

        # Allow drag and drop reordering of rows
        self.treeview.set_reorderable(True)

        self.scrollwindow.add_with_viewport(self.treeview)

        #self.window.add(self.treeview)
        
        self.window.vbox.pack_start(self.scrollwindow, True, True, 0)
        
        self.scrollwindow.show_all()
        self.window.show_all()


class Campyon(object):
    def _parsekwargs(self, key, default, kwargs):
        if key in kwargs:
            return kwargs[key]
        else:
            return default
    
    def parsecolumnindex(self, s):
        if isinstance(s,int):
            return s
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
	        opts, args = getopt.getopt(args, "f:k:d:e:D:o:is:SH:TC:nNM:1x:y:A:Z:a:vVg:G:",["bar","plotgrid","plotxlog","plotylog","plotconf=","plotfile=","scatterplot","lineplot","plottitle","copysuffix=","nl","html","latex"])
        except getopt.GetoptError, err:
	        # print help information and exit:
	        print str(err)
	        usage()
	        sys.exit(2)           
        
        self.filenames = self._parsekwargs('filenames',"",kwargs)
        self.encoding = self._parsekwargs('encoding',"utf-8",kwargs)
        self.delete = self._parsekwargs('delete',[],kwargs)
        self.keep = self._parsekwargs('keep',[],kwargs)
        self.delimiter = self._parsekwargs('delimiter',"",kwargs)
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

        self.copysuffix = self._parsekwargs('copysuffix',"",kwargs)
        
        self.plotgrid = self._parsekwargs('plotgrid',False,kwargs)
        self.plotxlog = self._parsekwargs('plotxlog',False,kwargs)
        self.plotylog = self._parsekwargs('plotylog',False,kwargs)
        self.plotconf = self._parsekwargs('plotconf',['r.-','g.-','b.-','y.-','m.-','c.-'],kwargs)
        self.plotfile = self._parsekwargs('plotfile',"",kwargs)
        self.plottitle = self._parsekwargs('plottitle',"",kwargs)
        
        self.prettyview = False
        self.extranewline = False
        self.guiview = False
        self.html = False
        self.latex = False
        
        self.fieldcount = 0
        self.header =  {}
        self.sortreverse = False
        self.inmemory = False
        self.xs = []
        self.ys = {}        
        
        self.keepsettings = ""
        self.deletesettings = ""
        self.histsettings = ""
        self.highlightsettings = ""
        self.sortsettings = ""
        self.plotxsettings = ""
        self.plotysettings = ""
    
        for o, a in opts:
            if o == "-e":	
                self.encoding = a
            elif o == "-f":	
                self.filenames = [a]              
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
            elif o == '-A':
                self.sortsettings = a
            elif o == '-Z':
                self.sortsettings = a
                self.sortreverse = True     
            elif o == '-x':           
                self.plotxsettings = a                
            elif o == '-y':           
                self.plotysettings = a
            elif o == '-v':           
                self.prettyview = True
            elif o == '-V':
                self.guiview = True
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
            elif o == "--lineplot":
                self.plotconf = self._parsekwargs('plotconf',['r-','g-','b-','y-','m-','c-'],kwargs)
            elif o == "--scatterplot":
                self.plotconf = self._parsekwargs('plotconf',['ro ','go ','bo ','yo ','mo ','co '],kwargs)
            elif o == '--copysuffix':
                self.copysuffix = a
            elif o == '--nl':
                self.extranewline = True
            elif o == '-g':
                self.select = 'A() == "' + a.replace('"','\\"') + '"'
            elif o == '-G':
                self.select = 'not (A() == "' + a.replace('"','\\"') + '")'                
            elif o == '-a':
                raise NotImplementedError
            else:
                raise Exception("invalid option: " + o)

        if args:
            self.filenames = args
                        
        if not self.filenames:    
            usage()
            sys.exit(2)

        for filename in self.filenames:
            if not os.path.exists(filename):
                print >>sys.stderr,"No such file: " + filename
                sys.exit(2)
 
           
        if self.sort or self.sortsettings or self.prettyview or self.guiview:
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
        
        
    def init(self, filename):       
        f = codecs.open(filename,'r',self.encoding)
        for line in f:
            if line.strip() and (not self.commentchar or line[:len(self.commentchar)] != self.commentchar):                    
                if not self.delimiter:
                    if "\t" in line:
                        self.delimiter = "\t"
                        print >>sys.stderr,"Guessed delimiter: TAB"
                    elif ";" in line:       
                        self.delimiter = ";"
                        print >>sys.stderr,"Guessed delimiter: SEMICOLON"
                    elif ":" in line:       
                        self.delimiter = ":"
                        print >>sys.stderr,"Guessed delimiter: COLON"                                                
                    elif "," in line:       
                        self.delimiter = ","
                        print >>sys.stderr,"Guessed delimiter: COMMA"
                    elif " " in line:                        
                        self.delimiter = " " 
                        print >>sys.stderr,"Guessed delimiter: SPACE"
                                           
                fields = line.strip().split(self.delimiter)
                self.fieldcount = len(fields)
                print >>sys.stderr,"Number of fields: ", self.fieldcount
                if self.DOHEADER:
                    self.header = dict([ (x+1,y.strip()) for x,y in enumerate(fields) ])
                    for col, name in self.header.items():
                        print >>sys.stderr,"Column #"+str(col)+":", name.encode('utf-8')
                break            
        f.close()
        
        
        if self.keepsettings: self.keep = self.parsecolumns(self.keepsettings)
        if self.deletesettings: self.delete = self.parsecolumns(self.deletesettings)
        if self.histsettings: self.hist = self.parsecolumns(self.histsettings)
        if self.highlightsettings: self.highlight = self.parsecolumns(self.highlightsettings)
        if self.sortsettings: self.sort = self.parsecolumns(self.sortsettings)
        if self.plotxsettings: self.x = self.parsecolumnindex(self.plotxsettings)
        if self.plotysettings: self.y = self.parsecolumns(self.plotysettings)                   
        
    def __call__(self):        
        self.memory = []    
        self.sumdata = {}
        self.nostats = set()
        self.freq = {}
        self.rowcount_in = 0
        self.rowcount_out = 0
        f_out = None
        
        if self.outputfile and not self.overwriteinput:
            f_out = codecs.open(self.outputfile, 'w',self.encoding)
                    
        if not self.overwriteinput and not self.copysuffix:                    
            self.init(self.filenames[0]) #initialise one, assume same column config for all!

        for filename in self.filenames:

            if self.overwriteinput or self.copysuffix:                
                self.memory = []
                self.nostats = set()                
                self.rowcount_in = 0
                self.rowcount_out = 0
                self.init(filename)
                if self.overwriteinput:
                  f_out = codecs.open(filename+".tmp", 'w',self.encoding)
                elif self.copysuffix:
                  f_out = codecs.open(filename + '.' + self.copysuffix, 'w',self.encoding)
            
            for line, fields, linenum in self.process(filename):
                if f_out:                                          
                    if self.numberlines: f_out.write(str(linenum) + self.delimiter)
                    f_out.write(line + "\n")
                else:
                    if self.numberlines: print green(str(linenum)) + self.delimiter,
                    print line.encode(self.encoding)
        
            if self.overwriteinput or self.copysuffix:             
                if self.prettyview:
                    margin = 2
                    colsize = {}
                    for line, fields, linenum in self.processmemory():                    
                        for i,field in enumerate(fields):
                            if not i in colsize:
                                colsize[i] = 0
                            if len(unicode(field))+margin > colsize[i]:
                                colsize[i] = len(unicode(field))+margin
                    for line, fields, linenum in self.processmemory():                                                        
                        for field in fields:
                            spaces = " " * (colsize[i] - len(unicode(field)))
                            if f_out:
                                f_out.write(unicode(field) + spaces)
                            else:
                                print unicode(field).encode(self.encoding) + spaces,
                        if f_out:
                            f_out.write("\n")
                        else:
                            print
                elif self.guiview: 
                    v = CampyonViewer(self, filename)
                    gtk.main()
                    del v
                    
                                     
            if f_out and (self.overwriteinput or self.copysuffix):
                if self.inmemory and not self.prettyview and not self.guiview:                    
                    for line, fields, linenum in self.processmemory():                                          
                        if self.numberlines: f_out.write(str(linenum) + self.delimiter)
                        f_out.write(line + "\n")                           
                f_out.close()
                f_out = None
                if self.overwriteinput:
                    os.rename(filename+".tmp",filename)
                
        if self.inmemory and not self.overwriteinput and not self.copysuffix:            
            
            if self.prettyview:
                margin = 2
                colsize = {}
                for line, fields, linenum in self.processmemory():                    
                    for i,field in enumerate(fields):
                        if not i in colsize:
                            colsize[i] = 0
                        if len(unicode(field))+margin > colsize[i]:
                            colsize[i] = len(unicode(field))+margin
                for line, fields, linenum in self.processmemory():                                                        
                    for field in fields:
                        spaces = " " * (colsize[i] - len(unicode(field)))
                        if f_out:
                            f_out.write(unicode(field) + spaces)
                        else:
                            print unicode(field).encode(self.encoding) + spaces,
                    if f_out:
                        f_out.write("\n")
                    else:
                        print
            elif self.guiview: 
                v = CampyonViewer(self, filename)
                gtk.main()
                del v
            else:            
                for line, fields, linenum in self.processmemory():
                    if f_out:                                          
                        if self.numberlines: f_out.write(str(linenum) + self.delimiter)
                        f_out.write(line + "\n")
                    else:
                        if self.numberlines: print green(str(linenum)) + self.delimiter,
                        print line.encode(self.encoding)            
              
        if f_out:
            f_out.close()        
            
        if self.DOSTATS:            
            self.printstats()
            
        if self.hist:
            for fieldnum in sorted(self.freq):
                print >>sys.stderr, "Histogram for column #" + str(fieldnum) + "\ttypes=" + str(self.types(fieldnum)) + "\ttokens=" + str(self.tokens(fieldnum)) + "\tttr=" +  str(self.ttr(fieldnum)) + "\tentropy=" + str(self.entropy(fieldnum))
                print >>sys.stderr,"------------------------------------------------------------------------"
                self.printhist(fieldnum)
                
        if self.x and self.y:
            self.plot()
            
    def __iter__(self):            
        self.memory = []    
        self.sumdata = {}
        self.nostats = set()
        self.freq = {}
        self.rowcount_in = 0
        self.rowcount_out = 0
                
        self.init(self.filenames[0]) #initialise one, assume same column config for all! 
                    
        for filename in self.filenames:        
            for line, fields, linenum in self.process(filename):
                yield line, fields, linenum
                 
    def __len__(self):        
        return self.rowcount_out
        
        
        
    def process(self, f):                                    
        if self.overwriteinput:
            self.rowcount_in = 0
            self.rowcount_out = 0 
        
        if self.keep: 
            default = 'delete'
        else:        
            default = 'keep'        
            
        headerfound = False

        if isinstance(f, str) or isinstance(f, unicode):                           
            f = codecs.open(f,'r',self.encoding)
            
        for line in f:
            if not isinstance(line, unicode):
                line = unicode(line, self.encoding)
            isheader = False
            self.rowcount_in += 1            
            
            if not line.strip() or (self.commentchar and line[:len(self.commentchar)] == self.commentchar):
                self.rowcount_out += 1
                if not self.inmemory:
                    yield line.strip(), [], self.rowcount_out
                continue
            
            
            

                
            fields = line.strip().split(self.delimiter)
            if len(fields) != self.fieldcount:
                raise CampyonError("Number of columns in line " + str(self.rowcount_in) + " deviates, expected " + str(self.fieldcount) + ", got " + str(len(fields))) 
                
            if self.DOHEADER and not headerfound:
                headerfound = True
                isheader = True
    
        
            
            if self.select and not isheader:
                c = lambda x: fields[self.parsecolumnindex(x)-1].strip()
                C = lambda x: ConjunctionSelector(c, *x)
                D = lambda x: DisjunctionSelector(c, *x)  
                r = lambda x,y: re.search(x,c(y))
                A = lambda: D(range(1,len(fields)+1))
                if not eval(self.select):
                    continue
        
            self.rowcount_out += 1
            

            if self.hist and not isheader:
                for fieldnum in self.hist:
                    if not fieldnum in self.freq:                
                        self.freq[fieldnum] = {}
                    if not fields[fieldnum-1] in self.freq[fieldnum]:
                        self.freq[fieldnum][fields[fieldnum-1]] = 0
                    self.freq[fieldnum][fields[fieldnum-1]] += 1
                
            if self.DOSTATS and not isheader:                    
                for i,field in enumerate(fields):
                    fieldnum = i+1
                    if not fieldnum in self.nostats: 
                        if '.' in field:
                            try:
                               x = float(field)
                            except:
                               self.nostats.add(fieldnum)
                               if fieldnum in self.sumdata: del self.sumdata[fieldnum]
                               continue
                        else:
                            try:
                               x = int(field)
                            except:
                               self.nostats.add(fieldnum)
                               if fieldnum in self.sumdata: del self.sumdata[fieldnum]
                               continue                        
                     
                        if not fieldnum in self.sumdata:
                            self.sumdata[fieldnum] = 0
                        self.sumdata[fieldnum] += x
                         
                
            
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
                if fieldnum in self.highlight and not self.guiview:
                    field = bold(red(field))                
                if self.numberfields:
                    if self.guiview:
                        field = str(fieldnum) + '=' + field
                    else:
                        field = magenta(str(fieldnum)) + '=' + field
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

                if self.x == fieldnum and not isheader:
                    self.xs.append(field)
                
                if fieldnum in self.y and not isheader:
                    if not isinstance(field, float) and not isinstance(field,int):
                        raise CampyonError("Can not plot non-numeric values: " + field) 
                        
                    if not fieldnum in self.ys:
                        self.ys[fieldnum] = []
                    self.ys[fieldnum].append(field)

                if action == 'keep':
                    newfields.append(field)
                    
            s = self.delimiter.join([ unicode(x) for x in newfields ])            
            if self.inmemory:
                if not isheader:
                    self.memory.append( (newfields, self.rowcount_out) )
            else:                
                yield s, newfields, self.rowcount_out

            
        print >>sys.stderr,"Read " + str(self.rowcount_in) + " lines, outputted " + str(self.rowcount_out)

    def processmemory(self):    
        if self.sort:
           self.memory = sorted(self.memory, key=lambda x: tuple([ x[0][i-1] for i in self.sort ]), reverse=self.sortreverse)
                          
        if self.header:    
            s = self.delimiter.join( self.headerfields()  )
            yield s, self.headerfields(), 0

        for fields, linenum in self.memory:
            s = self.delimiter.join([ unicode(x) for x in fields])
            yield s, fields, linenum

        

        
    def plot(self, show=True):        
        barcolors = 'rgbymc'
        
        #fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.clf()
        if self.plotgrid:
            matplotlib.pyplot.grid(True)
        else:
            matplotlib.pyplot.grid(False)            

        if self.plottitle:
            matplotlib.pyplot.title(self.plottitle)
            
        if all([ isinstance(x,float) or isinstance(x,int) for x in self.xs ]):
            if self.plotylog:            
                matplotlib.pyplot.set_yscale('log')            
            
            if self.x in self.header:            
                matplotlib.pyplot.xlabel(self.header[self.x])            
            if len(self.y) == 1 and self.y[0] in self.header:            
                matplotlib.pyplot.ylabel(self.header[self.y[0]])            
            
            if self.plotxlog:
                matplotlib.pyplot.set_xscale('log')

            l = []
            for i, field in enumerate(self.y):
                l.append(self.xs)
                l.append(self.ys[field])
                l.append(self.plotconf[i])
                
                            
            matplotlib.pyplot.plot(*l)
        else:            
           #do a horizontal barplot    
           
           if self.plotylog:            
                matplotlib.pyplot.set_xscale('log')
           
           if self.x in self.header:            
                matplotlib.pyplot.ylabel(self.header[self.x])            
           if len(self.y) == 1 and self.y[0] in self.header:            
                matplotlib.pyplot.xlabel(self.header[self.y[0]])   
            
             
           hbarheight = 0.2
           locations = numpy.arange(len(self.xs))
           for i, field in enumerate(self.y): #TODO support multiple graphs
                matplotlib.pyplot.barh(locations ,  self.ys[field], align='center', color=barcolors[i])
           matplotlib.pyplot.yticks(locations+hbarheight/2., self.xs)
           
            
        if self.plotfile:
            print >>sys.stderr, "Saving plot in " + self.plotfile
            matplotlib.pyplot.savefig(self.plotfile, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        elif show:
            print >>sys.stderr, "Showing plot"
            matplotlib.pyplot.show()
        
    def headerfields(self):          
        return [x[1] for x in sorted(self.header.items()) ]
        
    def printstats(self, out=sys.stderr):            
        out.write("COLUMN\tSUM\tAVERAGE\n")
        for colnum, colname, s, average in self.stats():
            if colname == str(colnum): 
                out.write(str(colnum) + "\t"+ str(s) + "\t" + str(average) + "\n")
            else:
                out.write(colname.encode(self.encoding) + "\t"+ str(s) + "\t" + str(average) + "\n")                
            
    def stats(self):
        for i in sorted(self.sumdata):
            if self.header:
                colname = self.header[i]
            else:
                colname = str(i)   
            
            yield i, colname, self.sumdata[i], self.sumdata[i] / float(self.rowcount_out)
                         
     
    def printhist(self, columnindex, out=sys.stderr):         
        for i, (word, count, f) in enumerate(self.histdata(columnindex)):
            print >>sys.stderr, str(i+1) + ")\t" + word.encode(self.encoding) + "\t" + str(count) + "\t" + str(f * 100) + '%'
        
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
        for i,key in self.header.items():
            if key == colname:
                return i
        raise KeyError("Column " + colname + " not found")


class ConjunctionSelector(object):
    def __init__(self, c, *args):    
        self.args = [ c(x) for x in args ]

    def __eq__(self, y):
        return all([ x == y for x in self.args ])
    
    def __ne__(self, y):
        return all([ x != y for x in self.args ])    
    
    def __gt__(self, y):
        return all([ x > y for x in self.args ])

    def __lt__(self, y):
        return all([ x < y for x in self.args ])
    
    def __ge__(self, y):
        return all([ x >= y for x in self.args ])
    
    def __le__(self, y):
        return all([ x <= y for x in self.args ])    
    
    def __contains__(self, y):
        return all([ x in y for x in self.args ])    
    
class DisjunctionSelector(object):
    def __init__(self, c, *args):
        self.args = [ c(x) for x in args ]    
    
    def __eq__(self, y):
        return any([ x == y for x in self.args ])
    
    def __ne__(self, y):
        return any([ x != y for x in self.args ])    
    
    def __gt__(self, y):
        return any([ x > y for x in self.args ])

    def __lt__(self, y):
        return any([ x < y for x in self.args ])
    
    def __ge__(self, y):
        return any([ x >= y for x in self.args ])
    
    def __le__(self, y):
        return any([ x <= y for x in self.args ])    

    def __contains__(self, y):
        return any([ x in y for x in self.args ])        
    

        
if __name__ == "__main__":         
    campyon = Campyon(*sys.argv[1:])
    campyon()
