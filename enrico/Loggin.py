import sys

class Message(object):
    def __init__(self):
        self.classname = self.__class__.__name__
        self.errorcolor = "\033[31m"#red
        self.infocolor = "\033[34m"#blue
        self.warningcolor = "\033[33m"#yellow
        self.successcolor = "\033[32m"#green
        self.endcolor = "\033[0m"#reset
        
    def error(self,message, functionname = ""):
        printstring = ""
        if functionname == "":
            printstring = "\n"+self.errorcolor+"*** Error ["+self.classname+"]: "+message+" ***\n"+self.endcolor
        else:
            printstring = "\n"+self.errorcolor+"*** Error ["+self.classname+"::"+functionname+"]: "+message+" ***\n"+self.endcolor
        sys.exit(printstring)
    
    def info(self,message,newline=True):
        printstring = self.infocolor+"["+self.classname+"]: "+message+self.endcolor     
        if newline:
            print printstring
        else:
            print self.infocolor+message+self.endcolor,
            sys.stdout.flush()  
        
    def warning(self,message,functionname = ""):
        printstring = ""
        if functionname == "":
            printstring = self.warningcolor+"["+self.classname+"] Warning: "+message+self.endcolor
        else:
            printstring = self.warningcolor+"["+self.classname+"::"+functionname+"] Warning: "+message+self.endcolor
        print printstring
        
    def success(self,message):
        printstring = self.successcolor+"["+self.classname+"]: "+message+self.endcolor
        print printstring
        
