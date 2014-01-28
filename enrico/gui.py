import gtk
import sys,os
import logging
from math import log10,pow
from enrico.extern.configobj import ConfigObj, flatten_errors
from enrico.config import get_config

class EnricoGui:

    def Launch(self,widget,data=None):
       self.save(None,None)
       os.system(data+" "+self.infile)

    def LaunchEbin(self,widget,data=None):
       self.save(None,None)
       if self.config["Ebin"]["NumEnergyBins"]>0 :
           os.system("enrico_sed "+self.config["out"]+'/Ebin'+str(self.config["Ebin"]["NumEnergyBins"])+"/*conf" )

    def Sync(self, widget, data=None):
       self.x.set_value(self.ra.get_value())
       self.y.set_value(self.dec.get_value())

    def set_active(self, widget, data):
       widget.set_active(0)
       if data == 'yes':
           widget.set_active(1)

    def fct_yesno(self, widget, data=None):
       if data=="Spectrum":
         self.config["Spectrum"]["FitsGeneration"] = "no"
         if widget.get_active():
           self.config["Spectrum"]["FitsGeneration"] = "yes"

       elif data=="LightCurve":
         self.config["LightCurve"]["FitsGeneration"] = "no"
         if widget.get_active():
           self.config["LightCurve"]["FitsGeneration"] = "yes"

       elif data=="AppLC":
         self.config["AppLC"]["FitsGeneration"] = "no"
         if widget.get_active():
           self.config["AppLC"]["FitsGeneration"] = "yes"

       elif data=="Ebin":
         self.config["Ebin"]["FitsGeneration"] = "no"
         if widget.get_active():
           self.config["Ebin"]["FitsGeneration"] = "yes"

       elif data=="findsrc":
         self.config["findsrc"]["FitsGeneration"] = "no"
         if widget.get_active():
           self.config["findsrc"]["FitsGeneration"] = "yes"

       elif data=="verbose":
         self.config["verbose"] = "no"
         if widget.get_active():
           self.config["verbose"] = "yes"

       elif data=="clobber":
         self.config["clobber"] = "no"
         if widget.get_active():
           self.config["clobber"] = "yes"

       elif data=="Submit":
         self.config["Submit"] = "no"
         if widget.get_active():
           self.config["Submit"] = "yes"

       elif data=="refit":
         self.config["findsrc"]["Refit"] = "no"
         if widget.get_active():
           self.config["findsrc"]["Refit"] = "yes"

       elif data=="ULenvelope":
         self.config["UpperLimit"]["envelope"] = "no"
         if widget.get_active():
           self.config["UpperLimit"]["envelope"] = "yes"

       elif data=="conffile":
         self.config["LightCurve"]["MakeConfFile"] = "no"
         if widget.get_active():
           self.config["LightCurve"]["MakeConfFile"] = "yes"

       elif data=="compvar":
         self.config["LightCurve"]["ComputeVarIndex"] = "no"
         if widget.get_active():
           self.config["LightCurve"]["ComputeVarIndex"] = "yes"

       elif data=="lcdiagplot":
         self.config["LightCurve"]["DiagnosticPlots"] = "no"
         if widget.get_active():
           self.config["LightCurve"]["DiagnosticPlots"] = "yes"

       elif data=="binfromdata":
         self.config["AppLC"]["binsFromData"] = "no"
         if widget.get_active():
           self.config["AppLC"]["binsFromData"] = "yes"

       elif data=="tsmap":
         self.config["TSMap"]["Re-Fit"] = "no"
         if widget.get_active():
           self.config["TSMap"]["Re-Fit"] = "yes"

       elif data=="removetgr":
         self.config["TSMap"]["RemoveTarget"] = "no"
         if widget.get_active():
           self.config["TSMap"]["RemoveTarget"] = "yes"

       elif data=="Diffresp":
         self.config["analysis"]["ComputeDiffrsp"] = "no"
         if widget.get_active():
           self.config["analysis"]["ComputeDiffrsp"] = "yes"

       elif data=="roicut":
         self.config["analysis"]["roicut"] = "no"
         if widget.get_active():
           self.config["analysis"]["roicut"] = "yes"


    def fct_rappel(self, widget, data=None):
       print "Le %s a ete %s." % (data, ("desactive", "active")[widget.get_active()])

    def delete(self, widget, event=None):
        gtk.main_quit()
        return False

    def change(self,widget,data=None):
        try :
          pixbuf = gtk.gdk.pixbuf_new_from_file(data)
          pixbuf = pixbuf.scale_simple(600, 300, gtk.gdk.INTERP_BILINEAR)
          self.image.set_from_pixbuf(pixbuf)

        except:
          logging.error('picture file '+data+' not found.')

    def save(self, widget, data=None):

        self.config['out'] = self.fout.get_current_folder()

        self.config["UpperLimit"]["SpectralIndex"] = self.ULindex.get_value()
        self.config["Spectrum"]["FrozenSpectralIndex"] = self.index.get_value()
        self.config["Spectrum"]["cl"] = self.ULcl.get_value()

        self.config["file"]["xml"] = self.fxml.get_filename()
        self.config["file"]["spacecraft"] = self.fsc.get_filename()
        self.config["file"]["event"] = self.fevent.get_filename()
        self.config["file"]["tag"] = self.ftag.get_text()


        self.config["target"]["name"] = self.fname.get_text()
        self.config["target"]["ra"] = self.ra.get_value()
        self.config["target"]["dec"] = self.dec.get_value()
        self.config["space"]["spectrum"] = self.listSpec.entry.get_text()

        self.config["space"]["xref"] = self.x.get_value()
        self.config["space"]["yref"] = self.y.get_value()
        self.config["space"]["phibins"] = int(self.phibin.get_value())

        self.config["space"]["rad"] = self.rad.get_value()
        self.config["space"]["binsz"] = self.binsz.get_value()

        self.config["space"]["coordsys"] = self.listSys.entry.get_text()
        self.config["space"]["proj"] = self.listProj.entry.get_text()

        self.config["energy"]["emin"] = self.emin.get_value()
        self.config["energy"]["emax"] = self.emax.get_value()
        self.config["energy"]["enumbins_per_decade"] = int(self.nbindec.get_value())

        self.config["time"]["tmin"] = self.tmin.get_value()
        self.config["time"]["tmax"] = self.tmax.get_value()
        if self.ftime.get_filename()==None:
           self.config['time']['file'] = ""
        else: 
           self.config['time']['file'] = self.ftime.get_filename()


        self.config["time"]["type"] = self.listtime.entry.get_text()

        self.config["Ebin"]["NumEnergyBins"] = int(self.nebin.get_value())
        self.config["Ebin"]["TSEnergyBins"] = int(self.tsebin.get_value())

        self.config["LightCurve"]["NLCbin"] = int(self.nlcbin.get_value())
        self.config["LightCurve"]["index"] = self.lcindex.get_value()

        self.config["AppLC"]["index"] = self.applcindex.get_value()
        self.config["AppLC"]["NLCbin"] = int(self.applcNbin.get_value())

        self.config["FoldedLC"]["NLCbin"] =int(self.follcNbin.get_value())
        self.config["FoldedLC"]["epoch"] =(self.folepoch.get_value())
        self.config["FoldedLC"]["Period"] =(self.folperiod.get_value())

        self.config["TSMap"]["npix"] = int(self.tsmapnpix.get_value())
        self.config["TSMap"]["method"] = self.listtsmethod.entry.get_text()

        self.config["analysis"]["likelihood"] = self.listchain.entry.get_text()
        self.config["analysis"]["evclass"] = int(self.evclass.get_value())
        self.config["analysis"]["zmax"] = self.zmax.get_value()
        self.config["analysis"]["filter"] = self.filter.get_text()
        self.config["analysis"]["irfs"] = self.irfs.get_text()
        self.config["analysis"]["convtype"] = self.convtype.get_text()

        self.config["fitting"]["optimizer"] = self.listopt.entry.get_text()
        self.config["fitting"]["ftol"] = pow(10,self.ftol.get_value())

        self.config.write(open(self.infile, 'w'))


#    def load(self, widget, data=None):
#        print "loading file ",self.fc.get_filename()
#        self.config = get_config(self.fc.get_filename()) 

    def AddBlocNotePage(self,lab='bn'):
        frame = gtk.Table(6, 5, True)
        frame.set_size_request(600, 400)
        frame.show()

        label = gtk.Label(lab)
        self.notebloc.append_page(frame, label)
        return frame

    def AddSpinButton(self,val,mini,maxi,incre,decimal):
        button = gtk.SpinButton(gtk.Adjustment(val, mini, maxi, incre), .5,decimal)
        button.set_numeric(True)
        button.show()
        return button

#    def _addconfigBN(self):
#        BNPage = self.AddBlocNotePage("Config file")
#        self.fc = gtk.FileChooserButton("Config file")
#        self.fc.set_title("config file")
#        self.fc.set_size_request(600, 400)
#        self.fc.show()
#        BNPage.attach(self.fc, 4, 8, 0, 1)

#        button = gtk.Button("Load config File")
#        button.connect("clicked", self.load, None)
#        button.show()
#        BNPage.attach(button,  5, 8, 6, 8)

    def _addFrame(self,title,x,y):
        frame = gtk.Frame(title)
        table = gtk.Table(x, y, True)
        table.show()
        frame.add(table)
        frame.show() 
        return frame,table

    def _addTargetSpace(self):
        BNPage = self.AddBlocNotePage("Target/Space")
        frameTarget,tableTarget = self._addFrame("Target",2,3)
        BNPage.attach(frameTarget, 0, 8, 0, 2)

        self.ra = self.AddSpinButton(self.config["target"]["ra"], 0, 360, .01,2)
        self.dec = self.AddSpinButton(self.config["target"]["dec"], -90, 90, .01,2)

        self.fname = gtk.Entry()
        self.fname.set_text(self.config['target']['name'])
        self.fname.show()

        label = gtk.Label("Target name")
        label.show()
        tableTarget.attach(label, 0,1,0,1)
        tableTarget.attach(self.fname, 1,2,0,1)

        label = gtk.Label("RA")
        label.show()
        tableTarget.attach(label, 0,1,1, 2)
        tableTarget.attach(self.ra,1,2,1, 2)

        label = gtk.Label("Dec")
        label.show()
        tableTarget.attach(label,2, 3,1, 2)
        tableTarget.attach(self.dec,3,4,1, 2)


        self.listSpec = gtk.Combo()
        self.listSpec.entry.set_text("Spectrum")
        listoption = [ "PowerLaw", "PowerLaw2", "LogParabola", "PLExpCutoff", "Generic" ]
        listoption.remove(self.config['target']['spectrum'])
        listoption.insert(0,self.config['target']['spectrum'])
        self.listSpec.set_popdown_strings(listoption)
        self.listSpec.show()

        label = gtk.Label("Model")
        label.show()
        tableTarget.attach(label, 0,1,2,3)
        tableTarget.attach(self.listSpec, 1,2,2,3)

        frameSpace,tableSpace = self._addFrame("Space",2,3)
        BNPage.attach(frameSpace, 0, 8, 2, 6)

        self.x = self.AddSpinButton(self.config["space"]["xref"], 0, 360, .01,2)
        self.y = self.AddSpinButton(self.config["space"]["yref"], -90, 90, .01,2)

        label = gtk.Label("Xref")
        label.show()
        tableSpace.attach(label, 0,1,1, 2)
        tableSpace.attach(self.x,1,2,1, 2)

        label = gtk.Label("Yref")
        label.show()
        tableSpace.attach(label,2, 3,1, 2)
        tableSpace.attach(self.y,3,4,1,2)

        syncButton = gtk.Button("sync with RA-Dec")
#        syncButton.set_size_request(20,20)
        syncButton.connect("clicked", self.Sync, "")
        syncButton.show()
        tableSpace.attach(syncButton,2,4,2,3)


        self.rad = self.AddSpinButton(self.config["space"]["rad"],0, 360, .1, 1)
        self.binsz = self.AddSpinButton(self.config["space"]["binsz"],  0, 1, .01, 2)
        self.phibin = self.AddSpinButton(self.config["space"]["phibins"], 0, 40, 1, 0)

        label = gtk.Label("ROI")
        label.show()
        tableSpace.attach(label, 0,1,0, 1)
        tableSpace.attach(self.rad,1,2,0,1)

        label = gtk.Label("Bin size")
        label.show()
        tableSpace.attach(label,2, 3,0,1)
        tableSpace.attach(self.binsz,3,4,0,1)

        label = gtk.Label("Number of bin in phi (default 0)")
        label.show()
        tableSpace.attach(label,0, 3,4,5)
        tableSpace.attach(self.phibin,3,4,4,5)

        self.listProj = gtk.Combo()
        self.listProj.entry.set_text("Projection")
        listoption = [ "AIT","ARC","CAR","GLS","MER","NCP","SIN","STG","TAN" ]
        listoption.remove(self.config['space']['proj'])
        listoption.insert(0,self.config['space']['proj'])
        self.listProj.set_popdown_strings(listoption)
        self.listProj.show()

        label = gtk.Label("Projection")
        label.show()
        tableSpace.attach(label,0,1,3,4)
        tableSpace.attach(self.listProj,1,2,3,4)

        self.listSys = gtk.Combo()
        self.listSys.entry.set_text("Systeme")
        listoption = [ "CEL", "GAL" ]
        listoption.remove(self.config['space']['coordsys'])
        listoption.insert(0,self.config['space']['coordsys'])
        self.listSys.set_popdown_strings(listoption)
        self.listSys.show()

        label = gtk.Label("Systeme")
        label.show()
        tableSpace.attach(label,2,3,3,4)
        tableSpace.attach(self.listSys,3,4,3,4)
        

    def _addMainOption(self):
        BNPage = self.AddBlocNotePage("Main options")
        frame,table = self._addFrame("General options",2,3)

        self.fout = gtk.FileChooserButton("Config file")
        self.fout.set_current_folder(self.config['out'])
        self.fout.set_action(gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER)
        self.fout.show()
#        self.fout.set_visible(True)

        VerboseButton = gtk.CheckButton("verbose")
        VerboseButton.connect("toggled", self.fct_yesno, "verbose")
        self.set_active(VerboseButton,self.config["verbose"])
        VerboseButton.show()

        ClobberButton = gtk.CheckButton("Clobber")
        ClobberButton.connect("toggled", self.fct_yesno, "clobber")
        self.set_active(ClobberButton,self.config["clobber"])
        ClobberButton.show()

        SubmitButton = gtk.CheckButton("Submit jobs")
        SubmitButton.connect("toggled", self.fct_yesno, "Submit")
        self.set_active(SubmitButton,self.config["Submit"])
        SubmitButton.show()


        label = gtk.Label("Output directory")
        label.show()
        table.attach(label, 0,1, 0, 1)
        table.attach(self.fout, 1,3, 0, 1,gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        table.attach(VerboseButton, 0,1, 1, 2)
        table.attach(ClobberButton,1,2,1, 2)
        table.attach(SubmitButton,2,3,1, 2)
        BNPage.attach(frame, 0, 8, 0, 2)

        framebut,tablebut = self._addFrame("Launcher",2,5)

        sedbutton = gtk.Button("run enrico_ sed")
        sedbutton.connect("clicked", self.Launch, 'enrico_sed')
        sedbutton.show()
        tablebut.attach(sedbutton,0, 2,0,1)

        lcbutton = gtk.Button("run enrico_ lc")
        lcbutton.connect("clicked", self.Launch, 'enrico_lc')
        lcbutton.show()
        tablebut.attach(lcbutton, 2,4,0,1)
 
        tsbutton = gtk.Button("run enrico_ tsmap")
        tsbutton.connect("clicked", self.Launch, 'enrico_tsmap')
        tsbutton.show()
        tablebut.attach(tsbutton, 4,6,0,1)


        sedplotbutton = gtk.Button("run enrico_ plot_ sed")
        sedplotbutton.connect("clicked", self.Launch, 'enrico_plot_sed')
        sedplotbutton.show()
        tablebut.attach(sedplotbutton,0, 2,1,2)

        lcplotbutton = gtk.Button("run enrico_ plot_ lc")
        lcplotbutton.connect("clicked", self.Launch, 'enrico_plot_lc')
        lcplotbutton.show()
        tablebut.attach(lcplotbutton, 2,4,1,2)
 
        tsplotbutton = gtk.Button("run enrico_ plot_ tsmap")
        tsplotbutton.connect("clicked", self.Launch, 'enrico_plot_tsmap')
        tsplotbutton.show()
        tablebut.attach(tsplotbutton, 4,6,1,2)

        sep2 = gtk.HSeparator()
        sep2.show()
        tablebut.attach(sep2, 0,6,2,3)

        xmlbutton = gtk.Button("run enrico_ xml")
        xmlbutton.connect("clicked", self.Launch, 'enrico_xml')
        xmlbutton.show()
        tablebut.attach(xmlbutton,0, 2,3,4)

        modelbutton = gtk.Button("run enrico_ testmodel")
        modelbutton.connect("clicked", self.Launch, 'enrico_testmodel')
        modelbutton.show()
        tablebut.attach(modelbutton, 2,4,3,4)
 
        findsrcbutton = gtk.Button("run enrico_ findsrc")
        findsrcbutton.connect("clicked", self.Launch, 'enrico_findsrc')
        findsrcbutton.show()
        tablebut.attach(findsrcbutton, 4,6,3,4)

        sep2 = gtk.HSeparator()
        sep2.show()
        tablebut.attach(sep2, 0,6,4,5)

        applcbutton = gtk.Button("run enrico_ applc")
        applcbutton.connect("clicked", self.Launch, 'enrico_applc')
        applcbutton.show()
        tablebut.attach(applcbutton,0, 2,5,6)

        foldedlcbutton = gtk.Button("run enrico_ foldedlc")
        foldedlcbutton.connect("clicked", self.Launch, 'enrico_foldedlc')
        foldedlcbutton.show()
        tablebut.attach(foldedlcbutton, 2,4,5,6)
 
        foldedlcplotbutton = gtk.Button("run enrico_ plot_ foldedlc")
        foldedlcplotbutton.connect("clicked", self.Launch, 'enrico_plot_foldedlc')
        foldedlcplotbutton.show()
        tablebut.attach(foldedlcplotbutton, 4,6,5,6)

        BNPage.attach(framebut, 0, 8, 2,6)
 
    def _addAnalysis(self):
        BNPage = self.AddBlocNotePage("Analysis")
        frameAna,tableAna = self._addFrame("Analysis options",2,4)
        frameFit,tableFit = self._addFrame("Fitting options",2,4)

        BNPage.attach(frameAna, 0, 8, 0, 4)
        self.listchain = gtk.Combo()
        self.listchain.entry.set_text("Chain")
        listoption = [ "unbinned","binned" ]
        listoption.remove(self.config['analysis']['likelihood'])
        listoption.insert(0,self.config['analysis']['likelihood'])
        self.listchain.set_popdown_strings(listoption)
        self.listchain.show()

        label = gtk.Label("Chain")
        label.show()
        tableAna.attach(label,0,1,0,1)
        tableAna.attach(self.listchain,1,2,0,1)

        self.evclass = self.AddSpinButton(self.config["analysis"]["evclass"], 1, 4, 1, 0)
        label = gtk.Label("Event Class")
        label.show()
        tableAna.attach(label, 2,3,0,1)
        tableAna.attach(self.evclass,3,4,0,1)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "Diffresp")
        self.set_active(FitButton,self.config["analysis"]["ComputeDiffrsp"])
        FitButton.show()

        label = gtk.Label("Compute diffuse response (only for UNBINNED)")
        label.show()
        tableAna.attach(label,0, 3,1,2)
        tableAna.attach(FitButton,3,4,1,2)

        self.zmax = self.AddSpinButton(self.config["analysis"]["zmax"], 50, 180, 1, 0)
        label = gtk.Label("zmax")
        label.show()
        tableAna.attach(label, 0,1,2,3)
        tableAna.attach(self.zmax,1,2,2,3)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "roicut")
        self.set_active(FitButton,self.config["analysis"]["roicut"])
        FitButton.show()

        label = gtk.Label("ROI cut")
        label.show()
        tableAna.attach(label,2, 3,2,3)
        tableAna.attach(FitButton,3,4,2,3)

        self.filter = gtk.Entry()
        self.filter.set_text(self.config['analysis']['filter'])
        self.filter.show()

        label = gtk.Label("filter")
        label.show()
        tableAna.attach(label, 0,1,3,4)
        tableAna.attach(self.filter, 1,4,3,4)

        self.irfs = gtk.Entry()
        self.irfs.set_text(self.config['analysis']['irfs'])
        self.irfs.show()

        label = gtk.Label("irfs")
        label.show()
        tableAna.attach(label, 0,1,4,5)
        tableAna.attach(self.irfs, 1,2,4,5)

        self.convtype = self.AddSpinButton(self.config["analysis"]["convtype"], -1,1, 1, 0)
        label = gtk.Label("convtype")
        label.show()
        tableAna.attach(label, 2,3,4,5)
        tableAna.attach(self.convtype,3,4,4,5)

        BNPage.attach(frameFit, 0, 8, 4, 6)

        self.listopt = gtk.Combo()
        self.listopt.entry.set_text("Optimizer")
        listoption = [ 'MINUIT', 'DRMNGB', 'DRMNFB', 'NEWMINUIT' ]
        listoption.remove(self.config['fitting']['optimizer'])
        listoption.insert(0,self.config['fitting']['optimizer'])
        self.listopt.set_popdown_strings(listoption)
        self.listopt.show()

        label = gtk.Label("Optimizer")
        label.show()
        tableFit.attach(label,0,1,0,1)
        tableFit.attach(self.listopt,1,2,0,1)

        self.ftol = self.AddSpinButton(log10(self.config["fitting"]["ftol"]), -9,-3, 1, 0)
        label = gtk.Label("log(tolerance)")
        label.show()
        tableFit.attach(label,0,1,1,2)
        tableFit.attach(self.ftol,1,2,1,2)

    def _addFiles(self):
        BNPage = self.AddBlocNotePage("Files")
        frame,table = self._addFrame("File options",4,2)

        BNPage.attach(frame, 0, 8, 0, 4)

        self.fevent = gtk.FileChooserButton("Event file")
        self.fevent.set_title("Event file")
        self.fevent.set_filename(self.config['file']['event'])
        self.fevent.show()

        self.fsc = gtk.FileChooserButton("Spacecraft file")
        self.fsc.set_title("Spacecraft file")
        self.fsc.set_filename(self.config['file']['spacecraft'])
        self.fsc.show()

        self.fxml = gtk.FileChooserButton("XML file")
        self.fxml.set_title("XML file")
        self.fxml.set_filename(self.config['file']['xml'])
        self.fxml.show()


        label = gtk.Label("Event file")
        label.show()
        table.attach(label, 0,1, 0, 1)
        table.attach(self.fevent, 1,3, 0, 1,gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        label = gtk.Label("spacecraft file")
        label.show()
        table.attach(label, 0,1, 1,2)
        table.attach(self.fsc, 1,3, 1,2,gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        label = gtk.Label("XML file")
        label.show()
        table.attach(label, 0,1, 2,3)
        table.attach(self.fxml, 1,3, 2,3,gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        self.ftag = gtk.Entry()
        self.ftag.set_text(self.config['file']['tag'])
        self.ftag.show()

        label = gtk.Label("tag")
        label.show()
        table.attach(label, 0,1,3,4)
        table.attach(self.ftag, 1,2,3,4)

    def _addEnergyTime(self):
        BNPage = self.AddBlocNotePage("Energy/Time")
        frameEner,tableEner = self._addFrame("Energy",2,2)
        frameTime,tableTime = self._addFrame("Time",3,2)

        BNPage.attach(frameEner, 0, 8, 0, 2)

        self.emin = self.AddSpinButton(self.config["energy"]["emin"], 0, 1e6, 100, 0)
        self.emax = self.AddSpinButton(self.config["energy"]["emax"], 0, 1e6, 100, 0)
        self.nbindec = self.AddSpinButton(self.config["energy"]["enumbins_per_decade"], 0, 40, 1, 0)

        label = gtk.Label("Emin (MeV)")
        label.show()
        tableEner.attach(label, 0,1,0, 1)
        tableEner.attach(self.emin,1,2,0,1)

        label = gtk.Label("Emax (MeV)")
        label.show()
        tableEner.attach(label,2, 3,0,1)
        tableEner.attach(self.emax,3,4,0,1)


        label = gtk.Label("Number of bin per decade (default 10)")
        label.show()
        tableEner.attach(label, 0,3,1, 2)
        tableEner.attach(self.nbindec,3,4,1,2)

        BNPage.attach(frameTime, 0, 8, 2, 5)

        self.tmin = self.AddSpinButton(self.config["time"]["tmin"], 239557418., 1e10, 1, 1)
        self.tmax = self.AddSpinButton(self.config["time"]["tmax"], 239557418., 1e10, 1, 1)

        label = gtk.Label("Tmin")
        label.show()
        tableTime.attach(label, 0,1,0, 1)
        tableTime.attach(self.tmin,1,2,0,1)

        label = gtk.Label("Tmax")
        label.show()
        tableTime.attach(label,2, 3,0,1)
        tableTime.attach(self.tmax,3,4,0,1)


        self.ftime = gtk.FileChooserButton("Time definition file")
        self.ftime.set_title("Time definition file")
        self.ftime.set_filename(self.config['time']['file'])
        self.ftime.show()

        label = gtk.Label("Time definition file (optinal)")
        label.show()
        tableTime.attach(label,0, 2,1,2)
        tableTime.attach(self.ftime,2,4,1,2)


        self.listtime = gtk.Combo()
        self.listtime.entry.set_text("time")
        listoption = [ "MET", "MJD", "JD" ]
        listoption.remove(self.config['time']['type'])
        listoption.insert(0,self.config['time']['type'])
        self.listtime.set_popdown_strings(listoption)
        self.listtime.show()

        label = gtk.Label("Unit (type)")
        label.show()
        tableTime.attach(label, 0,2, 2,3)
        tableTime.attach(self.listtime, 2,3, 2,3)

    def _addSpectrum(self):
        BNPage = self.AddBlocNotePage("Spectrum/Ebin")
        frame,table = self._addFrame("Spectrum options",2,4)
        BNPage.attach(frame, 0, 8, 0, 4)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "Spectrum")
        self.set_active(FitButton,self.config["Spectrum"]["FitsGeneration"])
        FitButton.show()

        PlotButton = gtk.CheckButton("")
        PlotButton.connect("toggled", self.fct_yesno, "Spectrum")
        self.set_active(PlotButton,self.config["Spectrum"]["ResultPlots"])
        PlotButton.show()

        self.index = self.AddSpinButton(self.config["Spectrum"]["FrozenSpectralIndex"], 0, 5, .1, 2)

        SummedButton = gtk.CheckButton("")
        SummedButton.connect("toggled", self.fct_yesno, "Spectrum")
        self.set_active(SummedButton,self.config["Spectrum"]["SummedLike"])
        SummedButton.show()

        labfits = gtk.Label("Generation of the fits files")
        labfits.show()
        table.attach(labfits,0, 2,0, 1)
        labplots = gtk.Label("Generation of the plots")
        labplots.show()
        table.attach(labplots,0, 2 ,1, 2)
        labindex = gtk.Label("Frozen Spectral Index value (no effect if 0)")
        labindex.show()
        table.attach(labindex,0, 3 ,2,3)
        labsumm = gtk.Label("Used the summed likelihood")
        labsumm.show()
        table.attach(labsumm,0, 2 ,3, 4)

        table.attach(FitButton, 3,4,0, 1)
        table.attach(PlotButton, 3,4,1, 2)
        table.attach(self.index, 3,4,2,3)
        table.attach(SummedButton, 3,4,3, 4)


        frameebin,tableebin = self._addFrame("Energy bins options",2,4)
        BNPage.attach(frameebin, 0, 8, 4, 7)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "Ebin")
        self.set_active(FitButton,self.config["Ebin"]["FitsGeneration"])
        FitButton.show()

        labfits = gtk.Label("Generation of the fits files")
        labfits.show()
        tableebin.attach(labfits,0, 2,0, 1)
        tableebin.attach(FitButton, 3,4,0, 1)

        self.nebin = self.AddSpinButton(self.config["Ebin"]["NumEnergyBins"], 0, 30, 1, 0)

        label = gtk.Label("Number of bins")
        label.show()
        tableebin.attach(label,0, 2,1,2)
        tableebin.attach(self.nebin,3,4,1,2)

        self.tsebin = self.AddSpinButton(self.config["Ebin"]["TSEnergyBins"],0, 1000, .1, 1)

        label = gtk.Label("Minimal TS")
        label.show()
        tableebin.attach(label,0, 2,2,3)
        tableebin.attach(self.tsebin,3,4,2,3)

        label = gtk.Label("Re-run the Ebin calculation only")
        label.show()
        ebinbutton = gtk.Button("Re-run Ebin")  
        ebinbutton.connect("clicked", self.LaunchEbin, '')
        ebinbutton.show()
        tableebin.attach(label,0, 2,3,4)
        tableebin.attach(ebinbutton,3,4,3,4)

    def _addUL(self):
        BNPage = self.AddBlocNotePage("Upper Limits")
        frame,table = self._addFrame("Upper Limits options",2,4)
        BNPage.attach(frame, 0, 8, 0, 4)

        self.ULindex = self.AddSpinButton(self.config["UpperLimit"]["SpectralIndex"], 0, 5, .1,2)

        MethodButton = gtk.CheckButton("")
        MethodButton.connect("toggled", self.fct_yesno, "ULenvelope")
        self.set_active(MethodButton,self.config["UpperLimit"]["envelope"])
        MethodButton.show()

        self.ULcl = self.AddSpinButton(self.config["UpperLimit"]["cl"], 0, 1, .01,2)
        self.TSlim = self.AddSpinButton(self.config["UpperLimit"]["TSlimit"], 0, 1000, .1, 1)

        labindex = gtk.Label("Assumed Spectral Index value")
        labindex.show()
        table.attach(labindex,0, 2 ,0,1)

        labts = gtk.Label("Minimal TS")
        labts.show()
        table.attach(labts,0, 2 ,1,2)

        labts = gtk.Label("Confidence level")
        labts.show()
        table.attach(labts,0, 2 ,2,3)

        labts = gtk.Label("Envelope UL")
        labts.show()
        table.attach(labts,0, 2 ,3,4)

        table.attach(self.ULindex, 3,4,0,1)
        table.attach(self.TSlim, 3,4,1,2)
        table.attach(self.ULcl, 3,4,2,3)
        table.attach(MethodButton, 3,4,3,4)

    def _addLC(self):
        BNPage = self.AddBlocNotePage("Light Curves")
        frame,table = self._addFrame("Light Curves options",2,4)
        BNPage.attach(frame, 0, 8, 0, 4)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "LightCurve")
        self.set_active(FitButton,self.config["LightCurve"]["FitsGeneration"])
        FitButton.show()

        labfits = gtk.Label("Generation of the fits files")
        labfits.show()
        table.attach(labfits,0, 2,0, 1)
        table.attach(FitButton, 3,4,0, 1)

        self.nlcbin = self.AddSpinButton(self.config["LightCurve"]["NLCbin"], 0, 1e6, 1, 0)
        label = gtk.Label("Number of bins")
        label.show()
        table.attach(label,0, 2,1,2)
        table.attach(self.nlcbin,3,4,1,2)

        self.lcindex = self.AddSpinButton(self.config["LightCurve"]["SpectralIndex"], 0, 5, .1, 2)
        label = gtk.Label("Spectral Index")
        label.show()
        table.attach(label,0, 2,2,3)
        table.attach(self.lcindex,3,4,2,3)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "conffile")
        self.set_active(FitButton,self.config["LightCurve"]["MakeConfFile"])
        FitButton.show()

        label = gtk.Label("Make config file")
        label.show()
        table.attach(label,0, 2,3,4)
        table.attach(FitButton,3,4,3,4)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "compvar")
        self.set_active(FitButton,self.config["LightCurve"]["ComputeVarIndex"])
        FitButton.show()

        label = gtk.Label("Compute variability index")
        label.show()
        table.attach(label,0, 2,4,5)
        table.attach(FitButton,3,4,4,5)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "lcdiagplot")
        self.set_active(FitButton,self.config["LightCurve"]["DiagnosticPlots"])
        FitButton.show()

        label = gtk.Label("Make diagnostic plots")
        label.show()
        table.attach(label,0, 2,5,6)
        table.attach(FitButton,3,4,5,6)



#    def _addFoldedLC(self):
#        BNPage = self.AddBlocNotePage("Folded LC")
#        frame,table = self._addFrame("Folded Light Curves options",2,4)
#        BNPage.attach(frame, 0, 8, 0, 4)

    def _addAppFoldedLC(self):
        BNPage = self.AddBlocNotePage("Apperture/Folded LC")
        frameapp,tableapp = self._addFrame("Apperture photometry options",2,4)
        BNPage.attach(frameapp, 0, 8, 0, 4)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "AppLC")
        self.set_active(FitButton,self.config["AppLC"]["FitsGeneration"])
        FitButton.show()

        labfits = gtk.Label("Generation of the fits files")
        labfits.show()
        tableapp.attach(labfits,0, 2,0, 1)
        tableapp.attach(FitButton, 3,4,0, 1)

        self.applcindex = self.AddSpinButton(self.config["AppLC"]["index"], 0, 5, .1, 2)
        label = gtk.Label("Spectral Index")
        label.show()
        tableapp.attach(label,0, 2,1,2)
        tableapp.attach(self.applcindex,3,4,1,2)

        self.applcNbin = self.AddSpinButton(self.config["AppLC"]["NLCbin"], 0, 1e6, 1, 0)
        label = gtk.Label("Number of bins")
        label.show()
        tableapp.attach(label,0, 2,2,3)
        tableapp.attach(self.applcNbin,3,4,2,3)


        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "binfromdata")
        self.set_active(FitButton,self.config["AppLC"]["binsFromData"])
        FitButton.show()

        labfits = gtk.Label("Make bin from data")
        labfits.show()
        tableapp.attach(labfits,0, 2,3,4)
        tableapp.attach(FitButton, 3,4,3,4)

        framefol,tablefol = self._addFrame("Folded LightCurves options",2,4)
        BNPage.attach(framefol, 0, 8, 4, 8)

        self.follcNbin = self.AddSpinButton(self.config["FoldedLC"]["NLCbin"], 0, 1e6, 1, 0)
        label = gtk.Label("Number of bins")
        label.show()
        tablefol.attach(label,0, 2,0,1)
        tablefol.attach(self.follcNbin,3,4,0,1)

        self.folepoch = self.AddSpinButton(self.config["FoldedLC"]["epoch"], 0, 1e10, 1, 0)
        label = gtk.Label("Epoch")
        label.show()
        tablefol.attach(label,0, 2,1,2)
        tablefol.attach(self.folepoch,3,4,1,2)

        self.folperiod = self.AddSpinButton(self.config["FoldedLC"]["Period"], 0, 1e10, 1, 0)
        label = gtk.Label("Period")
        label.show()
        tablefol.attach(label,0, 2,2,3)
        tablefol.attach(self.folperiod,3,4,2,3)

    def _addTSMap(self):
        BNPage = self.AddBlocNotePage("TS Map")
        frame,table = self._addFrame("TS Map options",2,4)
        BNPage.attach(frame, 0, 8, 0, 3)

        reFitButton = gtk.CheckButton("")
        reFitButton.connect("toggled", self.fct_yesno, "tsmap")
        self.set_active(reFitButton,self.config["TSMap"]["Re-Fit"])
        reFitButton.show()

        label = gtk.Label("Re-fit")
        label.show()
        table.attach(label,0, 2,0,1)
        table.attach(reFitButton, 3,4,0,1)

        self.tsmapnpix = self.AddSpinButton(self.config["TSMap"]["npix"], 0, 1000, 1, 0)
        label = gtk.Label("Number of pixel")
        label.show()
        table.attach(label,0, 2,1,2)
        table.attach(self.tsmapnpix,3,4,1,2)


        reFitButton = gtk.CheckButton("")
        reFitButton.connect("toggled", self.fct_yesno, "removetgr")
        self.set_active(reFitButton,self.config["TSMap"]["RemoveTarget"])
        reFitButton.show()

        label = gtk.Label("Remove Target")
        label.show()
        table.attach(label,0, 2,2,3)
        table.attach(reFitButton, 3,4,2,3)

        self.listtsmethod = gtk.Combo()
        self.listtsmethod.entry.set_text("method")
        listoption = [ "row", "pixel" ]
        listoption.remove(self.config['TSMap']['method'])
        listoption.insert(0,self.config['TSMap']['method'])
        self.listtsmethod.set_popdown_strings(listoption)
        self.listtsmethod.show()

        label = gtk.Label("Method to use")
        label.show()
        table.attach(label, 0, 2,3,4)
        table.attach(self.listtsmethod, 3,4,3,4)

    def _addfindsrc(self):
        BNPage = self.AddBlocNotePage("Findsrc")
        frame,table = self._addFrame("Find source options",2,4)
        BNPage.attach(frame, 0, 8, 0, 2)

        FitButton = gtk.CheckButton("")
        FitButton.connect("toggled", self.fct_yesno, "findsrc")
        self.set_active(FitButton,self.config["findsrc"]["FitsGeneration"])
        FitButton.show()

        labfits = gtk.Label("Generation of the fits files")
        labfits.show()
        table.attach(labfits,0, 2,0, 1)
        table.attach(FitButton, 3,4,0, 1)

        reFitButton = gtk.CheckButton("")
        reFitButton.connect("toggled", self.fct_yesno, "findsrc")
        self.set_active(FitButton,self.config["findsrc"]["FitsGeneration"])
        reFitButton.show()

        label = gtk.Label("Re-fit")
        label.show()
        table.attach(label,0, 2,1,2)
        table.attach(reFitButton, 3,4,1,2)

    def _addPlot(self):

        BNPage = self.AddBlocNotePage("Plots")
        frame,table = self._addFrame("Results and debug plots",2,4)
        BNPage.attach(frame, 0, 8, 0, 8)


        filebase= self.config['out'] + '/Spectrum/SED_' + self.config['target']['name'] +'_'+ self.config['target']['spectrum']
        self.image = gtk.Image()
        pixbuf = gtk.gdk.pixbuf_new_from_file(os.environ.get('ENRICO_DIR', '')+'/enrico/enrico.jpg')
        pixbuf = pixbuf.scale_simple(250, 300, gtk.gdk.INTERP_BILINEAR)
        self.image.set_from_pixbuf(pixbuf)
        self.image.show()
        table.attach(self.image,0,8,1,7)

        sedbutton = gtk.Button("SED")
        sedbutton.connect("clicked", self.change, filebase+".png")
        table.attach(sedbutton, 0,1,0,1)
        sedbutton.show()

        countbutton = gtk.Button("count")
        countbutton.connect("clicked", self.change, filebase+"_CountsPlot.png")
        table.attach(countbutton, 1,2,0,1)
        countbutton.show()

        resbutton = gtk.Button("residuals")
        resbutton.connect("clicked", self.change, filebase+"_ResPlot.png")
        table.attach(resbutton, 2,3,0,1)
        resbutton.show()

        sep = gtk.VSeparator()
        sep.show()
        table.attach(sep, 3,4,0,1)

        lcfiles = self.config["out"]+"/LightCurve_"+str(self.config["LightCurve"]["NLCbin"])+"bins/"
        lcbutton = gtk.Button("LC")
        lcbutton.connect("clicked", self.change, lcfiles+"_LC.png")
        table.attach(lcbutton, 4,5,0,1)
        lcbutton.show()

        tsbutton = gtk.Button("TS vs time")
        tsbutton.connect("clicked", self.change, lcfiles+"_TS.png")
        table.attach(tsbutton, 5,7,0,1)
        tsbutton.show()

        npredbutton = gtk.Button("Npred")
        npredbutton.connect("clicked", self.change, lcfiles+"_Npred.png")
        table.attach(npredbutton, 7,8,0,1)
        npredbutton.show()


    def __init__(self,infile):
        self.infile = infile
        try :
            self.config = get_config(infile)
        except :
            config = ConfigObj(indent_type='\t')
            config['out'] = os.getcwd()
            self.config = get_config(config) 

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.connect("delete_event", self.delete)
        window.set_border_width(10)

        table = gtk.Table(3,6,False)
        window.add(table)

        self.notebloc = gtk.Notebook()
        self.notebloc.set_tab_pos(gtk.POS_LEFT)
        table.attach(self.notebloc, 0,6,0,1)
        self.notebloc.show()
        self._addMainOption()
        self._addFiles()
        self._addTargetSpace()
        self._addAnalysis()
        self._addEnergyTime()
        self._addSpectrum()
        self._addUL()
        self._addLC()
#        self._addFoldedLC()
        self._addAppFoldedLC()
#        self._addEbin()
        self._addTSMap()
        self._addfindsrc()
        self._addPlot()

        self.notebloc.set_current_page(0)

        SaveButton = gtk.Button("Save file")
        SaveButton.connect("clicked", self.save, "")
        table.attach(SaveButton, 1,2,1,2)
        SaveButton.show()

        CloseButton = gtk.Button("Close")
        CloseButton.connect("clicked", self.delete)
        table.attach(CloseButton, 4,5,1,2)
        CloseButton.show()

        table.show()
        window.show()


if __name__ == "__main__":

    try:
        infile = sys.argv[1]
    except:
        logging.error('Config file not found.')
        print('Usage: '+sys.argv[0]+' <output config file name>')
        sys.exit(1)

    g=EnricoGui(infile)
    gtk.main()

