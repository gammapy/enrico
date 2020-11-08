#!/usr/bin/env python

import time,sys



epoch_ref = 1331679600+3600
mjd_ref = 56000.
met_ref = 353376002.000

def MJDtoJD(MJD):
  return MJD+2400000.5

def JDtoMJD(JD):
  return JD-2400000.5

def EpochToMeT(epoch):
# Work only on Linux Machine
  depoch = epoch-epoch_ref
  return depoch + met_ref+3600

def DateToMet(date):
# Work only on Linux Machine
  pattern = '%Y-%m-%d %H:%M:%S'
  try :
    epoch = int(time.mktime(time.strptime(date, pattern)))
  except:
    pattern = '%Y-%m-%d'
    epoch = int(time.mktime(time.strptime(date, pattern)))
  return EpochToMeT(epoch)


def MetToDate(met):
# Work only on Linux Machine
  epoch = met-met_ref+epoch_ref
  pattern = '%Y-%m-%d %H:%M:%S'
  return time.strftime(pattern,time.gmtime(epoch))


def MJDToMeT(MJD):
  return (MJD-mjd_ref)*3600*24+met_ref

def MeTToMJD(met):
  return (met-met_ref)/3600/24.+mjd_ref

def Print(date,mjd,met,jd):
  print(("ISO date : ",date))
  print(("MJD : ",mjd))
  print(("JD : ",jd))
  print(("MET : ",met))

def _log_():
    print(("Usage: ",sys.argv[0]," Type Date"))
    print("Type can be : MJD, JD, MET or ISO")
    print("Data is the date to convert")
    print(("exemple:\n\tpython ",sys.argv[0]," MJD 56101"))
    print(("\tpython ",sys.argv[0]," ISO 2012-02-05"))
    print(("\tpython ",sys.argv[0]," ISO 2012-02-05 16:15:18"))
if __name__=="__main__":

  try :
    Date_type = sys.argv[1]
  except :
    _log_()
    exit()


  if Date_type=="ISO":
    try :
      iso = sys.argv[2]+" "+sys.argv[3]
    except :
      iso = sys.argv[2]
    Met = DateToMet(iso)
    MJD = MeTToMJD(Met)
    JD = MJDtoJD(MJD)

  elif Date_type=="MJD":
    MJD = float(sys.argv[2])
    JD = MJDtoJD(MJD)
    Met = MJDToMeT(MJD)
    iso = MetToDate(Met)

  elif Date_type=="MET":
    Met = float(sys.argv[2])
    MJD = MeTToMJD(Met)
    JD = MJDtoJD(MJD)
    iso = MetToDate(Met)

  elif Date_type=="JD":
    JD = float(sys.argv[2])
    MJD = JDtoMJD(JD)
    Met = MJDToMeT(MJD)
    iso = MetToDate(Met)


  else:
    _log_()
    exit()

Print(iso,MJD,Met,JD)
