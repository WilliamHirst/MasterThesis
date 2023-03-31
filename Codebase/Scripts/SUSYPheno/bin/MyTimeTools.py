####################################################################
### 2010 July28  Borge Kile Gjelsten (b.k.gjelsten@fys.uio.no)
####################################################################

#from time import time
import time
#from math import log,fmod
import datetime


# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
def DatetimeFrom1970secs(sec):
  return datetime.datetime.utcfromtimestamp(sec)

def StampFrom1970secs(sec, mode=-1, what=[]):
  dt = DatetimeFrom1970secs(sec)
  stamp = StampFromDatetime(dt, mode=mode, what=what)
  return stamp


def TimeSince(stampT0, unit='d', stampT1=''):
  T0 = DatetimeFromStamp(stampT0)
  if stampT1: T1 = DatetimeFromStamp(stampT1)
  else: T1 = datetime.datetime.now()
  tdel = T1 - T0
  tdel_s = tdel.days * 24*3600 + tdel.seconds + tdel.microseconds
  if unit == 's': return int(tdel_s)
  if unit == 'm': return int(tdel_s/(60))
  if unit == 'h': return int(tdel_s/(60*60))
  if unit == 'd': return int(tdel_s/(60*60*24))
  if unit == 'w': return int(tdel_s/(60*60*24*7))
  if unit == 's_dec': return float(tdel_s)
  if unit == 'm_dec': return float(tdel_s/(60.))
  if unit == 'h_dec': return float(tdel_s/(60.*60))
  if unit == 'd_dec': return float(tdel_s/(60.*60*24))
  if unit == 'w_dec': return float(tdel_s/(60.*60*24*7))
  print('WARNING: Unknown unit: %s   Returning in seconds') 
  return int(tdel_s)



def DatetimeFromStamp(stamp, alter=[], mode=0):
  #print 'DEBUG: ', stamp
  # Ex:  DatetimeFromStamp('2011-10-06_12:13', alter={'w':2, 's':-1})

  # 1. Dechiffer the stamp
  # assume here stamp has format of type 2011-12-12_12:14:13 (mode0) (and h, m and s as well)
  W = stamp.replace('-',' ').replace('_',' ').replace(':',' ').replace('h',' ').replace('m',' ').replace('s','').split()
  
  z = []
  nw = len(W)
  for iw in range(6):
    if iw<nw: z.append(int(W[iw]))
    else: z.append(0)

  # 2. allows to add to / subtract from the stamp using datetime.timedelta
  w=0; d=0; h=0; m=0; s=0;
  if 'w' in alter: w = alter['w']
  if 'd' in alter: d = alter['d']
  if 'h' in alter: h = alter['h']
  if 'm' in alter: m = alter['m']
  if 's' in alter: s = alter['s']

  # 3. the result
  dt0 = datetime.datetime(z[0],z[1],z[2],z[3],z[4],z[5])
  diff = datetime.timedelta(weeks=w, days=d, hours=h, minutes=m, seconds=s)
  dt = dt0 + diff
  return dt


def fromTxtStamp(stamp):  # 2013-10-30
  res = {}
  res['datetime'] = DatetimeFromStamp(stamp)
  res['secs1970'] = time.mktime(res['datetime'].timetuple())

  return res


def StampFromDatetime(da, mode=-1, what=[]):
  if mode==-1 and not what: mode = 0
  if mode == 0: what = ['Y','M','D','h','m','s']
  if mode == 1: what = ['Y','M','D','h','m']
  if mode == 2: what = ['Y','M','D']

  res = ''
  if 'Y' in what: res += "%4i"  %(da.year)
  if 'M' in what: res += "-%2s" %(str(da.month).zfill(2))
  if 'D' in what: res += "-%2s" %(str(da.day).zfill(2))
  if 'h' in what: res += "_%2s" %(str(da.hour).zfill(2))
  if 'm' in what: res += ":%2s" %(str(da.minute).zfill(2))
  if 's' in what: res += ":%2s" %(str(da.second).zfill(2))

  if res.startswith(':'): res = res[1:]
  if res.startswith('-'): res = res[1:]
  if res.startswith('_'): res = res[1:]

  return res
  

def GetTimeStamp(mode=0):
  if mode == 0: #ex: Jan15_155950'
    t = time.ctime().split()
    stamp = t[1]+t[2].zfill(2)+"_"+t[3].replace(":","")
    return stamp
  
  if mode == 1: #ex: 2011-01-14_15h45 (used by SusySpaceTools)
    t = list(time.localtime())
    for iT in range(len(t)): t[iT] = str(t[iT])
    stamp = t[0].zfill(4)+'-'+t[1].zfill(2)+'-'+t[2].zfill(2)+'_'+t[3].zfill(2)+'h'+t[4].zfill(2)
    return stamp

  if mode == 2: #ex: 2011-01-14_15-45 (newer 'standard' used many places)
    t = list(time.localtime())
    for iT in range(len(t)): t[iT] = str(t[iT])
    stamp = t[0].zfill(4)+'-'+t[1].zfill(2)+'-'+t[2].zfill(2)+'_'+t[3].zfill(2)+'-'+t[4].zfill(2)
    return stamp

  if mode == 3: #ex: 2011-01-14_15-45-32
    t = list(time.localtime())
    for iT in range(len(t)): t[iT] = str(t[iT])
    stamp = t[0].zfill(4)+'-'+t[1].zfill(2)+'-'+t[2].zfill(2)+'_'+t[3].zfill(2)+'-'+t[4].zfill(2)+'-'+t[5].zfill(2)
    return stamp

  if mode == 4: #ex: 2011-01-14_15
    t = list(time.localtime())
    for iT in range(len(t)): t[iT] = str(t[iT])
    stamp = t[0].zfill(4)+'-'+t[1].zfill(2)+'-'+t[2].zfill(2)+'_'+t[3].zfill(2)
    return stamp

  if mode == 5: #ex: 2011-01-14
    t = list(time.localtime())
    for iT in range(len(t)): t[iT] = str(t[iT])
    stamp = t[0].zfill(4)+'-'+t[1].zfill(2)+'-'+t[2].zfill(2)
    return stamp


  t = list(time.localtime())
  for iT in range(len(t)): t[iT] = str(t[iT])
  stamp = t[0].zfill(4)+t[1].zfill(2)+t[2].zfill(2)+'_'+t[3].zfill(2)+t[4].zfill(2)+t[5].zfill(2)
  print('illegal mode: %i. Using stamp = %s' %(mode,stamp))
  return stamp


#class MyTime():
class MyTime:
  def __init__(self,txt):
    self.txt = txt
    self.ctime = time.ctime()
    self.localtime = time.localtime()
    self.time = time.time()
    

    
class MyTimer:
  def __init__(self):
    self.timer = []


  def Add(self,txt="",showwhich=-9, mode=0, modeprev=3, get=0):
    if(txt==""): txt="time_"+str(len(self.timer))
    self.timer.append(MyTime(txt))
    #self.timer.append([time(), txt])
    if showwhich>-9:
      return self.Show(which=showwhich, mode=mode, modeprev=modeprev, get=get)
        

  def ShowEstimate(self,fDone,txt="", mode=2):  # Method added 2010,July28
    # fDone is the fraction done, typically (float)thisEvent/totEvent
    which = len(self.timer)-1
    tStart = self.timer[0].time
    tPrev  = self.timer[len(self.timer)-2].time
    tNow   = self.timer[len(self.timer)-1].time
    
    dtLast  = tNow-tPrev
    dtSofar = tNow-tStart
    dtEnd   = dtSofar/fDone

    if mode == 0:
      return dtLast, dtSofar, dtEnd

    if mode == 1: 
      return FormatTime(dtLast), FormatTime(dtSofar), FormatTime(dtEnd)

    if mode == 2:
      return "delT(Last): %s   delT(SoFar): %s    delT(estEnd): %s" %(FormatTime(dtLast), FormatTime(dtSofar), FormatTime(dtEnd))
    
    
    


  def Show(self, which=-1, mode=0, modeprev=3, get=0):
    # which: -2: all
    # which: -1: latest
    # which: 0,...: the one asked for
    # mode: 0  [1 for earthtime too, seldom used]

    outs = []
    for iT in range(len(self.timer)):
      out_time = self.timer[iT].ctime
      #out = "TIME%4s: %-20s" %('['+str(iT)+']',self.timer[iT].txt)  #preDec12
      out = ""

      if mode==0: out += "TIME%4s: " %('['+str(iT)+']')  #hack

      if mode==1: out += "  %s" %out_time #show earth time as well

      # out += "   tot = %s" %FormatTime(self.Diff(iT,0))
      out += "   tot = %s" %FormatTime(self.Diff(iT,0), mode=2, space=0, ON=['m','s'])
      
      if iT >= 0:  ### >  # 2013-10-17 changed from >0 to >=0 Hmm
        # stot = self.Diff(iT)
        # print "stot = "+str(stot)
        out_diff = FormatTime(self.Diff(iT), mode=2, space=0, ON=['m','s'])
        #out +="      diff = %s   [prev = %-s]" %(out_diff, self.timer[iT-1].txt) #preDec12

        if modeprev & 1: out += "      [diff = %s]" %(out_diff)           # show diff
        out += "  %-20s" %(self.timer[iT].txt)        
        if modeprev & 2: out += "   [prev = %-s]" %(self.timer[iT-1].txt) # show prev

    
      if which==-2 or (which==-1 and iT==len(self.timer)-1) or which==iT:
        if get: outs.append(out)
        else: print(out)
        #if get: outs.append(out)
        #else: print out

    if get:
      if len(outs) == 1: return outs[0]  # legacy
      else: return outs  # new (e.g. for log->file)



  def Diff(self,iT1,iT0=-9):
    if iT1>=len(self.timer) or iT1<=0:
      diff = 0
    else:
      if iT0==-9: iT0=iT1-1
      diff = self.timer[iT1].time - self.timer[iT0].time
      # print "from Diff will return: %i" %diff
      
    return diff
      
# <><><><><> <><><><><> <><><><><> <><><><><> <><><><><> <><><><><>
  
def FormatTimeBack(ftime, unit='s'):
  # Ex.: FormatTimeBack(ftime='13h04m17s', unit='m')
  units = ['s','m','h','d'] #valid return units
  # First find all in seconds, then turn into preferred unit
  s = 0
  if 'd' in ftime:
    z = ftime.split('d')
    s += 60*60*24 * int(z[0])
    ftime = z[1]
  if 'h' in ftime:
    z = ftime.split('h')
    s += 60*60 * int(z[0])
    ftime = z[1]
  if 'm' in ftime:
    z = ftime.split('m')
    s += 60 * int(z[0])
    ftime = z[1]
  if 's' in ftime:
    z = ftime.split('s')
    s += 1 * int(z[0])
    ftime = z[1]

  if unit == 's': return s
  elif unit == 'm': return s / (60.)
  elif unit == 'h': return s / (60.*60)
  elif unit == 'd': return s / (60.*60*24)
  else: return -1
  

def FormatTime(stot,mode=0,delim=" ",ON = ['h','m','s'],space=1):  
  # print "FormatTime: stot = "+str(stot)
  h = int(stot/3600)
  m = int((stot-h*3600)/60)
  s = int(stot-h*3600-m*60)
  lineout = ""
  if mode==0:
    if h>0 and 'h' in ON: lineout += str(h)+"h"+delim
    if (h>0 or m>0) and 'm' in ON: lineout += str(m).zfill(2)+"m"+delim
    if 's' in ON: lineout += str(s).zfill(2)+"s"
  elif mode==1: 
    # lineout = txt+str(h)+"h "+str(m).zfill(2)+"m "+str(s).zfill(2)+"s" #missing txt...
    if 'h' in ON: lineout += str(h)+"h"
    if 'm' in ON: lineout += delim + str(m).zfill(2)+"m"
    if 's' in ON: lineout += delim + str(s).zfill(2)+"s"
    lineout = lineout.strip()
  elif mode==2: 
    # lineout = txt+str(h)+"h "+str(m).zfill(2)+"m "+str(s).zfill(2)+"s" #missing txt...
    if 'h' in ON or stot>=3600: lineout += str(h)+"h"
    if 'm' in ON or stot>=60: lineout += delim + str(m).zfill(2)+"m"
    if 's' in ON or stot>=0: lineout += delim + str(s).zfill(2)+"s"
    lineout = lineout.strip()

  if not space: lineout = lineout.replace(' ','')
  return lineout
  
  # <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>


class MyClock:
  
  def __init__(self):
    self.now = datetime.datetime.now()
    #self.dayno0 = ["noday", "Mon","Tue","Wed","Thu","Fri", "Sat","Sun"]
    #self.day0 = ["Mon","Tue","Wed","Thu","Fri", "Sat","Sun"]
    self.dayI0 = {"Mon":0,"Tue":1,"Wed":2,"Thu":3,"Fri":4, "Sat":5,"Sun":6}
    self.dayI0.update({"mon":0,"tue":1,"wed":2,"thu":3,"fri":4, "sat":5,"sun":6})
    self.dayI1 = {"Mon":1,"Tue":2,"Wed":3,"Thu":4,"Fri":5, "Sat":6,"Sun":7}
    self.dayI1.update({"mon":1,"tue":2,"wed":3,"thu":4,"fri":5, "sat":6,"sun":7})
    self.day1 = self.now.isoweekday()
    self.day0 = self.now.isoweekday()-1

  #def day0(self): return self.localtime.tm_day
  #def day1(self): return self.now.isoweekday()
  #def day0(self): return self.now.isoweekday()-1

  def dateofDayInWeek(self,dayT,week=0, mode=0):
    day_now = self.day1
    day_seek = self.dayI1[dayT]
    ndays = day_seek-day_now
    dt_seek = self.now + datetime.timedelta(days=ndays, weeks=week)

    if mode == 0:
      txt = "%4i-%2s-%2s" %(dt_seek.year, str(dt_seek.month).zfill(2), str(dt_seek.day).zfill(2))
      return txt
    
    
####################################################
    
def GetTimeOverview(TIME, TOT='tot', keys=[], mode=[]):
  '''
  TIME is a dict with times (typically the time per subroutine/task)
  Here return string showing total time (in s) and the various contributions (relative)
  '''
  
  if 'tot' not in TIME: sys.exit('FATAL::GetTimeString: no key %s in TIME' %(TOT))
  if not keys: keys = list(TIME.keys())

  txt = ''
  if 'notot' not in mode: txt = '%s: %.1f' %(TOT, TIME[TOT])
  for key in keys:
    if key == TOT: continue
    txt += '  %s:%.2f' %(key, TIME[key]/TIME[TOT])

  return txt


####################################################
daystxt3 = ['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
daystxt3_lower = ['mon','tue','wed','thu','fri','sat','sun']

#####
def timearrayFromStamp(stamp):
  T = stamp.replace('-',' ').replace('_',' ').replace(':',' ').replace('h',' ').replace('m',' ').replace('s','').split()
  timearray = [0,0,0,0,0,0]
  for iT in range(len(T)): timearray[iT] = int(T[iT])
  return timearray


#####
def CleanDT(dt, optD={}, optL=[]):
      if 'year'        in optL: dt = dt.replace(year=0)
      if 'month'       in optL: dt = dt.replace(month=0)
      if 'day'         in optL: dt = dt.replace(day=0)
      if 'hour'        in optL: dt = dt.replace(hour=0)
      if 'minute'      in optL: dt = dt.replace(minute=0)
      if 'second'      in optL: dt = dt.replace(second=0)
      if 'microsecond' in optL: dt = dt.replace(microsecond=0)
      return dt


#####
def GetDTlastXXXday(dtref, dayI=-1, dayT='', clean=['hour','minute','second','microsecond']):
  # E.g. GetLastXXXday(dayT='mon') ; finds DT of last monday (including today if monday)
  if dayI == -1 and dayT == '':
    print('Warning: GetDTlast')
  if dayI == -1: dayI = daystxt3_lower.index(dayT.lower())

  deltadays = dtref.day - dayI
  if deltadays < 0: deltadays += 7
  
  dtlast = dtref - datetime.timedelta(days=deltadays)
  
  # Now clean off minutes etc. if
  dtlast = CleanDT(dtlast, optL=clean)
  
  return dtlast


#####
def GetDTmonthincrement(dtref, monthincrement): 
  year = dtref.year
  month = dtref.month

  step = 1
  if monthincrement < 0: step = -1

  # increment (or decrement)
  for i in range(abs(monthincrement)):
    month += step
    if month == 13:
      month = 1
      year += 1
    if month == 0:
      month = 12
      year -= 1
    
  dt = dtref.replace(year=year).replace(month=month)

  return dt


####################################################
class MyTimeCalc:  # not (yet) in use

  #####
  def ___init___(self, optD={}, optL=[]):

    s = self

    # constants

    # -------

    # Set time "now"
    if 'time' in optD: s.dt = datetime.datetime(optD['time'])
    elif 'timetxt' in optD: s.dt = datetime.datetime(timearrayFromStamp(optD['timetxt']))
    else: s.dt = datetime.datetime.now()


  #####
  

####################################################
