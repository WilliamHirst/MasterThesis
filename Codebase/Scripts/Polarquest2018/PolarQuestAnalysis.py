import ROOT
from array import array
import os,sys
from os import listdir
from os.path import isfile, join
from datetime import datetime
import math
from ROOT import gROOT
gROOT.ProcessLine("gErrorIgnoreLevel = kError")

ROOT.TH1.SetDefaultSumw2()

# Function to retrieve difference in seconds between two time stamps
def getDuration(start,stop):
    timestamp_start = "%d %06d" %(t.GetDate(True,start),t.GetTime(True,start))
    timestamp_stop  = "%d %06d" %(t.GetDate(True,stop),t.GetTime(True,stop))
    
    datetime_object_start = datetime.strptime(timestamp_start, '%Y%m%d %H%M%S')
    datetime_object_stop  = datetime.strptime(timestamp_stop,  '%Y%m%d %H%M%S')
    difference = datetime_object_stop - datetime_object_start
    seconds_in_day = 24 * 60 * 60
    if difference.days < 0:
        return difference.seconds
    else:
        return difference.days * seconds_in_day + difference.seconds
    

# Folder where the data is stored
indir = sys.argv[1]

# The three trees containing different data
c_trend = ROOT.TChain("Trending")
c_headr = ROOT.TChain("Header")
c_weath = ROOT.TChain("Weather")

# Define time span (default: defined as the expedition with Nanuq (PolarquEEEst)
start_time = ROOT.TTimeStamp(2018,7,21,0,0,0)
stop_time  = ROOT.TTimeStamp(2018,9,5,0,0,0)

# Loop over files in input directory and select those
# which are inside the time span defined above
onlyfiles = [f for f in listdir(indir) if isfile(join(indir, f))]
files_to_use = []
for fname in sorted(onlyfiles):
    if not fname.endswith("summary.root"): continue

    # Get start date and stop date from file name
    st_y,st_m,st_d = (fname.split("_")[1]).split("-")
    start_dt = ROOT.TTimeStamp(int(st_y),int(st_m),int(st_d),0,0,0).AsDouble()
    en_y,en_m,en_d = (fname.split("_")[2]).split("-")
    stop_dt  = ROOT.TTimeStamp(int(en_y),int(en_m),int(en_d),0,0,0).AsDouble()

    # check if within the specified time period 
    if start_dt > start_time.AsDouble() and stop_dt < stop_time.AsDouble(): 
        files_to_use.append(join(indir,fname))
        #print fname

# Print summary and add files to TChains
print("INFO \t Found %i files in wanted time span from %s to %s" %(len(files_to_use),start_time.AsString("s"),stop_time.AsString("s")))
for ftu in files_to_use:
    c_trend.Add(ftu)
    c_headr.Add(ftu)
    c_weath.Add(ftu)

# The time offset used in the POLA data is 1st of January 2007
t = ROOT.TTimeStamp(2007,1,1,0,0,0)
da = ROOT.TDatime(2007,1,1,0,0,0)
ROOT.gStyle.SetTimeOffset(da.Convert());

# Arrays to store values
x_rawrate = array( 'd' )
y_rawrate = array( 'd' )

x_event = array('d')
y_pressure = array('d')
y_temp_in = array('d')

# Define some histograms
h_nev_temp = ROOT.TH1F("h_nev_temp","",60,0,60)
h_time_temp = ROOT.TH1F("h_time_temp","",60,0,60)
h_err_temp = ROOT.TH1F("h_err_temp","",60,0,60)

h_nev_pressure = ROOT.TH1F("h_nev_pressure","",1500,0,1500)
h_time_pressure = ROOT.TH1F("h_time_pressure","",1500,0,1500)
h_err_pressure = ROOT.TH1F("h_err_pressure","",1500,0,1500)

# Get the number of entries
nentries_headr = c_headr.GetEntries()
nentries_weath = c_weath.GetEntries()
nentries_trend = c_trend.GetEntries()

print("Number of runs   :", nentries_headr)
print("Number of events :", nentries_trend)

if nentries_headr != nentries_weath:
    sys.exit("ERROR \t Number of events in header and weather do not match!. Exiting")


# Set vb to non-zero value if you want some extra print-out
vb = 0
storeStartTime_run = True
tot_numevents = 0.0
tot_duration = 0
# Loop over runs
for i in range(nentries_headr):
    if i%100 == 0: print("Run %i/%i" %(i,nentries_headr))

    # Get entry i from TTrees (important!)
    c_headr.GetEntry(i)
    c_weath.GetEntry(i)

    # Exclude events with few events
    if c_headr.NumEvents < 12000: continue
    
    # Get duration, start and end time for run (in seconds from 1/1/2007) 
    start = int(c_headr.RunStart)
    stop  = int(c_headr.RunStop)

    run_duration  = getDuration(start,stop)

    # Some runs seem to have some bad start and stop values,
    # exclude those
    if run_duration > 1000000:
        print("skipping", fname)
        continue
    # Skip short runs (less than 10 minutes)
    if (run_duration/60.) < 10: continue

    # store start of time bin
    if storeStartTime_run: 
        start_new_block = start
        storeStartTime_run = False
        
    if vb > 0:
        print("Run duration = %d" %(run_duration))
        print("Start: %d %d" %(t.GetDate(True,start),t.GetTime(True,start)))
        print("Stop : %d %d" %(t.GetDate(True,stop),t.GetTime(True,stop)))
        

    # Accumulate number of events and time
    tot_duration  += run_duration
    tot_numevents += c_headr.NumEvents

    # If reached the wanted time span for adding a new point (default: 12 hours)
    if tot_duration >= 12.*60.*60.:
        # Fill x-values with the time (in seconds) in the middle of the range
        x_rawrate.append((start_new_block+(tot_duration/2.)))
        # Fill the raw rate in y-value 
        y_rawrate.append(tot_numevents/tot_duration)
        # reset variables
        tot_duration = 0.0
        tot_numevents = 0.0
        
        storeStartTime_run = True

    # Fill histograms for every event
    h_nev_temp.Fill(c_weath.IndoorTemperature,c_headr.NumEvents)
    h_time_temp.Fill(c_weath.IndoorTemperature,run_duration)
    
    # This histogram is needed to get the proper statistical uncertainty
    # (i.e. same histograms as above, but storing the number of entries, N)
    # Error is given by sqrt(N)
    h_err_temp.Fill(c_weath.IndoorTemperature)

# Loop over all the events in this run
tot_duration_ev = 0
tot_intemp_ev = 0.0
tot_pressure_ev = 0.0
nev_in_block = 0
storeStartTime_ev = True
for j in range(nentries_trend):
    c_trend.GetEntry(j)

    if j%1000 == 0: print("Event %i/%i" %(j,nentries_trend))

    # Start and stop time of events
    ev_start = int(c_trend.BinStart)
    ev_stop  = int(c_trend.BinEnd)

    # If start of a new time interval
    if storeStartTime_ev:
        ev_start_block = ev_start
        storeStartTime_ev = False

    # Get length of run (in seconds)
    ev_duration = getDuration(ev_start,ev_stop)

    # Some runs seem to have invalid start/stop times
    # Exclude those which are outside the time range defined at the beginning of the program
    # (i.e. the time of the PolarquEEEst experiment)
    if ev_stop > stop_time.AsDouble():
        continue
    # Skip short runs
    if ev_duration > 61: continue

    tot_duration_ev  += ev_duration
    tot_intemp_ev    += c_trend.IndoorTemperature
    tot_pressure_ev  += c_trend.Pressure
    nev_in_block += 1

    # If total duration is larger then the chosen time interval (default: 10 min)
    if tot_duration_ev >= 10.*60.:
        if vb > 0:
            print("Tot duration for data period ", (tot_duration_ev/(60.)))

        # Fill x-axis (time in middle of interval)
        x_event.append((ev_start_block+(tot_duration_ev/2.)))

        # Fill y-values (temperature, pressure)
        y_pressure.append(tot_pressure_ev/nev_in_block)
        y_temp_in.append(tot_intemp_ev/nev_in_block)

        # Reset counters
        storeStartTime_ev = True
        tot_duration_ev = 0.0
        tot_intemp_ev = 0.0
        tot_pressure_ev = 0.0
        nev_in_block = 0    


C = ROOT.TCanvas("c", "c", 1200, 600)
g_RawRate = ROOT.TGraph(len(x_rawrate),x_rawrate,y_rawrate)
# The default marker in ROOT is horrible so let's change
g_RawRate.SetMarkerStyle(20)
# Get x-axis in date and time instead of seconds since 2007
g_RawRate.GetXaxis().SetTimeDisplay(1)
g_RawRate.GetXaxis().SetTitle("Date")
g_RawRate.GetYaxis().SetTitle("Raw Rate [Hz]")
# Format the axis (default is date and time split in two lines)
g_RawRate.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
g_RawRate.GetXaxis().SetLabelOffset(0.03)
# Draw with axis and points
g_RawRate.Draw("AP")
C.Draw()
C.SaveAs("./rawrate.png")


C1 = ROOT.TCanvas("c1", "c1", 1200, 600)
g_pressure = ROOT.TGraph(len(x_event),x_event,y_pressure)
# The default marker in ROOT is horrible so let's change
g_pressure.SetMarkerStyle(20)
# Get x-axis in date and time instead of seconds since 2007
g_pressure.GetXaxis().SetTimeDisplay(1)
g_pressure.GetXaxis().SetTitle("Date")
g_pressure.GetYaxis().SetTitle("Pressure [hPa]")
# Format the axis (default is date and time split in two lines)
g_pressure.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
g_pressure.GetXaxis().SetLabelOffset(0.03)
# Draw with axis and points
g_pressure.Draw("AP")
C1.Draw()
C1.SaveAs("./pressure.png")


C3 = ROOT.TCanvas("c3", "c3", 1200, 600)
h_div_temp = h_nev_temp.Clone("h_div_temp")
h_div_temp.SetMarkerStyle(20)
# Set the errors correctly
for i in range(1,h_div_temp.GetNbinsX()+1):
    h_nev_temp.SetBinError(i,h_err_temp.GetBinError(i))
    h_time_temp.SetBinError(i,h_err_temp.GetBinError(i))
h_div_temp.Divide(h_nev_temp,h_time_temp)
h_div_temp.Draw("e")
h_div_temp.GetXaxis().SetTitle("Temperature [^{o} C]")
h_div_temp.GetYaxis().SetTitle("Raw Rate [Hz]")
C3.Draw()
C3.SaveAs("./temperature.png")

