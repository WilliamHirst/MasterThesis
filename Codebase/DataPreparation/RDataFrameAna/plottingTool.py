import ROOT as R
import sys
from samples import configure_samples

d_samp,d_type,d_reg = configure_samples()#False,False,True,False,False)


astyle = "/home/wohirst/atlasrootstyle"
R.gROOT.SetMacroPath(astyle)
R.gROOT.LoadMacro("AtlasUtils.C")
R.gROOT.LoadMacro("AtlasStyle.C")
R.gROOT.LoadMacro("AtlasLabels.C")
R.SetAtlasStyle()

def myLine(x1, y1, x2, y2, color=R.kBlack, style = 1, width = 2) :
  
  l = R.TLine()
  l.SetNDC()
  l.SetLineColor(color)
  l.SetLineStyle(style)
  l.DrawLine(x1,y1,x2,y2)

def myText(x, y, text, tsize=0.05, color=R.kBlack, angle=0) :
  
  l = R.TLatex()
  l.SetTextSize(tsize)
  l.SetNDC()
  l.SetTextColor(color)
  l.SetTextAngle(angle)
  l.DrawLatex(x,y,'#bf{' + text + '}')
  l.SetTextFont(4)

def findKey(histname):
  sp = histname.split("_")
  for i in range(len(sp)):
    s = "_".join(sp[i:])
    if s in d_samp.keys():
      return s
  return ""

class Plot:

    def __init__(p,hdic,hname = "lepPt_ele_hT", bkgs = [], is1D = True, xtext = "Feature 1", doscale = False):

      p.doscale = doscale
      p.is1D = is1D
      p.is2D = not is1D
      p.xlabel = xtext
      p.Other = [ "higgs", "Wjets", "triboson"]
      p.Zjets = ["Zeejets","Zttjets", "Zmmjets" ]
      p.Data = ["data18","data15", "data16", "data17"]

      p.Neg = False
      if "Neg" in hname:
        p.Neg = True
        print()


      if p.is1D: p.plot1D(hdic,hname,bkgs)
      else: p.plot2D(hdic,hname,bkgs)

    def plot2D(p,hdic,hname = "lepPt_ele_hT", bkgs = []):

      p.isEff = False

      p.nTotBkg = 0.0
      p.stackorder = []
      p.dyield = {}
      
      #p.doScale(hdic,hname,["WmuHNL50_30G"],1.0/80000.)
      # Legend
      p.leg = R.TLegend(0.60,0.56,0.91,0.90)
      p.leg.SetBorderSize(0)
      p.leg.SetTextSize(0.02)
      p.leg.SetNColumns(2)

        
      p.getYields(hdic,hname,bkgs)
      p.can  = R.TCanvas('c1','c1',1000,1000)
      #p.customise_gPad()

      if p.doscale:
        p.scaleToUnity(hdic,hname,bkgs)

      p.hstack = R.THStack()
      p.fillStack(hdic,hname,bkgs)

      histarray = p.hstack.GetStack()

      print("All histograms: %i"%histarray.GetEntries())

      for bkg in p.dyield.keys():
        if d_samp[bkg]["type"] == "bkg" and (p.dyield[bkg]/p.nTotBkg) < 0.05:
          for ha in histarray:
            if bkg in ha.GetName():
              histarray.Remove(ha)
              break

      print("Histograms to plot (>0.05 of total): %i"%histarray.GetEntries())

      p.can.Divide(3,3)
      p.can.Update()
      i = 1
      lines = [R.TLine(-1,-1,1,1) for i in range(len(histarray))]
      for ha in histarray:
        p.can.cd(i)
        #p.customise_gPad()
        p.customise_gPad(top=0.08, bot=0.08, left=0.1, right=0.15)
        ha.GetZaxis().SetRangeUser(0,1)
        ha.Draw("colz")
        key = findKey(ha.GetName())
        myText(0.37, 0.87, 'N(%s) = %.1f'%(key,p.dyield[key]), 0.05, R.kRed+2)
        myLine(-1,-1,1,1,R.kRed+2,8,2)
        i += 1

      
      #p.can.Update()

    
      

      #p.hstack.Draw("padscolz")
        
    def plot1D(p,hdic,hname = "lepPt_ele_hT", bkgs = []):

        p.isEff = False
        if "_EF_" in hname:
          p.isEff = True
      
        # Define canvas and pads
        p.can  = R.TCanvas('','',1000,1000)
        p.customise_gPad()
        if not p.isEff:
          if p.Neg:
            p.pad1 = R.TPad('pad1', '', 0.0, 0., 1.0, 1.0)
          else:
            p.pad1 = R.TPad('pad1', '', 0.0, 0.40, 1.0, 1.0)
            p.pad2 = R.TPad('pad2', '', 0.0, 0.00, 1.0, 0.4)

        # Margins used for the pads
        gpLeft = 0.17
        gpRight = 0.05
        #-------
        # PAD1
        #-------
        if not p.isEff:
          p.pad1.Draw()
          p.pad1.cd()
          if p.Neg:
            p.customise_gPad(top=0.08, bot=0.25, left=gpLeft, right=gpRight)
          else:
            p.customise_gPad(top=0.08, bot=0.04, left=gpLeft, right=gpRight)

        # Legend
        if not p.isEff:
          if p.Neg:
            p.leg = R.TLegend(0.60,0.71,0.91,0.90)
          else:
            p.leg = R.TLegend(0.60,0.56,0.91,0.90)
        else:
          p.leg = R.TLegend(0.55,0.77,0.91,0.94)
        p.leg.SetBorderSize(0)
        p.leg.SetNColumns(2)
        if p.Neg:
          p.leg.SetTextSize(0.015)
        else:
          p.leg.SetTextSize(0.02)


        p.nTotBkg = 0.0
        p.stackorder = []
        p.dyield = {}
        

        p.getYields(hdic,hname,bkgs)

        p.hstack = R.THStack()
        p.fillStack(hdic,hname,reversed(p.stackorder))

        
        p.datastack = R.THStack()
        p.getData(hdic,hname,bkgs)

        p.signalstack = R.THStack()
        p.getSignal(hdic,hname,bkgs)

        print("-->",p.hstack.GetNhists())
        if p.hstack.GetNhists() > 0:
          if not p.isEff: p.hstack.Draw("hist")
          else: p.hstack.Draw("nostack")
          if p.datastack.GetNhists() > 0:
            p.datastack.Draw("same ep")
          if p.signalstack.GetNhists() > 0:
            p.signalstack.Draw("nostack same hist")
        elif p.datastack.GetNhists() > 0:
          p.datastack.Draw("ep")
          if p.signalstack.GetNhists() > 0: p.signalstack.Draw("nostack same hist")
        elif p.signalstack.GetNhists() > 0:
          p.signalstack.Draw("nostack hist")
        else:
          print("Sorry there's nothing there to plot")
        
        p.leg.Draw()

        ATL_status = "Internal"
        text_size = 0.045


        myText(0.22, 0.87, '#bf{#it{ATLAS}} ' + ATL_status, text_size*1.2, R.kBlack)

       

        xtitle = hname
        ytitle = 'Events' if not p.isEff else 'Efficieny'
        IsLogY = True
        enlargeYaxis = False
        scaling = "False"

        if not p.isEff:
          p.pad1.SetLogy(IsLogY)
        

        try:
          maxbin = p.hstack.GetStack().Last().GetBinContent(p.hstack.GetStack().Last().GetMaximumBin())
          p.customise_axes(p.hstack, xtitle, ytitle, 1.1, IsLogY, enlargeYaxis, maxbin, scaling == 'True')
        except:
          maxbin = 0
          p.customise_axes(p.datastack, xtitle, ytitle, 1.1, IsLogY, enlargeYaxis, maxbin, scaling == 'True')
        
        if p.Neg:
          p.can.Update()
        

        if not p.isEff:
          myText(0.77, 0.47, 'N(Bkg) = %.0f'%(p.nTotBkg), 0.025, R.kBlack)

        
        p.ratio = R.TH1D()
        # if (p.hstack.GetStack().Last() != None or p.datastack.GetStack().Last() != None):
        if not p.isEff:
          try:
            p.getRatio(p.hstack.GetStack().Last(),p.datastack.GetStack().Last())
          except:
            print("Could not get ratio!")
        
        #-------
        # PAD2
        #-------
        p.can.cd()
        if not p.isEff and not p.Neg:
          p.pad2.Draw()
          p.pad2.cd()
          p.customise_gPad(top=0.05, bot=0.39, left=gpLeft, right=gpRight)
          #customise_gPad(top=0, bot=0.39, left=gpLeft, right=gpRight) # joins upper and lower plot
          p.pad2.SetGridy()
          p.ratio.Draw("e0p")

          xtitle = hname
          ytitle = 'Data / SM'
          IsLogY = True
          enlargeYaxis = False

          maxbin = p.ratio.GetBinContent(p.ratio.GetMaximumBin())

          p.customise_axes(p.ratio, xtitle, ytitle, 1.1, IsLogY, enlargeYaxis, maxbin, scaling == 'True')

          p.can.Update()

    def scaleToUnity(p,histo,hkey,procs):
      for k in procs:
        if not hkey+"_%s"%k in histo.keys():
          continue
        if p.is1D:
          histo[hkey+"_%s"%k].Scale(1./(histo[hkey+"_%s"%k].Integral(0,histo[hkey+"_%s"%k].GetNbinsX()+1)))
        else:
          histo[hkey+"_%s"%k].Scale(1./histo[hkey+"_%s"%k].Integral(0,histo[hkey+"_%s"%k].GetNbinsX()+1,0,histo[hkey+"_%s"%k].GetNbinsY()+1))

    def getYields(p,histo,hkey,procs):
        for k in procs:
            if d_samp[k]["type"] == "sig": continue
            newkey = hkey+"_%s"%k
            if p.isEff:
              newkey = hkey.replace("_EF_","_SG_")+"_%s"%k
            if not newkey in histo.keys():
                continue
            if not (k in p.Other or k in p.Data or k in p.Zjets):
              histo_i = histo[hkey+"_%s"%k]
            elif k == p.Other[0]:
              histo_i, _ = p.mergeChannels(histo,p.Other,hkey,k=k)
            elif k == p.Data[0]:
              histo_i, _ = p.mergeChannels(histo,p.Data,hkey,k=k)
            elif k == p.Zjets[0]:
              histo_i, _ = p.mergeChannels(histo,p.Zjets,hkey,k=k)
            else:
              continue 
            if p.is1D:
              p.dyield[k] = histo_i.Integral(0,histo_i.GetNbinsX()+1)
            else:
              p.dyield[k] = histo_i.Integral(0,histo_i.GetNbinsX()+1,0,histo_i.GetNbinsY()+1)
            if d_samp[k]["type"] == "bkg":
                p.nTotBkg += p.dyield[k]
        newdict = p.dyield.copy()
        while True:
            maxi = -99999
            maxkey = ""
            for k in newdict.keys():
                if not d_samp[k]["type"] == "bkg": continue
                if newdict[k] > maxi:
                    maxkey = k
                    maxi = newdict[k]
            if not maxkey: break
            if maxkey in p.stackorder:
                break
            p.stackorder.append(maxkey)
            ret = newdict.pop(maxkey,None)
            if ret == None:
                break

    def doScale(p,histo,hkey,procs,fac = 1.0):
        for k in procs:
            if d_samp[k]["type"] == "data": continue
            histo[hkey+"_%s"%k].Scale(fac)

    def getRatio(p,h1,h2):
        p.ratio = h2.Clone("hRatio")
        p.ratio.Sumw2() 
        p.ratio.Divide(h1) 
        #p.ratio.SetLineColor(R.kGray+2)
        p.ratio.SetLineWidth(2)
        p.ratio.SetMarkerStyle(21)
      
    def mergeChannels(p,histo,group,hkey,name = None,pc_yield = None,k = None):
      histo_i = histo[hkey+"_%s"%group[0]].Clone(hkey+"_%s_SUM"%k)
      for i in range(1,len(group)):
        histo_i.Add(histo[hkey+"_%s"%group[i]].GetPtr())
      if pc_yield is None:
        return histo_i, None
      else:
        leg_txt = '{0} ({1:.1f}%)'.format(name, pc_yield)
        return histo_i, leg_txt


        
    def fillStack(p,histo,hkey,procs):
        for k in procs:
            if p.is1D and not d_samp[k]["type"] == "bkg": continue
            if not hkey+"_%s"%k in histo.keys():
                print("Could not find key %s in histo dctionary"%(hkey+"_%s"%k ))
                continue
            
            if float(p.nTotBkg) != 0: pc_yield = 100 * ( p.dyield[k] / float(p.nTotBkg) )
            if not p.isEff:
              leg_txt = '{0} ({1:.1f}%)'.format( d_samp[k]["leg"], pc_yield )
            else:
              print(k, pc_yield)
              if pc_yield < 10: continue
              leg_txt = '{0}'.format( d_samp[k]["leg"] )
            if not (k in p.Other or k in p.Data or k in p.Zjets):
              histo_i = histo[hkey+"_%s"%k]
            elif k == p.Other[0]:
              histo_i, leg_txt = p.mergeChannels(histo,p.Other,hkey,"Other",pc_yield,k)
            elif k == p.Zjets[0]:
              histo_i, leg_txt = p.mergeChannels(histo,p.Zjets,hkey,"Zjets",pc_yield,k)
            else:
              continue  
            try:
              p.hstack.Add(histo_i)
              p.leg.AddEntry(histo_i,leg_txt,"lpf")
            except:
              p.hstack.Add(histo_i.GetValue())
              p.leg.AddEntry(histo_i.GetValue(),leg_txt,"lpf")

                
    def getData(p,histo,hkey,procs):
      k = p.Data[0]
      histo_i, _ = p.mergeChannels(histo,p.Data,hkey,"Data",k = k)
      pc_yield = histo_i.Integral(0,histo_i.GetNbinsX()+1)

      if not p.isEff:
        leg_txt = '{0} ({1:.0f} Events)'.format("Data", pc_yield)
      else:
        leg_txt = '{0})'.format(d_samp[k]["leg"])
      try:
          p.datastack.Add(histo_i)
          p.leg.AddEntry(histo_i,leg_txt,"lp")
      except:
          p.datastack.Add(histo_i.GetValue())
          p.leg.AddEntry(histo_i.GetValue(),leg_txt,"lp")

    def getSignal(p,histo,hkey,procs):
      i = 1
      for k in procs:
        if not d_samp[k]["type"] == "sig": continue
        if not hkey+"_%s"%k in histo.keys():
            continue
        if i:
          histo_i = histo[hkey+"_%s"%k].Clone(hkey+"_%s_SUM"%k)
          i = 0
          continue
        histo_i.Add(histo[hkey+"_%s"%k].GetPtr())
      if not i:
        pc_yield = histo_i.Integral(0,histo_i.GetNbinsX()+1)
        leg_txt =  '{0} ({1:.0f})'.format("Signal", pc_yield)
        try:
          p.signalstack.Add(histo_i)
          p.leg.AddEntry(histo_i,leg_txt,"lp")
        except:
          p.signalstack.Add(histo_i.GetValue())
          p.leg.AddEntry(histo_i.GetValue(),leg_txt,"lp")
        p.signalstack.SetLineWidth(3)
    # Function for customising the gPad (gPad points to the current pad, and one can use gPad to set attributes of the current pad)

    def customise_gPad(p,top=0.03, bot=0.15, left=0.17, right=0.08):

        R.gPad.Update()

        R.gStyle.SetTitleFontSize(0.0)

        # gPad margins
        R.gPad.SetTopMargin(top)
        R.gPad.SetBottomMargin(bot)
        R.gPad.SetLeftMargin(left)
        R.gPad.SetRightMargin(right)

        R.gStyle.SetOptStat(0) # Hide usual stats box

        R.gPad.Update()

    # Funcion for customising axes
    def customise_axes(p,hist, xtitle, ytitle, scaleFactor=1.1, IsLogY=False, enlargeYaxis=False, maxbin = 10, scaling=False):

        # Set a universal text size
        text_size = 45

        R.TGaxis.SetMaxDigits(4)

        ##################################
        # X axis
        xax = hist.GetXaxis()

        # Precision 3 Helvetica (specify label size in pixels)
        xax.SetLabelFont(43)
        xax.SetTitleFont(43)
        # xax.SetTitleFont(13) # times

        xax.SetTitle(p.xlabel)
        xax.SetTitleSize(text_size)

        print("ytitle  = ",  ytitle)

        # Top panel
        if 'Events' in ytitle and not p.Neg:
            xax.SetLabelSize(0)
            xax.SetLabelOffset(0.02)
            xax.SetTitleOffset(2.0)
            xax.SetTickSize(0.04)
        # Bottom panel
        else:
            xax.SetLabelSize(text_size - 7)
            xax.SetLabelOffset(0.03)
            if p.Neg:
              xax.SetTitleOffset(1.5)
            else:
              xax.SetTitleOffset(3.5)
            xax.SetTickSize(0.08)

        # xax.SetRangeUser(0,2000)
        # xax.SetNdivisions(-505)

        R.gPad.SetTickx()

        ##################################
        # Y axis
        yax = hist.GetYaxis()

        # Precision 3 Helvetica (specify label size in pixels)
        yax.SetLabelFont(43)
        yax.SetTitleFont(43)

        yax.SetTitle(ytitle)
        yax.SetTitleSize(text_size)
        yax.SetTitleOffset(1.8)

        yax.SetLabelOffset(0.015)
        yax.SetLabelSize(text_size - 7)

        ymax = hist.GetMaximum()
        ymin = hist.GetMinimum()

        # print(yax.GetTitle()
        # print('ymax       = ', ymax)
        # print('SF         = ', scaleFactor )
        # print('ymax x SF  = ', ymax*scaleFactor)
        # print('ymin       = ', ymin)
        # print('ymin x 0.9 = ', ymin*scaleFactor)
        # if ymin == 0.0:
        #     print 'ymin = 0.0'

        # Top events panel
        if 'Events' in ytitle:
            yax.SetNdivisions(505)
            if IsLogY:
                if enlargeYaxis:
                    ymax *= 2 * 10 ** 10
                    ymin = 0.01
                else:
                    # ymax = 3 * 10 ** 4
                    # ymin = 0.5
                    ymax *= 10*10
                    if p.Neg:
                      ymin = 0.05
                      ymax = 1e6
                    else:
                      ymin = 5
                #if scaling:
                #    hist.SetMaximum(1.0)
                #else:
                hist.SetMaximum(ymax)
                hist.SetMinimum(ymin)
            else:
                #if scaling:
                #    hist.SetMaximum(1.0)
                #else:
                hist.SetMaximum(ymax*scaleFactor)
                hist.SetMinimum(0.0)
        elif 'Efficiency' in ytitle:
          yax.SetNdivisions(505)
          hist.SetMaximum(1.2)
          hist.SetMinimum(0.0)
        # Bottom panel
        elif 'Ratio' or 'Data / SM'  in ytitle:
            yax.SetNdivisions(505)
            # Dynamic 
            if ymax*scaleFactor > 5:
                hist.SetMaximum(5)
            else: hist.SetMaximum(ymax*scaleFactor)
            if ymin*0.9 < -1:
                hist.SetMinimum(-1)#ymin*0.9)
            else: hist.SetMinimum(ymin*0.9)
            # Fixed
            #hist.SetMinimum(-0.5) 
            #hist.SetMaximum(2.5)  

        R.gPad.SetTicky()

        R.gPad.Update()
    

        
