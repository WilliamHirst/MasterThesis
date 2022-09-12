import os
import sys

yr = int(sys.argv[1])

submit = int(sys.argv[2])
analy = sys.argv[3]

if analy == "FCL":
    if yr == 2018:
        #indir = "/scratch/eirikgr/ANAoutput/Thu_May_20_2021_02L_NTUP_data18_R21/"
        #indir = "/scratch/eirikgr/ANAoutput/Thu_May_20_2021_02L_NTUP_data18_REL22/"#
        #indir = "/scratch/eirikgr/ANAoutput/Wed_May_19_2021_02L_NTUP_data18_REL22/"#
        #indir = "/scratch/eirikgr/ANAoutput/Wed_May_19_2021_02L_NTUP_data18_R21/"#
        #indir = "/scratch/eirikgr/ANAoutput/Wed_Jun_02_2021_02L_NTUP_data18_R21/"
        indir = "/scratch/eirikgr/ANAoutput/Tue_Jun_08_2021_02L_NTUP_data18_R21/"
        #Tue_May_18_2021_02L_NTUP_data18_R21"#Wed_May_12_2021_02L_NTUP_data18_R21/"#"/scratch/eirikgr/ANAoutput/Wed_Feb_03_2021_02L_NTUP_data18_EST2/"
    elif yr == 2015 or yr == 2016:
        #indir = "/scratch/eirikgr/ANAoutput/Wed_Feb_03_2021_02L_NTUP_data1516_EST2/"
        indir = "/scratch/eirikgr/ANAoutput/Tue_Jun_08_2021_02L_NTUP_data1516_R21/"
    elif yr == 2017:
        #indir = "/scratch/eirikgr/ANAoutput/Wed_Feb_03_2021_02L_NTUP_data17_EST2/"
        indir = "/scratch/eirikgr/ANAoutput/Tue_Jun_08_2021_02L_NTUP_data17_R21/"
elif analy == "2L2J":
    if yr == 2018:
        indir = "/scratch/eirikgr/ANAoutput/Fri_Jul_09_2021_02L_NTUP_data18_OLD2L2J/"
    elif yr == 2015 or yr == 2016:
        indir = "/scratch/eirikgr/ANAoutput/Mon_Feb_08_2021_02L_NTUP_data1516_EST2/"
    elif yr == 2017:
        indir = "/scratch/eirikgr/ANAoutput/Mon_Feb_08_2021_02L_NTUP_data17_EST2/"
elif analy == "ZMET":
    if yr == 2018:
        indir = "/scratch/eirikgr/ANAoutput/Mon_Jul_25_2022_02L_NTUP_year18_ZMET/"
    elif yr == 2015 or yr == 2016:
        indir = "/scratch/eirikgr/ANAoutput/Mon_Jul_25_2022_02L_NTUP_year1516_ZMET/."
    elif yr == 2017:
        indir = "/scratch/eirikgr/ANAoutput/Mon_Jul_25_2022_02L_NTUP_year17_ZMET/."
leptons = ["el","mu"]
regions_to_use = ["real1d","light","heavy","conv"]#,"real1d","real1deta"]#"#"light","heavy","conv",
regions_to_use = ["real","heavy"]#,"real1d","real1deta"]#"#"light","heavy","conv",

#regions_to_use = ["real"]
uncert = []
uncert.append("\"\"")
'''
uncert.append("EL_EFF_ID_TOTAL_1down")            
uncert.append("EL_EFF_ID_TOTAL_1up")              
uncert.append("EL_EFF_Iso_TOTAL_1down")           
uncert.append("EL_EFF_Iso_TOTAL_1up")             
uncert.append("EL_EFF_Reco_TOTAL_1down")          
uncert.append("EL_EFF_Reco_TOTAL_1up")            
uncert.append("MUON_EFF_ISO_STAT_1down")          
uncert.append("MUON_EFF_ISO_STAT_1up")            
uncert.append("MUON_EFF_ISO_SYS_1down")           
uncert.append("MUON_EFF_ISO_SYS_1up")             
uncert.append("MUON_EFF_RECO_STAT_1down")         
uncert.append("MUON_EFF_RECO_STAT_1up")           
uncert.append("MUON_EFF_RECO_STAT_LOWPT_1down")   
uncert.append("MUON_EFF_RECO_STAT_LOWPT_1up")     
uncert.append("MUON_EFF_RECO_SYS_1down")          
uncert.append("MUON_EFF_RECO_SYS_1up")            
uncert.append("MUON_EFF_RECO_SYS_LOWPT_1down")    
uncert.append("MUON_EFF_RECO_SYS_LOWPT_1up")      
uncert.append("EL_EFF_Trigger_TOTAL_1down")       
uncert.append("EL_EFF_Trigger_TOTAL_1up")         
uncert.append("MUON_EFF_TrigStat_1down")          
uncert.append("MUON_EFF_TrigStat_1up")            
uncert.append("MUON_EFF_TrigSyst_1down")          
uncert.append("MUON_EFF_TrigSyst_1up")    
'''
for l in leptons:
    for reg in regions_to_use:
        if l == "mu" and reg in ["conv","light"]: continue
        for u in uncert:
            if ("EL" in u or "MUON" in u) and reg in ["conv","heavy"]: continue
            if l == "el" and "MUON" in u: continue
            if l == "mu" and "EL" in u: continue
            cmd = "python makeMMinput.py %s %s %s %s %s" %(indir,reg,l,u,analy)
            if submit: os.system(cmd)
            print(cmd)
