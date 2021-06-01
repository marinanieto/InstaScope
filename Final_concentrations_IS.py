#Define Thresholds for FL1, FL2 and FL3
FL1= 255.6947324
FL2= 82.92607225
FL3= 18.33222312

#Define collection time in minutes (hr*60min/1h)
time=425.65

#Define LowEnd Size [LE] and HighEnd Size [HE]
LEpollen=1.2
HEpollen=10

LEbacteria=0.4
HEbacteria=1.9

LEfungi=2
HEfungi=10

#Define IS Flow Rate (L/min)
flowrate=0.835

import pandas as pd
import glob
import os
import numpy as np


ROWS_TO_SKIP = 38
#Define name of csv files
CSV_PREFIX = "bsubtilis_"
COL_TO_UPDATE = "Time"

if __name__ == "__main__":
    merged = pd.DataFrame()
    last_time = 0
    for csv in sorted(glob.glob(CSV_PREFIX + "*.csv")):
        df = pd.read_csv(csv, skiprows=ROWS_TO_SKIP)
        merged = merged.append(df)
        merged[COL_TO_UPDATE] += last_time
        last_time = int(merged.iloc[-1][COL_TO_UPDATE])
    merged.to_csv("merged.csv", index=None)



#Combine all CSV files in a concatenated CSV
combined_csv=pd.read_csv('merged.csv')

#Show combined CSV
combined_csv

#Total number of nonfluorescent particles
total_particle = combined_csv.shape[0] + 1

#Calculate all pollen particles in each Fluorescence Channel FL based on Mark's paper Graph (correction factor INCLUDED)
ABC_pollen=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370>=FL3) & (combined_csv.Size>=LEpollen) & (combined_csv.Size<=HEpollen)]
BC_pollen=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370>=FL3) & (combined_csv.Size>=LEpollen) & (combined_csv.Size<=HEpollen)]
A_pollen=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370>FL3) & (combined_csv.Size>=LEpollen) & (combined_csv.Size<=HEpollen)]
B_pollen=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370>FL3) & (combined_csv.Size>=LEpollen) & (combined_csv.Size<=HEpollen)]
C_pollen=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370>=FL3) & (combined_csv.Size>=LEpollen) & (combined_csv.Size<=HEpollen)]

#Merge all pollen particles together    
pollen_merged=pd.concat([0.05*A_pollen,0.02*B_pollen,0.05*C_pollen,0.34*BC_pollen,0.5*ABC_pollen], ignore_index=True)
#Number of pollen particles
total_pollen = pollen_merged.shape[0] + 1

#Calculate all bacteria particles in each Fluorescence Channel FL based on Mark's paper (correction factor INCLUDED)
A_bacteria=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370<FL3) & (combined_csv.Size>=LEbacteria) & (combined_csv.Size<=HEbacteria)]
AB_bacteria=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370<FL3) & (combined_csv.Size>=LEbacteria) & (combined_csv.Size<=HEbacteria)]
#Merge all bacteria particles together
bacteria_merged=pd.concat([0.91*A_bacteria,0.09*AB_bacteria], ignore_index=True)  
#Number of bacteria particles
total_bacteria = bacteria_merged.shape[0] + 1   

#Calculate all fungi particles in each Fluorescence Channel FL based on Mark's paper (correction factors INCLUDED)
ABC_fungi=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370>=FL3) & (combined_csv.Size>=LEfungi) & (combined_csv.Size<=HEfungi)]
AB_fungi=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370<FL3) & (combined_csv.Size>=LEfungi) & (combined_csv.Size<=HEfungi)]
A_fungi=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370<FL3) & (combined_csv.Size>=LEfungi) & (combined_csv.Size<=HEfungi)]
#Merge all fungi together
fungi_merged=pd.concat([0.80*A_fungi,0.15*AB_fungi,0.05*ABC_fungi], ignore_index=True)
#Number of fungi particles
total_fungi = fungi_merged.shape[0] + 1

#Calculate all fluorescent particles in file with each of the 7 types of fluorescence. No size gating.
ABC_fluorescent=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370>=FL3)]
AB_fluorescent=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370<FL3)]
A_fluorescent=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370<FL3)]
B_fluorescent=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370<FL3)]
C_fluorescent=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370>=FL3)]
AC_fluorescent=combined_csv.loc[(combined_csv.FL1_280>=FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370>=FL3)]
BC_fluorescent=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280>=FL2) & (combined_csv.FL2_370>=FL3)]
#Combine all fluorescent particles in one data frame
all_fluorescent=pd.concat([A_fluorescent,B_fluorescent,C_fluorescent,AB_fluorescent,AC_fluorescent,BC_fluorescent,ABC_fluorescent], ignore_index=True)
#Number of fluorescent particles
total_fluorescent=all_fluorescent.shape[0] + 1

#Total nonfluorescent
all_nonfluorescent=combined_csv.loc[(combined_csv.FL1_280<FL1) & (combined_csv.FL2_280<FL2) & (combined_csv.FL2_370<FL3)]
#Number of fluorescent particles
total_nonfluorescent=all_nonfluorescent.shape[0] + 1

#Calculate concentrations
Conc_pollen=total_pollen/time/flowrate
Conc_bacteria=total_bacteria/time/flowrate
Conc_fungi=total_fungi/time/flowrate
Conc_fluorescent=total_fluorescent/time/flowrate
Conc_nonfluorescent=total_nonfluorescent/time/flowrate
Conc_total=total_particle/time/flowrate

#Now we need to include the missing particles
ALLparticles=combined_csv['TPCT2'].sum()
missingparticles=ALLparticles-total_particle
ratio_allmissingparticles=missingparticles/total_particle
miss_totalpollen=total_pollen+total_pollen/total_particle*missingparticles

miss_totalbacteria=total_bacteria+total_bacteria/total_particle*missingparticles
miss_totalfungi=total_fungi+total_fungi/total_particle*missingparticles
miss_totalfluoresc=total_fluorescent+total_fluorescent/total_particle*missingparticles
miss_totalnonfluoresc=total_nonfluorescent+total_nonfluorescent/total_particle*missingparticles

#Calculate concentrations including missing particles
MissedConc_pollen=miss_totalpollen/time/flowrate
MissedConc_bacteria=miss_totalbacteria/time/flowrate
MissedConc_fungi=miss_totalfungi/time/flowrate
MissedConc_fluorescent=miss_totalfluoresc/time/flowrate
MissedConc_nonfluorescent=miss_totalnonfluoresc/time/flowrate
MissedConc_total=ALLparticles/time/flowrate

#Create new dataframe with important values
Columns={'Particle #':[total_pollen, total_bacteria, total_fungi, total_fluorescent, total_nonfluorescent, total_particle], 'Concentration (L/min)':[Conc_pollen,Conc_bacteria,Conc_fungi, Conc_fluorescent, Conc_nonfluorescent, Conc_total], 'Particle # (including missed)':[miss_totalpollen,miss_totalbacteria,miss_totalfungi, miss_totalfluoresc,miss_totalnonfluoresc,ALLparticles],'Concentration including missed particles(L/min)':[MissedConc_pollen,MissedConc_bacteria,MissedConc_fungi,MissedConc_fluorescent,MissedConc_nonfluorescent,MissedConc_total]}
Table=pd.DataFrame(Columns, index=('Pollen','Bacteria','Fungi','Fluorescent','Non Fluorescent','Total Particle'))
Rounded_table=Table.round(3)
#Create a CSV with new table
Rounded_table.to_csv('Results.csv')




