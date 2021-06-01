
## Data analysis for InstaScope data with respect to time - Master File ##

## SET YOUR OWN PATH WHERE FILES TO BE ANALYZED ARE
getwd()
setwd("set_working_directory")
#Make sure path is correct
getwd()

#Get data from CSV files 
#library(readr)
raw_dataset <- read.csv("your_file.csv", skip = 38)

#Check that your top rows of the dataset are OK
head(raw_dataset)

#Define Parameters
###Define IS Noise thresholds using Force trigger data (3 standard deviations)
FL1 = 255.6947324
FL2 = 82.92607225
FL3 = 18.33222312
  
#Pollen size limits (um)
LEpollen = 1.2
HEpollen = 10

#Bacteria size limits (um)
LEbacteria = 0.5
HEbacteria = 1.9

#Fungi size limits (um)
LEfungi = 2
HEfungi = 10


#Make it a data table object
library(data.table)
library('plyr')
raw_table <- data.table(raw_dataset)
raw_table_minutes <- raw_table[, ':='(Time_minutes = round(raw_table$Time/(1000*60)))]
head(raw_table_minutes)

#Filter data for Biological Particles (all particles that ARE fluorescent) above 0.5 um
fluorescent_particles <- raw_table_minutes[FL1_280>=FL1 & FL2_280>=FL2 & FL2_370>=FL3 & Size >= 0.5]

#Filter data for Non-biological Particles (all particles that are NOT fluorescent) above 0.5 um
non_fluorescent_particles <- raw_table_minutes[FL1_280<FL1 & FL2_280<FL2 & FL2_370<FL3 & Size >= 0.5]

#Filter data for all types of particles (total particles) above 0.5 um (instrument threshold)
all_particles <- raw_table_minutes[Size>=0.5]

#Filter data for Pollen Particles
ABC_pollen <- raw_table_minutes[FL1_280>=FL1 & FL2_280>=FL2 & FL2_370>=FL3 & Size>=LEpollen & Size<=HEpollen]
freq_ABC_pollen_table <- setDT(count(ABC_pollen$Time_minutes))
ABC_pollen_percentage <- freq_ABC_pollen_table[,':='(Pollen_estimate = freq_ABC_pollen_table$freq *0.5)]

BC_pollen <- raw_table_minutes[FL1_280<FL1 & FL2_280>=FL2 & FL2_370>=FL3 & Size>=LEpollen & Size<=HEpollen]
freq_BC_pollen_table <- setDT(count(BC_pollen$Time_minutes))
BC_pollen_percentage <- freq_BC_pollen_table[,':='(Pollen_estimate = freq_BC_pollen_table$freq *0.34)]

A_pollen <- raw_table_minutes[FL1_280>=FL1 & FL2_280<FL2 & FL2_370<FL3 & Size>=LEpollen & Size<=HEpollen]
freq_A_pollen_table <- setDT(count(A_pollen$Time_minutes))
A_pollen_percentage <- freq_A_pollen_table[,':='(Pollen_estimate = freq_A_pollen_table$freq *0.06)]

B_pollen <- raw_table_minutes[FL1_280<FL1 & FL2_280>=FL2 & FL2_370<FL3 & Size>=LEpollen & Size<=HEpollen]
freq_B_pollen_table <- setDT(count(B_pollen$Time_minutes))
B_pollen_percentage <- freq_B_pollen_table[,':='(Pollen_estimate = freq_B_pollen_table$freq *0.04)]

C_pollen <- raw_table_minutes[FL1_280<FL1 & FL2_280<FL2 & FL2_370>=FL3 & Size>=LEpollen & Size<=HEpollen]
freq_C_pollen_table <- setDT(count(C_pollen$Time_minutes))
C_pollen_percentage <- freq_C_pollen_table[,':='(Pollen_estimate = freq_C_pollen_table$freq *0.06)]

pollen_particles <- rbind(ABC_pollen_percentage,BC_pollen_percentage,A_pollen_percentage,B_pollen_percentage,C_pollen_percentage)

#Filter data for Fungal Particles 
AB_fungi <- raw_table_minutes[FL1_280>=FL1 & FL2_280>=FL2 & FL2_370<FL3 & Size>=LEfungi & Size<HEfungi]
freq_AB_fungi_table <- setDT(count(AB_fungi$Time_minutes))
AB_fungi_percentage <- freq_AB_fungi_table[,':='(Fungi_estimate = freq_AB_fungi_table$freq *0.15)]

A_fungi <- raw_table_minutes[FL1_280>=FL1 & FL2_280<FL2 & FL2_370<FL3 & Size>=LEfungi & Size<HEfungi]
freq_A_fungi_table <- setDT(count(A_fungi$Time_minutes))
A_fungi_percentage <- freq_A_fungi_table[,':='(Fungi_estimate = freq_A_fungi_table$freq *0.8)]

ABC_fungi <- raw_table_minutes[FL1_280>=FL1 & FL2_280>=FL2 & FL2_370>=FL3 & Size>=LEfungi & Size<HEfungi]
freq_ABC_fungi_table <- setDT(count(ABC_fungi$Time_minutes))
ABC_fungi_percentage <- freq_ABC_fungi_table[,':='(Fungi_estimate = freq_ABC_fungi_table$freq *0.05)]

fungal_particles <- rbind(A_fungi_percentage, AB_fungi_percentage, ABC_fungi_percentage)

#Filter data for Bacterial Particles
A_bacteria <- raw_table_minutes[FL1_280>=FL1 & FL2_280<FL2 & FL2_370<FL3 & Size>=LEbacteria & Size<HEbacteria]
freq_A_bacteria_table <- setDT(count(A_bacteria$Time_minutes))
A_bacteria_percentage <- freq_A_bacteria_table[,':='(Bacteria_estimate = freq_A_bacteria_table$freq *0.91)]

AB_bacteria <- raw_table_minutes[FL1_280>=FL1 & FL2_280>=FL2 & FL2_370<FL3 & Size>=LEbacteria & Size<HEbacteria]
freq_AB_bacteria_table <- setDT(count(AB_bacteria$Time_minutes))
AB_bacteria_percentage <- freq_AB_bacteria_table[,':='(Bacteria_estimate = freq_AB_bacteria_table$freq *0.09)]

bacterial_particles <- rbind(A_bacteria_percentage, AB_bacteria_percentage)

#Consolidate those tables!
pollen_after_percentages <- as.data.table(pollen_particles)[,lapply(.SD, sum), by = .(x=tolower(x))]
fungi_after_percentages <- as.data.table(fungal_particles)[,lapply(.SD, sum), by = .(x=tolower(x))]
bacteria_after_percentages <- as.data.table(bacterial_particles)[,lapply(.SD, sum), by = .(x=tolower(x))]
freq_fluorescent_table <- count(fluorescent_particles$Time_minutes)
freq_nonfluorescent_table <- count(non_fluorescent_particles$Time_minutes)
freq_allparticles_table <- count(all_particles$Time_minutes)

#### Include missing particles ####
missing_particles <- sum(raw_table_minutes$TPCT2)
counted_particles <- nrow(all_particles)
total_partices = missing_particles + counted_particles
ratio_missing_to_total = missing_particles/total_partices

#### Calculate missing ratios for each fluorescence type ####
#Pollen
pollen_total_value = sum(pollen_after_percentages$Pollen_estimate)
ratio_pollen_measured = pollen_total_value/counted_particles
pollen_missed = ratio_pollen_measured * missing_particles
ratio_pollen_missed = pollen_missed/total_partices

#Fungi
fungi_total_value = sum(fungi_after_percentages$Fungi_estimate)
ratio_fungi_measured = fungi_total_value/counted_particles
fungi_missed = ratio_fungi_measured * missing_particles
ratio_fungi_missed = fungi_missed/total_partices


#Bacteria
bacteria_total_value = sum(bacteria_after_percentages$Bacteria_estimate)
ratio_bacteria_measured = bacteria_total_value/counted_particles
bacteria_missed = ratio_bacteria_measured * missing_particles
ratio_bacteria_missed = bacteria_missed/total_partices

#Fluorescent
fluorescent_total_value = nrow(fluorescent_particles)
ratio_fluorescent_measured = fluorescent_total_value/counted_particles
fluorescent_missed = ratio_fluorescent_measured * missing_particles
ratio_fluorescent_missed = fluorescent_missed/total_partices


#Non Fluorescent
nonfluorescent_total_value = nrow(non_fluorescent_particles)
ratio_nonfluorescent_measured = nonfluorescent_total_value/counted_particles
nonfluorescent_missed = ratio_nonfluorescent_measured * missing_particles
ratio_nonfluorescent_missed = nonfluorescent_missed/total_partices

#New values when including missing particles
pollen_totals <- pollen_after_percentages[, ':='(Pollen_total = pollen_after_percentages$Pollen_estimate*(1+ratio_pollen_missed))]
fungi_totals <- fungi_after_percentages[, ':='(Fungi_total = fungi_after_percentages$Fungi_estimate*(1+ratio_fungi_missed))]
bacteria_totals <- bacteria_after_percentages[, ':='(Bacteria_total = bacteria_after_percentages$Bacteria_estimate*(1+ratio_bacteria_missed))]
setDT(freq_fluorescent_table)
fluorescent_totals <- freq_fluorescent_table[, ':='(Fluorescent_total = freq_fluorescent_table$freq*(1+ratio_fluorescent_missed))]
setDT(freq_nonfluorescent_table)
nonfluorescent_totals <- freq_nonfluorescent_table[, ':='(Nonfluorescent_total = freq_nonfluorescent_table$freq*(1+ratio_nonfluorescent_missed))]
setDT(freq_allparticles_table)
total_particles_table <- freq_allparticles_table[, ':='(All_particles_total = freq_allparticles_table$freq*(1+ratio_missing_to_total))]



#Plot those particles (this is IN LOG SCALE, but it could be linear as well by removing "log="y""). Can also change Y axes values with ylim (first value before coma is minimum, second is maximum)
plot(total_particles_table$x, total_particles_table$All_particles_total, col='black', pch=19, ylim=c(1,11000),xlab = 'Time (minutes)', ylab = 'Particles',log="y")
points(fluorescent_totals$x, fluorescent_totals$Fluorescent_total, col='cyan', pch=19)
points(pollen_totals$x, pollen_totals$Pollen_total, col='blue', pch=19)
points(fungi_totals$x, fungi_totals$Fungi_total, col='green',pch=19)
points(bacteria_totals$x,  bacteria_totals$Bacteria_total, col='red', pch=19)
points(nonfluorescent_totals$x, nonfluorescent_totals$Nonfluorescent_total, col='magenta',pch=19)
legend('topleft',c("All Particles","Fluorescent","Pollen","Fungi","Bacteria","Non Fluorescent"), fill=c("black","cyan","blue","green","red","magenta"))



##CSV file preparation and generation
#Create frequency tables for each particle type
library('plyr')
#Pollen
Time_pollen <- pollen_totals
colnames(Time_pollen) <- c("Time (minutes)", "Frequency_Raw", "Pollen_estimate_counted", "Pollen_Total")
head(Time_pollen)

#Fungi
Time_fungi <- fungi_totals
colnames(Time_fungi) <- c("Time (minutes)", "Frequency_Raw", "Fungi_estimate_counted", "Fungi_Total")
head(Time_fungi)

#Bacteria
Time_bacteria <- bacteria_totals
colnames(Time_bacteria) <- c("Time (minutes)", "Frequency_Raw", "Bacteria_estimate_counted", "Bacteria_Total")
head(Time_bacteria)

#Fluorescent
Time_fluorescent <- fluorescent_totals
colnames(Time_fluorescent) <- c("Time (minutes)", "Frequency_Raw", "Fluorescent_Total")
head(Time_fluorescent)

#Non Fluorescent
Time_nonfluorescent <- nonfluorescent_totals
colnames(Time_nonfluorescent) <- c("Time (minutes)", "Frequency_Raw", "Nonfluorescent_Total")
head(Time_nonfluorescent)
 
#Generate separate CSVs for each particle type, so you can open them with Excel
write.csv(Time_pollen, file='Pollen_wrt_time.csv')
write.csv(Time_fungi, file='Fungi_wrt_time.csv')
write.csv(Time_bacteria, file='Bacteria_wrt_time.csv')
write.csv(Time_fluorescent, file='Fluorescent_wrt_time.csv')
write.csv(Time_nonfluorescent, file='Nonfluorescent_wrt_time.csv')
