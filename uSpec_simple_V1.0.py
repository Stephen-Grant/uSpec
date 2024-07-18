#uSPEC Data Reading and Processing
#S.Grant


#Import Libraries
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.integrate as integrate

#Today's date, the date you are processing data, not when it was collected!
date = datetime.now().strftime("%Y%m%d_%H%M%S")

#calculating energy of photons
c = 299792458; # speed of light, units m/s
h = 6.62607015E-34; # Planck's constant, units J*s

#Open the Log file and import all data in 
with open(r"C:\Users\stepheg\Downloads\uSpec_240415-16\LOG_0082.txt") as f: #Enter path for LOG file here <----
   lis = [line.split() for line in f]        
   wl=(lis[37][2])                                                                  #Extract Wavelengths and convert to array
   wlarray=(np.array(wl.split(",")))                                        
   wlarray[0]=wlarray[0].translate({ord('{'): None})
   wlarray[268]=wlarray[268].translate({ord('}'): None})
   wlarray=wlarray.astype(np.float64)

   #Empty variables for holding data, pre allocated for speed
   spectra_energy=[];spectra_photons=[];specdate=[];spectime=[];
   specsig=[];specinttime=[];specdepth=[];Eparlist=[];Epar_photons_list=[];timestamp=[];

   E_photon = ((h*c)/(wlarray))/1E-9; # Energy of photons, units J/photon
   
   #Saves spectra, date, integration time, signal %, depth.
   #Save integrated Epar in W/m^2
   for i in range(54,(len(lis)-1)):

      spec=np.array(lis[i][0].split(","))

      specdate.append(spec[2]);
      dati=spec[2]+" "+spec[3]
      timestamp.append(dati)
      spectime.append(spec[3]);

      specsig.append(spec[10]);

      specinttime.append(float(spec[9]))

      specdepth.append(float(spec[5]))

      spectrum=spec[11:280].astype(np.float64)
      spectrum_W=spectrum*1E-6;                                #Convert from uW to W
      spectrum_W_m2=spectrum_W*10000;                 #Convert from cm^-2 to m^-2
      spectrum_W_m2=spectrum*1E-2;
      spectrum_photons =spectrum_W_m2/E_photon; # number of photons, units photons/m^2/nm/s
      spectrum_photons_umol=spectrum_photons/6.022E17;
 
      
      Epar_energy=(integrate.trapezoid(spectrum_W_m2[15:146],wlarray[15:146]));
      Epar_photons=(integrate.trapezoid(spectrum_photons_umol[15:146],wlarray[15:146]));
      Eparlist.append(Epar_energy)
      Epar_photons_list.append(Epar_photons)


      spectra_energy.append(spectrum_W_m2)
      spectra_photons.append(spectrum_photons)
      

      #Spectra Plot
      plt.figure(0)
      plt.plot(wlarray,spectrum_W_m2)
      plt.xlabel("wavelength")
      plt.ylabel("W m^-2")

      #Spectra Plot
      plt.figure(1)
      plt.plot(wlarray,spectrum_photons_umol)
      plt.xlabel("wavelength")
      plt.ylabel("umol photons")


spectra_energy_array=np.vstack(spectra_energy[1:len(spectra_energy)])
spectra_photons_array=np.vstack(spectra_photons[1:len(spectra_photons)])

#First Epar Plot
plt.figure(2)
plt.plot(Eparlist,'k',marker=".")
plt.xlabel("Measurement Number")
plt.ylabel("Epar (Energy W m^-2)")

#Second Epar vs Depth Plot
plt.figure(3)
plt.scatter(Eparlist,specdepth)
plt.axhline(0,0,1)
plt.ylim(75,-5)
plt.xlabel("Epar (Energy W m^-2)")
plt.ylabel("Depth (m)")

#Create dataframe to export all energy spectra
df0=pd.DataFrame(spectra_energy_array-1)
df0.columns=wlarray
df0.insert(0,'timestamp', timestamp[0:len(spectra_energy_array)])
df0.to_csv(lis[4][3][0:8]+'Spectra_energy_'+date+'.csv')

#Create dataframe to export all photons spectra
df1=pd.DataFrame(spectra_photons_array-1)
df1.columns=wlarray
df1.insert(0,'timestamp', timestamp[0:len(spectra_photons_array)])
df1.to_csv(lis[4][3][0:8]+'Spectra_photons_'+date+'.csv')

#Create dataframe to export Epar data
result=[timestamp,specinttime,specdepth,Eparlist,Epar_photons_list,specsig]
data=np.array(result)
df = pd.DataFrame(data=data.T, index=[np.linspace(1,len(data.T) , len(data.T))], columns=["datetime","IntegrationTime","Depth","Epar Energy","Epar Photons","Signal %"])
df.to_csv(lis[4][3][0:8]+'DataFrame_'+date+'.csv')

#Remove empty data points
indices=np.argwhere(np.isnan(Eparlist))
timestamp=np.delete(timestamp,indices,0)
Eparlist=np.delete(Eparlist,indices,0)
Epar_photons_list=np.delete(Epar_photons_list,indices,0)


#Create dataframe for average of burst measurements
result2=[timestamp,Eparlist,Epar_photons_list]
data2=np.array(result2)
df2 = pd.DataFrame(data=data2.T, index=[np.linspace(1,len(data2.T) , len(data2.T))], columns=["datetime","Epar_Energy","Epar_Photons"])
df2['datetime'] = pd.to_datetime(df2['datetime'])
df2["Epar_Energy"] = pd.to_numeric(df2["Epar_Energy"], downcast="float")
df2["Epar_Photons"] = pd.to_numeric(df2["Epar_Photons"], downcast="float")
burstsumenergy=df2.resample('T', on='datetime').Epar_Energy.sum()
burstsumenergy = burstsumenergy.replace(0, np.nan).dropna()
burstsumphotons=df2.resample('T', on='datetime').Epar_Photons.sum()
burstsumphotons = burstsumphotons.replace(0, np.nan).dropna()
avg_energy=df2.resample('T', on='datetime').Epar_Energy.mean()
avg_energy= avg_energy.dropna()
avg_photons=df2.resample('T', on='datetime').Epar_Photons.mean()
avg_photons= avg_photons.dropna()

#Third Average Epar plot
plt.figure(4)
avg_energy.plot(y='Epar', use_index=True)
avg_photons.plot(y='Epar', use_index=True)
#burstsum.plot(y=Epar, use_index=True)
plt.show()
avg_energy.to_csv(lis[4][3][0:8]+'Average_energy_'+date+'.csv')
avg_photons.to_csv(lis[4][3][0:8]+'Average_photons_'+date+'.csv')

