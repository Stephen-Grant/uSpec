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

#Open the Log file and import all data in 
with open(r"C:\Users\stepheg\Downloads\uSpec_240415-16\LOG_0082.txt") as f: #Enter path for LOG file here <----
   lis = [line.split() for line in f]        

   wl=(lis[37][2])                                                                  #Extract Wavelengths and convert to array
   wlarray=(np.array(wl.split(",")))                                        
   wlarray[0]=wlarray[0].translate({ord('{'): None})
   wlarray[268]=wlarray[268].translate({ord('}'): None})
   wlarray=wlarray.astype(np.float64)

   #Empty variables for holding data, pre allocated for speed
   spectraEpar=[];specdate=[];spectime=[];specsig=[];specinttime=[];specdepth=[];Eparlist=[];
   Eparredlist=[];Epargreenlist=[];Eparbluelist=[];Eparbandstotallist=[];timestamp=[];

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
      Epar=(integrate.trapezoid(spectrum[15:146],wlarray[15:146]))*1E-6; #Convert from uW to W
      Epar=Epar*10000; #Convert from cm^-2 to m^-2
      Eparlist.append(Epar)
      spectraEpar.append(spectrum[15:146])
      
spectraEpararray=np.vstack(spectraEpar[1:len(spectraEpar)])

#First Epar Plot
plt.figure(0)
plt.plot(Eparlist,'k',marker=".")
plt.xlabel("Measurement Number")
plt.ylabel("Epar")

#Second Epar vs Depth Plot
plt.figure(1)
plt.scatter(Eparlist,specdepth)
plt.axhline(0,0,1)
plt.ylim(75,-5)
plt.xlabel("Epar")
plt.ylabel("Depth")

#Create dataframe to export all spectra
df0=pd.DataFrame(spectraEpararray)
df0.columns=wlarray[15:146]
df0.insert(0,'timestamp', timestamp[1:len(spectraEpar)])
df0.to_csv(lis[4][3][0:8]+'Spectra_'+date+'.csv')

#Create dataframe to export Epar data
result=[timestamp,specinttime,specdepth,Eparlist,specsig]
data=np.array(result)
df = pd.DataFrame(data=data.T, index=[np.linspace(1,len(data.T) , len(data.T))], columns=["datetime","IntegrationTime","Depth","Epar","Signal %"])
df.to_csv(lis[4][3][0:8]+'DataFrame_'+date+'.csv')

#Remove empty data points
indices=np.argwhere(np.isnan(Eparlist))
timestamp=np.delete(timestamp,indices,0)
Eparlist=np.delete(Eparlist,indices,0)

#Create dataframe for average of burst measurements
result2=[timestamp,Eparlist]
data2=np.array(result2)
df2 = pd.DataFrame(data=data2.T, index=[np.linspace(1,len(data2.T) , len(data2.T))], columns=["datetime","Epar"])
df2['datetime'] = pd.to_datetime(df2['datetime'])
df2["Epar"] = pd.to_numeric(df2["Epar"], downcast="float")
burstsum=df2.resample('T', on='datetime').Epar.sum()
burstsum = burstsum.replace(0, np.nan).dropna()
avg=df2.resample('T', on='datetime').Epar.mean()
avg = avg.dropna()

#Third Average Epar plot
plt.figure(2)
avg.plot(y='Epar', use_index=True)
#burstsum.plot(y=Epar, use_index=True)
plt.show()
avg.to_csv(lis[4][3][0:8]+'Average_'+date+'.csv')
