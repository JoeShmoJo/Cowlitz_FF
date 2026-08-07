"""
Code to calculate an inflow dataset for the POR to feed to REsSim for the Cowlitz
Takes daily average inflows at Mossyrock from 2020 modified flow report, and shapes
them to an hourly shape by assuming the peak is a bit higher than the max 1-day.
Other inflow datasets did not seem as reliable (e.g CWMS)

The local flow between Mossyrock and Castle Rock (includes into Mayfield) is calculated
as the 1-day releases (from Modified Flow Report), shaped to hourly, then routed to 
Castle Rock. Then, the USGS records at Castle Rock (hourly if available, otherwise
dailies shaped to hourlies) are deducted to get a local flow.

For 2018-present, CWMS data is used instead of Modified Flows

Outputs a dss file with two pathnames--the Mossyrock inflow and Castle Rock locals

Messy script, but it works
Ryan Cahill
8/6/2026

"""

import pandas as pd
import numpy as np
import pickle
from scipy import interpolate
from numpy import log10
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, PercentFormatter, FixedLocator
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
from matplotlib.patches import Patch
import datetime
from datetime import timedelta
import dateutil
import os, sys
# Run-from-anywhere: resolve paths from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "Modules"))
from utilsDSS import HecDss
import requests
from pyextremes import get_extremes
from scipy.interpolate import interp1d
from scipy import stats
import HydrologicRouting
import HourlyHydrographShaping
from HydrologicRouting import SsarrReach
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

###############################################################################
#CONFIGURATION

DSS_FILE = os.path.join(REPO_ROOT, "CAS_Unreg_FF", "data", "obsData.dss")

OUT_DSS_FILE = "ResSimInflows.dss" 

CFSDAYS_TO_KAF = 86400/43560/1000.


path_dict_peak = {"Cowlitz_atRandle":"/COWLITZ RIVER AT RANDLE/14231000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/", #1993-2025
                  "Cowlitz_nrRandle":"/COWLITZ RIVER NEAR RANDLE, WA/14233400/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/", #1987-1993
                  "Cowlitz_atMossyrock":"/COWLITZ RIVER AT MOSSYROCK/14235000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/", #1906-1959
                  "Castle Rock":"/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/"
                  } 

path_dict_daily = {"Mossyrock_Inflow": "/MOS6A/MOS/FLOW-IN//1DAY/2020_Modified_Flow_Report/",
                  "Mossyrock_Outflow":"/MOS6H/MOS/FLOW-OUT//1DAY/2020_Modified_Flow_Report/",
                  "Cowlitz_atRandle":"/COWLITZ RIVER AT RANDLE/14231000/FLOW//1DAY/USGS/", #1993-2025
                  "Cowlitz_nrRandle":"/COWLITZ RIVER NEAR RANDLE, WA/14233400/FLOW//1DAY/USGS/", #1969-1993
                  "Cowlitz_atMossyrock":"/COWLITZ RIVER AT MOSSYROCK/14235000/FLOW//1DAY/USGS/", #1906-1959
                  "Castle Rock":"/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW//1DAY/USGS/"
                  } 



path_dict_hourly = {"Castle Rock": "/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW//1HOUR/USGS/",
              "Mossyrock_Elev": "//MOS/ELEV//1HOUR/CWMS/",
              "Mayfield_Outflow":"//MAYFIELD OUTFLOW/FLOW//1HOUR/USGS 14238000/",
              "Cowlitz_atRandle":"/COWLITZ RIVER AT RANDLE/14231000/FLOW//1HOUR/USGS/", #1993-2025
              "Cowlitz_nrRandle":"/COWLITZ RIVER NEAR RANDLE, WA/14233400/FLOW//1HOUR/USGS/", #1987-1993
              "Mossyrock_Outflow": "//MOS/FLOW-OUT//1HOUR/CWMS/", #for 2018-present
              "Mossyrock_Inflow":"//MOS/FLOW-IN//1HOUR/CWMS/" #for 2018-present
              } 


# Table 3 calibration, lower Cowlitz River model. Only the Cowlitz mainstem
# chain Mayfield -> Castle Rock is used; the Toutle River's own reach
# ("Tower to Cowlitz+Toutle") is excluded since only the holdout signal is
# being routed, not the Toutle's inflow. Each tuple is (kts, n, numSubreaches/phases).
COWLITZ_REACH_PARAMS = [
    (5, 0.1, 5),  # Mayfield_OUT -> Cowlitz R above Toutle R
    (1, 0.2, 1),  # Cowlitz R above Toutle R -> Cowlitz+Toutle
    (1, 0.2, 5),  # Cowlitz+Toutle -> Castle Rock
]

ROUTING_TIMESTEP_HRS = 1

################################################################################
# CLASS DEFINITIONS
    


################################################################################
# FUNCTION DEFINITIONS

def getWY(datetimeObj):
    year = datetimeObj.year
    if datetimeObj.month>=10:
        return(year+1)
    else:
        return(year)

def getDate0000(datetimeObj):
    ''' returns a datetime object with no hour, minute, or second'''
    return datetimeObj.replace(hour=0, minute=0, second=0)

def getDate2000(datetimeObj):
    ''' returns a datetime object set to water year 2000 (leap year) (keeps month and day, but resets year)'''
    if datetimeObj.month>=10:
        return datetimeObj.replace(year=1999)
    else:
        return datetimeObj.replace(year=2000)

def getDate2001(datetimeObj):
    ''' returns a datetime object set to water year 2001 (no leap year) (keeps month and day, but resets year)'''
    if datetimeObj.month>=10:
        return datetimeObj.replace(year=2000)
    elif datetimeObj.month == 2 and datetimeObj.day == 29:
        return np.nan
    else:
        return datetimeObj.replace(year=2001)

def getMonthAbbr(monthNum):
    '''
    Given a month number (e.g. 4), return month abbreviation (e.g. "Apr")
    '''
    return datetime.datetime(2001,monthNum,1).strftime("%b")

def readAllData():
    '''
    Reads all time series data from dss and store it to a dictionary
    1st keys is timestep, 2nd key is location, value is a dataframe
    with index as datetime, and column names are the locations
    '''
    dataDict = {"Daily":{},"Peak":{}, "Hourly":{}}
    
    dssFile = HecDss.open(DSS_FILE)
    #Do daily
    for name, pathname in path_dict_daily.items():
        print("  Reading: %s" %pathname)        
        df = dssFile.readDF(pathname).dropna()
        df = df.rename(columns={"value":name})
        dataDict["Daily"][name] = df
    #Do hourly
    for name, pathname in path_dict_hourly.items():
        print("  Reading: %s" %pathname)        
        df = dssFile.readDF(pathname).dropna()
        df = df.rename(columns={"value":name})
        dataDict["Hourly"][name] = df
    #Do daily
    for name, pathname in path_dict_peak.items():
        print("  Reading: %s" %pathname)        
        df = dssFile.readDF(pathname).dropna()
        df = df.rename(columns={"value":name})
        dataDict["Peak"][name] = df
    dssFile.close()
    return dataDict
    
def getnDayMaxAroundPeak(peakDateList, locName="Castle Rock", nDays=3):
    #Gets n-day max around a list of peak dates
    #Will go + and minus 4 days
    rows = []
    for peakDate in peakDateList:
        dfDaily = dataDict["Daily"][locName]
        startDate = peakDate - timedelta(days=4)
        endDate = peakDate + timedelta(days=4)
        dfWindow = dfDaily.loc[((dfDaily.index >=startDate) & (dfDaily.index <= endDate))]
        if len(dfWindow) == 0:
            continue
        dfRolling = dfWindow.rolling(nDays, center=True).mean()
        rows.append([dfRolling.idxmax().iloc[0], dfRolling.max().iloc[0]])
    
    dfCoincident = pd.DataFrame(rows, columns=[f"Date_{locName}", locName])
    dfCoincident["WY"] = dfCoincident[f"Date_{locName}"].apply(getWY)
    #dfCoincident = pd.concat(dfs, ignore_index=True)
    return dfCoincident

def build_cowlitz_reaches():
    """Fresh SsarrReach objects for the three-reach Mayfield -> Castle
    Rock chain (COWLITZ_REACH_PARAMS). Built fresh per call so a route_reach()
    call never reuses another call's internal state (subreachOutflows)."""
    reaches = []
    for kts, n, num_subreaches in COWLITZ_REACH_PARAMS:
        reach = SsarrReach(timestepHrs=ROUTING_TIMESTEP_HRS)
        reach.buildWithKTS(numSubreaches=num_subreaches, n=n, kts=kts)
        reaches.append(reach)
    return reaches


def route_reach(flow_series):
    """
    Route an hourly flow series through the three chained SSARR reaches
    spanning Mayfield -> Castle Rock (Table 3 calibration). flow_series
    should already be non-negative (see cap_negative_flows) -- passed
    with allowNegatives=True since that capping is done explicitly
    upstream rather than left to SsarrReach's own floor-to-1-cfs
    behavior. Returns a pandas Series on the same index as flow_series.
    """
    values = flow_series.values.tolist()
    for reach in build_cowlitz_reaches():
        values = reach.routeHydrograph(values, allowNegatives=True)
    return pd.Series(values, index=flow_series.index)


#####  Main  ###################################################################
# MAIN CODE
if __name__ == "main" or __name__ == "__main__":
    #Get the data
    dataDict = readAllData()

#Generate relationship of highest daily average flow vs peaks
#Mash 3 different datasets together to represent Mossyrock inflow

#For the locations where we have hourly data, grab the annual instantaneous peak and 1-day peak
dict_dfPeaks = {}
for locName in ["Cowlitz_atRandle", "Cowlitz_nrRandle"]:
    df = dataDict["Hourly"][locName]
    dfDaily = df.resample("1d").mean()
    df["WY"] = df.index.map(getWY)
    dfPeaks = df.groupby("WY").max()
    dfPeaks = dfPeaks.rename(columns={locName:"Peak"})
    dfPeaks["PeakDate"] = df.groupby("WY").idxmax()
    dfPeaks = dfPeaks.reset_index()
    #Add the daily average max. Not 24-hour rolling, but midnight to midnight
    dfDaily["WY"] = dfDaily.index.map(getWY)
    dfDailyMax = dfDaily.groupby("WY").max().rename(columns={locName:"DailyMax"})
    dfDailyMax = dfDailyMax.reset_index()
    dfPeaks = pd.merge(dfPeaks, dfDailyMax, on="WY", how="inner")
    dfPeaks = dfPeaks.dropna()
    dfPeaks["Location"] = locName
    dict_dfPeaks[locName] = dfPeaks

#For the locations where we have daily avg and annual peak, mush them together
dict_dfAnnualPeaks = {}
for locName in ["Cowlitz_atMossyrock", "Cowlitz_atRandle", "Cowlitz_nrRandle"]:
    dfPeaks = dataDict["Peak"][locName].copy()
    dfPeaks.index.name = "PeakDate"
    dfPeaks["WY"] = dfPeaks.index.map(getWY)
    dfPeaks = dfPeaks.rename(columns={locName:"Peak"}).reset_index()
    #Now process the daily average to get 1-day max for each wy
    dfDaily = dataDict["Daily"][locName].copy()
    dfDaily["WY"] = dfDaily.index.map(getWY)
    dfDailyMax = dfDaily.groupby("WY").max()
    dfDailyMax = dfDailyMax.rename(columns={locName:"DailyMax"})
    dfPeaks = pd.merge(dfPeaks, dfDailyMax, on="WY", how="inner")
    dfPeaks = dfPeaks.dropna()
    dfPeaks["Location"] = locName
    dict_dfAnnualPeaks[locName] = dfPeaks

#Merge it all together
dfAll = pd.concat([dict_dfPeaks["Cowlitz_atRandle"],
                  dict_dfPeaks["Cowlitz_nrRandle"],
                  dict_dfAnnualPeaks["Cowlitz_nrRandle"],
                  dict_dfAnnualPeaks["Cowlitz_atRandle"],
                  dict_dfAnnualPeaks["Cowlitz_atMossyrock"]]
                  )

#Get rid of years with bad data
dfAll = dfAll.loc[~(dfAll["WY"]==2002)] #missing hourly peak
dfAll = dfAll.loc[~(dfAll["WY"].isin([1914,1915,1929]))] #old
dfAll = dfAll.loc[~(dfAll["WY"].isin([1979]))] #daily same as inst
dfAll = dfAll.loc[~((dfAll["WY"]==1993)&(dfAll["Location"]=="Cowlitz_atRandle"))] #overlap with nr Randle
#Some duplicates because using different methods
dfAll = dfAll.drop_duplicates(subset="WY")
#The main number--what is the ratio we need to factor up a daily max value to get a peak estimate?
avgRatio = ((dfAll["Peak"]-dfAll["DailyMax"])/dfAll["DailyMax"]).mean()

#Plot it up
plt.close('all')
fig, ax = plt.subplots()
for locName in ["Cowlitz_atMossyrock", "Cowlitz_atRandle", "Cowlitz_nrRandle"]:
    dfPlot = dfAll.loc[dfAll["Location"]==locName].copy()
    dfPlot["y"] = (dfPlot["Peak"]-dfPlot["DailyMax"])/dfPlot["DailyMax"]
    ax.scatter(dfPlot["DailyMax"],dfPlot["y"], label=locName)
    for i, rowPlot in dfPlot.iterrows():
        ax.text(rowPlot["DailyMax"], rowPlot["y"], rowPlot["WY"], ha="left", va="top")
ax.axhline(avgRatio, ls="dashed", color="gray")
ax.text(50000, avgRatio, f"Average Increase = {avgRatio*100:0.0f}%", ha="left", va="top")
ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
ax.set_ylabel("Difference between Instantaneous Peak and Daily Average (%)")
ax.set_xlabel("Daily Average maximum flow (cfs)")
ax.legend()



###########################################################################
#Now that we have the adjustment factor, we need to create an hourly inflow dataset from daily
#First, create a dataframe with faked in peaks.
#This is the 1-day flow scaled up by the adjustment ratio
#Then, apply the smoothing with the dailies and fake peaks supplied
tsName = "Mossyrock_Inflow"
dfDaily = dataDict["Daily"][tsName]
peakIndices, _ = find_peaks(dfDaily[tsName].values, distance=30)
dfFakePeaks = dfDaily.iloc[peakIndices]*(1 + avgRatio)

#Standardize the column names
dfDaily = dfDaily[[tsName]].rename(columns={tsName:"value"})
dfFakePeaks = dfFakePeaks[[tsName]].rename(columns={tsName:"value"}).copy()

dfFakePeaks.index = dfFakePeaks.index

#dfDaily = dfDaily.loc[(dfDaily.index > datetime.datetime(1996,2,1)) & (dfDaily.index < datetime.datetime(1996,2,15))]
#dfFakePeaks = dfFakePeaks.loc[(dfFakePeaks.index > datetime.datetime(1996,2,1)) & (dfFakePeaks.index < datetime.datetime(1996,2,15))]

dfHourlyInflow = HourlyHydrographShaping.createHourlyUsingSpline(dfDaily, dfFakePeaks)

#Tack on the data after 2018
df2018 = dataDict["Hourly"][tsName]
df2018 = df2018.loc[df2018.index >= datetime.datetime(2018,9,30,13)]
df2018 = df2018.rename(columns={tsName:"value"})
dfHourlyInflow = pd.concat([dfHourlyInflow, df2018])

#Check one event
fig, ax = plt.subplots()
ax.step(dfDaily.index, dfDaily["value"])
ax.scatter(dfFakePeaks.index, dfFakePeaks["value"])
ax.plot(dfHourlyInflow.index, dfHourlyInflow["value"])
ax.set_xlim([datetime.datetime(1996,2,1),datetime.datetime(1996,2,15)])


##################################################################
#Next step is to calculate the local between Mossyrock and Castle Rock
#For Mossyrock outflow, don't do peak data, just use a basic spline interpolator
#Since it bounces so much and is heavily regulated

tsName = "Mossyrock_Outflow"

df = dataDict["Daily"][tsName]
x_daily = np.arange(0, len(df), 1)

times_hourly = pd.date_range(df.index[0], df.index[-1], freq="1h")
x_hourly = np.arange(0, len(df)-.99, 1/24)

#Shift the date to noon instead of midnight to do the interpolation
yInterp = PchipInterpolator(x_daily-.5, df[tsName])(x_hourly)
dfHourlyOutflow = pd.DataFrame(index=times_hourly, data=yInterp, columns=["Release"])

#Tack on the data after 2018
df2018 = dataDict["Hourly"][tsName]
df2018 = df2018.loc[df2018.index >= datetime.datetime(2018,10,1,1)]
df2018 = df2018.rename(columns={tsName:"Release"})
dfHourlyOutflow = pd.concat([dfHourlyOutflow, df2018])

#Route the outflow down to Castle Rock
dfHourlyOutflow["Routed"] = route_reach(dfHourlyOutflow["Release"])

#Tack on the observed hourly castle rock
#Shape it
dfDaily = dataDict["Daily"]["Castle Rock"].rename(columns={"Castle Rock":"value"})
dfPeaks = dataDict["Peak"]["Castle Rock"].rename(columns={"Castle Rock":"value"})
print("taking a while to shape a hydrograph...")
dfCAS = HourlyHydrographShaping.createHourlyUsingSpline(dfDaily, dfPeaks)
dfCAS = dfCAS.rename(columns={"value":"Castle Rock"})
#Overwrite the shaped values with real values when we have them
#Fill in some small gaps
seriesCASObs = dataDict["Hourly"]["Castle Rock"]["Castle Rock"].interpolate(method='linear',limit=6)
dfCAS.loc[seriesCASObs.index, "Castle Rock"] = seriesCASObs
dfHourlyOutflow["Castle Rock"] = dfCAS["Castle Rock"]
dfHourlyOutflow["Local"] = dfHourlyOutflow["Castle Rock"] - dfHourlyOutflow["Routed"]
dfHourlyOutflow["Local"] = dfHourlyOutflow["Local"].rolling("3h", center=True).mean()
dfHourlyOutflow["Local"] = dfHourlyOutflow["Local"].clip(0)

#Write out the inflows and locals to DSS


dssFile = HecDss.open(OUT_DSS_FILE)
dssFile.writeDF(dfHourlyOutflow["Local"], "//CASTLE ROCK/FLOW-LOCAL//1HOUR/FOR_RESSIM/", "cfs", "INST-VAL")
dssFile.writeDF(dfHourlyInflow, "//MOSSYROCK/FLOW-IN//1HOUR/FOR_RESSIM/", "cfs", "INST-VAL")
dssFile.close()

fix, ax = plt.subplots()
ax.step(df.index, df[tsName])
ax.plot(times_hourly, yInterp)

ax.set_xlim([datetime.datetime(1996,2,1),datetime.datetime(1996,2,15)])
            