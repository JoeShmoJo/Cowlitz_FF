"""
Contains code to do some hydrology work in the Cowlitz Basin
The main point is to come up with annual maximum unregulated flow estimates at Castle Rock

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
from utilsDSS import HecDss
import requests
from pyextremes import get_extremes
from scipy.interpolate import interp1d
from scipy import stats
import HydrologicRouting

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

###############################################################################
#CONFIGURATION

CONFIG_FILE = "CowlitzConfig.xlsx"
DSS_FILE = "ObsData.dss" 

READ_FROM_DSS = True
DATA_PICKLE_FILE = "AllData.pickle" 

CFSDAYS_TO_KAF = 86400/43560/1000.

#pathnames in the dss file
#TODO need to load up annual observed peak only at Castle Rock
path_dict_peak = {"Castle Rock": ""}

#TODO need to load up daily flows and elevations
path_dict_daily = {"Castle Rock": "",
                  "Mossyrock_Elev":""}

path_dict_hourly = {"Castle Rock": "/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW//1HOUR/USGS/",
              "Mossyrock_Elev": "//MOS/ELEV//1HOUR/CWMS/",
              "Mayfield_Outflow":"//MAYFIELD OUTFLOW/FLOW//1HOUR/USGS 14238000/"}

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
    #Do hourly
    for name, pathname in path_dict_hourly.items():
        print("  Reading: %s" %pathname)        
        df = dssFile.readDF(pathname).dropna()
        df = df.rename(columns={"value":name})
        dataDict["Hourly"][name] = df
    '''
    TODO load up daily/peaks when available
    '''

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

def getStor_ChangeAroundPeak(peakDateList, window_days=1):
    #Given a list of dates of peak flow at Castle Rock,
    #Return a dataframe with the corresponding change in storage (converted to flow)
    #Assumes we only have end of day readings at Mossyrock, since we'll be using that for the regression
    #window_days is how far back to look before the date of peak. IF 1, then only the data on the day of peak is used
    #smoothing_days is if we don't want a 1-day flow, but maybe a 2-day or 3-day to get a better regression
    locName = "Castle Rock_FlowChange"
    rows = []
    for peakDate in peakDateList:
        dfDaily = dataDict["Daily"][locName]
        #Erase any hours, since daily data doesn't have it
        endDate = datetime.datetime(peakDate.year, peakDate.month, peakDate.day)
        startDate = endDate - timedelta(days=window_days)
        dfWindow = dfDaily.loc[((dfDaily.index >startDate) & (dfDaily.index <= endDate))]
        if len(dfWindow) == 0:
            continue
        maxVal = dfWindow["FlowChange"].mean()
        rows.append([peakDate, maxVal])
    
    dfCoincident = pd.DataFrame(rows, columns=[f"Date_{locName}", locName])
    dfCoincident["WY"] = dfCoincident[f"Date_{locName}"].apply(getWY)
    return dfCoincident

def getMaxAroundPeak(dfCowlitz, dfRegUnregEvents):
    #Get max increase from regulated to unregulated within a couple days of the time of peak at Castle Rock
    rows = []
    for i, rowDict in dfRegUnregEvents.iterrows():
        peakDate = rowDict["PeakDate"]
        startDate = rowDict["StartDate"]
        endDate = rowDict["EndDate"]
        dfWindow = dfCowlitz.loc[((dfCowlitz.index >=startDate) & (dfCowlitz.index <= endDate))]
        if len(dfWindow) == 0:
            continue
        #Get max increase
        maxUnreg = dfWindow["Castle Rock_Unreg"].max()
        maxObs = dfWindow["Castle Rock_Obs"].max()
        
        peakIncrease = maxUnreg - maxObs
        rows.append([peakDate, peakIncrease, maxUnreg, maxObs])
    
    dfCoincident = pd.DataFrame(rows, columns=["Date_RegToUnreg", "RegToUnreg", "UnregPeak", "ObsPeak"])
    dfCoincident["WY"] = dfCoincident["Date_RegToUnreg"].apply(getWY)
    #dfCoincident = pd.concat(dfs, ignore_index=True)
    return dfCoincident

def plotUnregvsObsForEveryEvent():
    #Plot each event
    plt.close("all")
    for i, dfGroup in dfRegUnregEvents.groupby("PlotGroup"):
        rowIdx = 0
        colIdx = 0
        fig, axes = plt.subplots(2,3, figsize=[11,7])
        for j, rowDict in dfGroup.iterrows():
            ax = axes[rowIdx, colIdx]
            startDate = rowDict["StartDate"]
            endDate = rowDict["EndDate"]
            dfWindow = dfCowlitz.loc[((dfCowlitz.index >=startDate) & (dfCowlitz.index <= endDate))]
            ax.plot(dfWindow.index, dfWindow["Castle Rock_Unreg"], label="Unreg")
            ax.plot(dfWindow.index, dfWindow["Castle Rock_Obs"], label="Obs")
            ax.set_ylabel("Flow (cfs)")
            ax.xaxis.set_major_formatter(DateFormatter("%b-%d"))
            ax.xaxis.set_major_locator(mdates.DayLocator())
            ax.text(0.5,0.05, datetime.datetime.strftime(rowDict["PeakDate"], "%b-%Y"), ha='center', transform = ax.transAxes, bbox=dict(facecolor='whitesmoke', edgecolor='dimgray', boxstyle='round'))
            ax.legend()    
            
            colIdx += 1 
            if colIdx > 2: 
                colIdx = 0
                rowIdx += 1
        plt.tight_layout()

def plotStorChangeVsPeakIncreaseRegression(dfRegToUnregCoincident):
    """
    Get the regression that relates change in storage at Mossyrock to the change in peak flow
    Returns the slope and intercept of the regression line

    """
    fig, ax = plt.subplots()
    durationStr  = "2" #2-day
    xColName = f"Mossyrock_FlowChangeMax{durationStr}"
    df = dfRegToUnregCoincident.copy()
    ax.scatter(df[xColName], df["RegToUnreg"], label=f"{durationStr}-day Max")
    #ax.scatter(df[xColName], df["RegToUnreg"], label=f"{durationStr}-day Max, r={corr:.2f}")
    for i, rowDict in df.iterrows():
        dateThing = datetime.datetime(rowDict["WY"], rowDict["Date_RegToUnreg"].month, rowDict["Date_RegToUnreg"].day)
        if dateThing.month >9:
            dateThing = dateThing - timedelta(days=365)
        dateStr = datetime.datetime.strftime(dateThing, "%b-%Y")
        ax.text(rowDict["Mossyrock_FlowChangeMax2"],rowDict["RegToUnreg"], dateStr, ha="center", va="bottom")

    #Plot Dec 2025 as a red dot
    ax.scatter(df[xColName].values[-1], df["RegToUnreg"].values[-1], color="tab:red", label=None)
    #Drop the December 2025 event (last one)
    df = df.iloc[:-1].copy()
    corr = df.corr()[xColName]["RegToUnreg"]
    
    #Get the regression equation
    slope, intercept, r_value, p_value, stdErrRegression = stats.linregress(df[xColName],df["RegToUnreg"])
    xLine = np.array(ax.get_xlim())
    yLine = xLine*slope + intercept
    ax.plot(xLine, yLine, color="tab:blue")
    
    ax.text(0.95,0.6,f"y={slope:.2f}x + {intercept:.0f}\nr={r_value:.2f}", ha='right', transform = ax.transAxes, bbox=dict(facecolor='whitesmoke', edgecolor='dimgray', boxstyle='round'))
    
    #ax.legend()
    ax.set_xlabel("Mossyrock 2-day max storage change (converted to flow) near date of Castle Rock peak (cfs)")
    ax.set_ylabel("Unregulated peak minus regulated peak (cfs)")
    ax.set_title("Castle Rock Peak Flow Increase Regression")
    return slope, intercept

def applyRegressionToGetUnregPeaks(slope, intercept):
    '''
    We have a regression on change in storage in Mossyrock Lake to the increase in peak flow from obs to unreg
    Let's apply this for early so that we have a complete unregulated dataset
    Use a 2-day change in storage

    '''
    dfCastleRock = dataDict["Peak"]["Castle Rock"].copy()
    dfCastleRock["WY"] = dfCastleRock.index.map(getWY)
    dfMossyrock2Day = getnDayMaxAroundPeak(dfCastleRock.index, "Mossyrock_FlowChange", 2)
    dfMerge = pd.merge(dfCastleRock, dfMossyrock2Day, how="outer", on="WY")
    dfMerge["Unreg_Castle Rock"] = dfMerge["Flow"] + slope*dfMerge["Mossyrock_FlowChange"] + intercept
    
    #override for specific water years
    #dfMerge.loc[dfMerge["WY"].isin([1996,1991,1976]), "Unreg_Castle Rock"] = dfMerge["Flow"]
    dfMerge = dfMerge.set_index("WY")
    
    #We have good data after 1990, merge in the older estimates with the new estimates
    dfRecent = dfCowlitz.copy()
    dfRecent["WY"] = dfRecent.index.map(getWY) 
    dfRecent = dfRecent.loc[dfRecent["WY"]>1990].copy()
    dfRecent = dfRecent.groupby("WY")["Castle Rock_Unreg"].max()
    
    dfMerge.loc[dfRecent.index, "Unreg_Castle Rock"] = dfRecent
    dfOut = dfMerge[["Unreg_Castle Rock", "Flow"]]
    dfOut = dfOut.rename(columns={"Flow":"Obs_Castle Rock"})
    
    fig, ax = plt.subplots()
    df = dfOut.loc[dfOut.index>=1970]
    ax.scatter(df["Obs_Castle Rock"], df["Unreg_Castle Rock"]-df["Obs_Castle Rock"])
    ax.set_xlabel("Observed Annual Peak Flow at Castle Rock (cfs)")
    ax.set_ylabel("Increase from Observed to Unregulated (cfs)")
    
    fig, ax = plt.subplots()
    ax.scatter(dfOut.index, dfOut["Unreg_Castle Rock"], label="Unregulated")
    ax.scatter(dfOut.index, dfOut["Obs_Castle Rock"], label="Observed")
    ax.set_xlabel("Water Year")
    ax.set_ylabel("Annual Peak Flow (cfs)")
    ax.legend()
    plt.tight_layout()
    
    return dfOut

def getNDayUnregAnnualMax(dfUnregDaily, durations=[3,5]):
    '''
    Use the daily average unregulated flow estimate for the POR to come up with
    n-day unregulated flow estimates
    Ignores the effect of routing, should only be used for durations longer than 3 days
    '''

    outDict = {}
    for durDays in durations:
        duration = f"{durDays}-Day"
        #Look at unreg
        dfUnregMax = dfUnregDaily.rolling(f"{durDays}D", center=True).mean()
        dfUnregMax["WY"] = dfUnregMax.index.map(getWY)
        maxFlowSeries = dfUnregMax.groupby("WY")["Unregulated"].max()
        dateMaxSeries = dfUnregMax.groupby("WY")["Unregulated"].idxmax()
        dfAMSUnreg = pd.DataFrame()
        dfAMSUnreg["WY"] = maxFlowSeries.index
        dfAMSUnreg["Unreg_Orig"] = maxFlowSeries.values
        dfAMSUnreg["Unreg_Date"] = dateMaxSeries.values
        
        #Repeat for observed
        dfObsMax = dataDict["Daily"]["Castle Rock"].rolling(f"{durDays}D", center=True).mean()
        dfObsMax["WY"] = dfObsMax.index.map(getWY)
        maxFlowSeries = dfObsMax.groupby("WY")["Flow"].max()
        dateMaxSeries = dfObsMax.groupby("WY")["Flow"].idxmax()
        dfAMSObs = pd.DataFrame()
        dfAMSObs["WY"] = maxFlowSeries.index
        dfAMSObs["Obs"] = maxFlowSeries.values
        dfAMSObs["Obs_Date"] = dateMaxSeries.values
        
        #Merge them together
        dfAMS = pd.merge(dfAMSUnreg, dfAMSObs, how="outer", on="WY")
        #Create a "final" unreg dataset that merges the unreg and observed
        #Start it as the original calculated unreg
        dfAMS["Unreg"] = dfAMS["Unreg_Orig"]
        #Mossyrock Lake had no data in WY 1994
        dfAMS.loc[dfAMS["WY"].isin([1994]), "Unreg"] = dfAMS["Obs"]
        #Mossyrock Lake was much smaller before it was enlarged. WY 1984 is the first year with the big pool
        dfAMS.loc[dfAMS["WY"]<1984, "Unreg"] = dfAMS["Obs"]
        #Never let the unreg estimate be lower than the observed
        dfAMS.loc[dfAMS["Unreg"]<dfAMS["Obs"], "Unreg"] = dfAMS["Obs"]
        #Drop WY 1963, the gage started in January but the peak was probably in Fall 1962
        dfAMS = dfAMS.loc[dfAMS["WY"]>1963].copy()
        #
        dfAMS["Diff"] = dfAMS["Unreg"] - dfAMS["Obs"]
        
        outDict[duration] = dfAMS
        
        #Make a plot
        fig, ax = plt.subplots()
        ax.scatter(dfAMS["WY"], dfAMS["Unreg"], color="tab:blue", label="Unregulated")
        ax.scatter(dfAMS["WY"], dfAMS["Obs"], color="tab:orange", label="Observed")
        ax.legend()
        ax.set_xlabel("Water Year")
        ax.set_ylabel(f"Annual Max {duration} Flow (cfs)")
        plt.tight_layout()
    return outDict


#####  Main  ###################################################################
# MAIN CODE
if __name__ == "main" or __name__ == "__main__":
    dfElevStor = pd.read_excel(CONFIG_FILE,sheet_name="ElevStor",header=3)
    dfElevStor = dfElevStor.sort_values("Elev")
    #Get the data
    if READ_FROM_DSS:
        #Read all data from dss and store to a pickle file
        dataDict = readAllData()
        with open(DATA_PICKLE_FILE, 'wb') as f:
            pickle.dump(dataDict, f)
            
    with open(DATA_PICKLE_FILE, 'rb') as f:
        dataDict = pickle.load(f) 
        
    #Get daily change in storage of Mossyrock Lake (in cfs)
    # Create the interpolation function based on reference data
    interpolator = interp1d(dfElevStor["Elev"], dfElevStor["Stor"], kind='linear', fill_value="extrapolate")
    '''
    TODO add this daily data in when appropriate
    
    dfMossyrock = dataDict["Daily"]["Mossyrock_Elev"]
    dfMossyrock["Stor"] = interpolator(dfMossyrock["Elev"])
    dfMossyrock["FlowChange"] = dfMossyrock["Stor"].diff()*43560/86400
    dataDict["Daily"]["Mossyrock_FlowChange"] = dfMossyrock[["FlowChange"]].dropna()
    
    #Estimate daily unreg flows at Castle Rock (assuming no routing)
    dfUnregDaily = dataDict["Daily"]["Castle Rock"].rename(columns={"Flow":"Regulated"})
    dfUnregDaily["FlowChange"] = dfMossyrock["FlowChange"]
    dfUnregDaily = dfUnregDaily.dropna()
    dfUnregDaily["Unregulated"] = dfUnregDaily["Regulated"] + dfUnregDaily["FlowChange"]
    dfUnregDaily = dfUnregDaily.drop(columns="FlowChange")
    '''
    
    #Estimate unreg timeseries at Cowlitz downstream of Mayfield using short-duration data
    dfMossyrockHourly = dataDict["Hourly"]["Mossyrock_Elev"]
    dfMossyrockHourly["Stor"] = interpolator(dfMossyrockHourly["Mossyrock_Elev"])
    dfMossyrockHourly["FlowChange"] = dfMossyrockHourly["Stor"].diff()*43560/(60*60)
    dfCowlitz = dataDict["Hourly"]["Mayfield_Outflow"]
    dfCowlitz["FlowChange"] = dfMossyrockHourly["FlowChange"]
    dfCowlitz["Unreg"] = (dfCowlitz["Mayfield_Outflow"] + dfCowlitz["FlowChange"]).clip(0)
    dfCowlitz = dfCowlitz.drop(columns=["FlowChange"])
    dfCowlitz = dfCowlitz.dropna()
    
    #Route both observed and unreg down to Castle Rock
    #TODO Chop to a short period for demonstration
    dfCowlitz = dfCowlitz.loc[(dfCowlitz.index > datetime.datetime(2019,1,5))].copy()
    reachObj = HydrologicRouting.SsarrReach(timestepHrs=1)
    numSubreaches = 8
    n = 0.2
    kts = 4
    reachObj.buildWithKTS(numSubreaches, n, kts)
    print("Routing flows")
    dfCowlitz["Unreg_Routed"] = reachObj.routeHydrograph(dfCowlitz["Unreg"].values)
    dfCowlitz["Reg_Routed"] = reachObj.routeHydrograph(dfCowlitz["Mayfield_Outflow"].values)
    
    #Get the difference to be added to the Castle Rock record
    dfCowlitz["Castle Rock_Increase"] = dfCowlitz["Unreg_Routed"] - dfCowlitz["Reg_Routed"]
    dfCowlitz["Castle Rock_Obs"] = dataDict["Hourly"]["Castle Rock"]
    dfCowlitz["Castle Rock_Unreg"] = dfCowlitz["Castle Rock_Obs"] + dfCowlitz["Castle Rock_Increase"]
    
    #Do the analysis/plots
    plt.close('all')
    dfRegUnregEvents = pd.read_excel(CONFIG_FILE, header=3, sheet_name="RegUnregEvents")
    plotUnregvsObsForEveryEvent()
    
    #Now shift to reg vs unreg regression relationship
    #When looking at reg vs unreg, want more than just the annual max to enrich the dataset
    '''
    dfRegUnregEvents = pd.read_excel(CONFIG_FILE, header=3, sheet_name="RegUnregEvents")
    dfRegToUnregCoincident = getMaxAroundPeak(dfCowlitz, dfRegUnregEvents)
    dfMossyrockFlowChange = getnDayMaxAroundPeak(dfRegToUnregCoincident["Date_RegToUnreg"], "Mossyrock_FlowChange", 1)
    dfRegToUnregCoincident["Mossyrock_FlowChangeMax1"] = dfMossyrockFlowChange["Mossyrock_FlowChange"].values
    dfMossyrockFlowChange = getnDayMaxAroundPeak(dfRegToUnregCoincident["Date_RegToUnreg"], "Mossyrock_FlowChange", 2)
    dfRegToUnregCoincident["Mossyrock_FlowChangeMax2"] = dfMossyrockFlowChange["Mossyrock_FlowChange"].values
    dfMossyrockFlowChange = getnDayMaxAroundPeak(dfRegToUnregCoincident["Date_RegToUnreg"], "Mossyrock_FlowChange", 3)
    dfRegToUnregCoincident["Mossyrock_FlowChangeMax3"] = dfMossyrockFlowChange["Mossyrock_FlowChange"].values
    dfRegToUnregCoincident["Mossyrock_FlowChange1"] = getStor_ChangeAroundPeak(dfRegToUnregCoincident["Date_RegToUnreg"], window_days=1)["Mossyrock_FlowChange"].values
    dfRegToUnregCoincident["Mossyrock_FlowChange2"] = getStor_ChangeAroundPeak(dfRegToUnregCoincident["Date_RegToUnreg"], window_days=2)["Mossyrock_FlowChange"].values
    dfRegToUnregCoincident["Mossyrock_FlowChange3"] = getStor_ChangeAroundPeak(dfRegToUnregCoincident["Date_RegToUnreg"], window_days=3)["Mossyrock_FlowChange"].values

    #Print out correlations of candidate explanatory variables
    print(dfRegToUnregCoincident.corr()["RegToUnreg"])
    plotUnregvsObsForEveryEvent()
    slope, intercept = plotStorChangeVsPeakIncreaseRegression(dfRegToUnregCoincident)
    dfUnregPeaks = applyRegressionToGetUnregPeaks(slope, intercept)
    '''
    
    
            