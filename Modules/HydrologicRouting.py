"""
Contains code to simulate routing reaches.
All routing reaches have a method called "routeHydrograph" which will route a timeseries container.
Example usage::

    myReach = initialize the reach
    routedFlows = myReach.routeHydrograph(inputFlows)
    
Much of this code was ripped from the cRouting.py module in the NWDJyLib.
But it was written in Jython, so some modifications had to be made to make it work in Python. 
"""

import pandas as pd
import numpy as np
from scipy import interpolate
from numpy import log10
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, PercentFormatter
from matplotlib.dates import DateFormatter
import math
import copy

###############################################################################
#CONFIGURATION

AC_FT_TO_CF = 43560.0 #acre-feet to cubic feet

################################################################################
# CLASS DEFINITIONS

class NullReach:
    """
    Represents Null Routing
    """
    def __init__(self):
        return None
    
    def routeHydrograph(self, inflowList):
        """Does nothing, returns input"""
        return inflowList
    
class ModPulsReach:
    """
    Represents a Modified Puls routing reach.
    Cannot handle reach losses--just does the routing.
    Typically just init the reach, and then route a timeseries with :meth:`routeHydrograph`.
    author: Ryan Cahill.
    date modified: September 5, 2014
    
    :param list storageList: list or array of storages in storage-outflow table
    :param list dischargeList: list or array of discharges in storage-outflow table
    :param int numSubreaches: the number of subreaches (must be 1 or greater)
    
    As the number of routing steps (subreaches) increases, the amount of attenuation decreases. 
    The maximum attenuation corresponds to one step. This is opposite of SSARR.
    """
    def __init__(self, storageList, dischargeList, numSubreaches) :
        self.numSubreaches = int(numSubreaches)
        self.storages = storageList
        self.discharges = dischargeList
        self.timeStepSeconds = 0 # compute interval in seconds (set depending on tsc)
        self.workingDischarge = [0 for i in range(self.numSubreaches)] # array to keep track of discharge for each subreach
        self.finalSubreachOutflow = [0 for i in range(self.numSubreaches)] # array to keep track of final outflow for each subreach
        self.storIndTbl = None # PairedValuesExt (x is flow, y is storage indication)
        self.storIndTblReverse = None # PairedValuesExt (x is storage indication, y is flow)
    
    def setTimeStepSeconds(self, timeStepSeconds):
        """Sets time step for routing and initializes storage indication table. Must be done before routing"""
        self.timeStepSeconds = float(timeStepSeconds)
        self.initStorIndicationTbl()
    
    def initStorIndicationTbl(self):
        """
        Initialize storage indication table.
        Storage indication table converts the storage values to cfs, based
        on the time step so that the storage and flow data can be easily added.
        timeStepSeconds and numSubreaches must already be defined .      
        Assumes the units of storage are ac-ft and outflow are cfs.
        """
        if self.timeStepSeconds == 0: 
            raise AssertionError("Need to set the timestep of the Modified Puls Routing")
            return False
        storages = self.storages
        outflows = self.discharges
        # manipulate the storage indication values
        storIndications = []
        for i in range(len(storages)):
            storIndications.append(0.5*outflows[i] + storages[i]*43560./self.timeStepSeconds/self.numSubreaches)     
        self.storIndTbl = interpolate.interp1d(outflows, storIndications)
        # generate the table in reverse as well
        self.storIndTblReverse = interpolate.interp1d(storIndications, outflows)
    
    def routeHydrograph(self, inflowList) :
        """
        Method to route a whole hydrograph using Modified Puls.
        inflowList is just an array or list of the inflows to the reach.
        The timestep should have already been defined before calling this method.
        Return is a list of routed flows
        """
        flowsRouted = list(inflowList) # routed flows
        # set the initial flows
        # qInit = inflowList[0]
        qInit = inflowList.iloc[0] # updated by sam
        self.finalSubreachOutflow = [qInit for i in range(self.numSubreaches)]
        self.workingDischarge = [qInit for i in range(self.numSubreaches)]
        for j in range(len(inflowList)) :
            q1 = inflowList.iloc[j]
            if j == 0 :
                # if the first timestep, just use first flow value
                q0 = inflowList.iloc[j]  
            else :
                q0 = inflowList.iloc[j-1]
            # Route the flows
            flowsRouted[j] = float(self.routeOneStep(q0, q1))
        return flowsRouted                                 
        
    def routeOneStep(self, qin0, qin1) :
        """
        Route one timestep.
        The following variables must already be defined:
        
        * workingDischarge:
        * finalSubreachOutflow:
        * storIndTbl: the time step is wrapped up in this table, don't need to know timestep
              
        :param float qin0: reach inflow for previous timestep
        :param float qin1: reach inflow for current timestep
        """
        if self.timeStepSeconds == 0: 
            raise AssertionError("Need to set the timestep of the Modified Puls Routing")
        # Loop through each subreach
        for i in range(self.numSubreaches):
            # average the beginning and end of period outflows
            qin_avg = .5 * (qin0 + qin1)
            # get beginning of period working discharge (outflow from the subreach from previous timestep)
            d0 = self.workingDischarge[i]
            # calculate storage indication
            storageIndication0 = self.storIndTbl(d0)
            # compute storage for current step
            storageIndication1 = max(0, storageIndication0 - d0 + qin_avg)
            # Interpolate outflow (working discharge) from storage indication
            d1 = max(0, self.storIndTblReverse(storageIndication1))
            # save current subreach values for next time step
            self.workingDischarge[i] = d1
            # transfer values to do next subreach
            qin0 = self.finalSubreachOutflow[i]
            qin1 = d1               
            self.finalSubreachOutflow[i] = d1
        return d1

class MuskingumReach:
    """
    Represents a Muskingum routing reach.
    Cannot handle reach losses--just does the routing.
    Typically just init the reach, and then route a timeseries with :meth:`routeHydrograph`.
    author: Ryan Cahill.
    date modified: September 12, 2014
    
    :param float muskingumK: Slope of the storage outflow relationship (hours)
    :param float muskingumX: Weighting factor on inflow. 
          If X=0, inflow doesn't affect storage and there is maximum attenuation. 
    :param int numSubreaches: the number of subreaches (must be 1 or greater)
    """
    def __init__(self, muskingumK, muskingumX, numSubreaches, timestepHrs=0):
        assert 0 <= muskingumX <= 0.5 #x cannot be above 0.5
        self.numSubreaches = numSubreaches
        self.muskingumK = muskingumK
        self.muskingumX = muskingumX
        self.initialSubreachOutflow = [0 for i in range(numSubreaches)] # array to keep track of initial outflow for each subreach
        self.finalSubreachOutflow = [0 for i in range(numSubreaches)] # array to keep track of final outflow for each subreach
        self.timestepHrs = timestepHrs #Can supply timestep at reach creation or when supplying a hydrograph
        
    def routeHydrograph(self, inflowList, timestepHrs = None):
        """
        Method to route a whole TimeSeriesContainer.
        The TimeSeriesContainer must have the "interval" property defined (minutes in each timestep).
        """
        flowsRouted = list(inflowList) # routed flows
        # Define the storage indication table (if not already defined)
        if timestepHrs:
            self.timestepHrs = timestepHrs
        if self.timestepHrs <= 0:
            raise AssertionError("Need to feed in a hydrograph with regular interval data defined")
            return None
        dtHours = self.timestepHrs #timestep in hours (routing equations work in hours)
        # set the initial flows
        qInit = inflowList.iloc[0]
        self.finalSubreachOutflow = [qInit for i in range(self.numSubreaches)]
        self.initialSubreachOutflow = [qInit for i in range(self.numSubreaches)]
        
        for j in range(len(inflowList)):
            q1 = inflowList.iloc[j]
            if j == 0 :
                # if the first timestep, just use first flow value
                q0 = inflowList.iloc[j]  
            else :
                q0 = inflowList.iloc[j-1]
            # Route the flows
            flowsRouted[j] = self.routeOneStep(dtHours, q0, q1)
        return flowsRouted                                 
        
    def routeOneStep (self, deltaT, q0, q1):
        """
        Route one timestep.
        The following variables must already be defined:
              
        * initialSubreachOutflow
        * finalSubreachOutflow
        
        :param int deltaT: timestep (hours)
        :param float qin0: reach inflow for previous timestep
        :param float qin1: reach inflow for current timestep
        """
        muskingumX = self.muskingumX
        muskingumK = self.muskingumK
        numberSubreaches = self.numSubreaches
        subK = muskingumK/float(numberSubreaches)
        # calculate muskingum coeffients
        c0 = (-subK*muskingumX + .5*deltaT)/(subK-subK*muskingumX + .5*deltaT)
        c1 = (+subK*muskingumX + .5*deltaT)/(subK-subK*muskingumX + .5*deltaT)
        c2 = (subK-subK*muskingumX - .5*deltaT)/(subK-subK*muskingumX + .5*deltaT)
        # define initial estimates
        qin0 = q0
        qin1 = q1
        # Loop through each subreach
        for i in range(numberSubreaches):
            qout0 = self.initialSubreachOutflow[i]
            qout1 = c0*qin1 + c1*qin0 + c2*qout0  #apply the Muskingum coefficients
            qin1 = qout1
            self.finalSubreachOutflow[i] = qout1
            # outflow0 is inflow0 for next subreach
            qin0 = qout0
        # save values for next time step
        for i in range(numberSubreaches):
            self.initialSubreachOutflow[i] = self.finalSubreachOutflow[i]
        # return reach outflow for this step
        return qout1

class SsarrReach:
    """
    Class to simulate a SSARR routing reach 
    Originally written by John McCoskery, modified by Ryan Cahill
    Recommended usage: Initialize the reach, 
    apply one of the build methods, 
    initialize the starting flows
    to use :func:`buildReach` rather than creating this object directly.
    As the number of phases (subreaches) increases with all other factors constant the
    peak decreases and the total time of storage increases.This is opposite of Modified Puls subreaches. 
    """
    def __init__(self, timestepHrs = 0):
        self.numSubreaches = 0 #Also known as "phases" in SSARR
        self.ncoeff = 0  # SSARR n-coefficient
        self.kts = 0  # SSARR Time of Storage coefficient
        self.timestepHrs = timestepHrs  # compute interval (always in hours)
        self.flow_ts_Table = None  # flow-Ts table (interpolate.interp1d object)
        self.usesTable = False                    # flag if table is used or not
        self.subreachOutflows = []  # list with outflow for each subreach. Updated each timestep
    
    def clone (self) :
        """clone the object (useful if doing trial releases but don't want to reset subreachOutflows)"""
        return copy.copy(self)
        
    def buildWithTable(self, flow_ts_Table, numSubreaches):
        """method to initialize w/ outflow-time of storage table"""
        self.flow_ts_Table = flow_ts_Table
        self.usesTable = True
        self.numSubreaches = numSubreaches
    
    def buildWithKTS(self, numSubreaches, n, kts):
        """method initialize w/ kts method
        Supply the number of subreaches, n coefficient, and KTS coefficient"""
        self.numSubreaches = numSubreaches
        self.ncoeff = n
        self.kts = kts
    
    def init_flows(self, initialFlow):
        """initialize starting flows for first step"""
        self.subreachOutflows = [initialFlow for i in range(self.numSubreaches)] 
    
    def routeOneStep(self, qIn1, qIn2):
        """
        Route flows for one step (sets subreachOutflows as the resulting outflows)
        
        :param float qIn1: reach inflow for previous timestep
        :param float qIn2: reach inflow for current timestep
        """
        #get the flow value for lookup on Ts, average starting flow over all subreaches
        qTS = sum(self.subreachOutflows)/float(self.numSubreaches)
        # set a floor on the flow for interpoolation, for Ts = KTS or close to bottom of QR-Ts table
        # this is what ResSim does to get around (-) subreachOutflows values
        qTS = max(1.0, qTS)  
        # compute time of storage value, Ts
        if self.usesTable: #uses w/ qr-ts table
            ts = self.flow_ts_Table(qTS) 
        else : # w/ kts method
            ts = self.kts / qTS**self.ncoeff
        q_Outflows = [] #list of resulting outflows for each subreach
        if ts < 0.05 : # if ts is very small, don't route, set outflows as inflows
            q_Outflows = [qIn2 for i in range(self.numSubreaches)]
        else :
            q_Outflows = list(self.subreachOutflows)# copy subreachOutflows from previous step to next step as starting point. 
            
            if ts < self.timestepHrs/2. :
                # if ts is smaller than timestepHrs/2, need to break into sub-period
                nper = int(self.timestepHrs / 2. / ts) + 1
            else :
                nper = 1
            #TODO this step of dt and tsr is different in the original McCoskery script file and ResSim SsarrRouting.java file.
            dt = self.timestepHrs/2./float(nper) #Sub-period timestep 
            tsr = dt / (ts + dt)  # compute tsr
            qInAvgIncr = (qIn2 - qIn1) / float(nper)   # compute incremental average flow over sub-interval
            pqi = qIn1 - qInAvgIncr / 2  # compute average qin to reach
            for j in range(nper): # compute phase flows for each substep
                pqi = pqi + qInAvgIncr
                qi = pqi   # initialze inflow for first sub-reach
                for i in range(self.numSubreaches):
                    dq = qi - q_Outflows[i]  # compute difference between current qi and previous qi for sub-reach
                    dq = dq * tsr  # multiply by tsr for sub interval
                    qi = q_Outflows[i] + dq  # adjust initial qi to sub-reach
                    q_Outflows[i] = qi + dq  # compute final flow for sub-reach
        # set new subreachOutflows to computed outflows
        self.subreachOutflows = q_Outflows
        # return the reach outflow (from the last subreach)
        return q_Outflows[-1]
    
    def routeHydrograph (self, inflowList, allowNegatives=False) :
        """
        Method to route a whole hydrograph
        
        :param list inflowList: The upstream inflow hydrograph
        :param boolean allowNegatives: (Optional, default is False). If False,
            any of the input hydrograph values that are negative will be set to 1.
            If True, any negative values in the input will be left as is and routed.
            When ResSim does SSARR routing, it uses False. 
            SSARR routing will still function with negative flows, but it is not 
            applied in ResSim.
        """
        if self.timestepHrs == 0:
            raise AssertionError("Need to set the timestep of the SSARR Routing first")
        qrr = list(inflowList) # routed flow array
        for j in range(len(inflowList)) :
            if j == 0 :
                #initialize subreachOutflows array
                self.init_flows(inflowList[j])
            else :
                #Route the flows
                if allowNegatives:
                    self.routeOneStep(inflowList[j-1], inflowList[j])
                else:
                    #floor the values to 1 cfs, don't allow negatives
                    self.routeOneStep(max(1.,inflowList[j-1]), max(1.,inflowList[j]))
            # store final phase flow to qrr array
            qrr[j] = self.subreachOutflows[-1]
        return qrr

################################################################################
# FUNCTION DEFINITIONS
       
def convertSSARRtoModPuls(SsarrReachObj, discharges):
    """
    Attempts to develop a storage-release table for Modified Puls that mimics
    the routing that would occur under the SSARR routing reach provided.
    Assumes that the desire is a Modified Puls table with 1 subreach
    
    :param SsarrReachObj SsarrReach: The SSARR routing reach to replicate (already has timestep embedded in it)
    :param list discharges: A list (or numpy array) of all the discharges desired for the output storage-discharges table, ascending order
    
    Returns a pandas dataframe with the Modified Puls discharge-storage table, with units in cfs and acre-feet
    """
    
    #First, do SSARRR routings for the entire set of previous release (O1) and mean inflow (Im)
    #Previous releases are columns, and mean inflows are rows
    numDischarges = len(discharges)
    outflowArray = np.zeros((numDischarges,numDischarges))
    for j, meanInflow in enumerate(discharges):
        for i, prevDischarge in enumerate(discharges):
            SsarrReachObj.init_flows(prevDischarge)
            outflow = reachObj.routeOneStep(prevDischarge,2*meanInflow - prevDischarge)
            if outflow < 0: outflow = np.nan
            outflowArray[i,j] = outflow
    dfSSARR = pd.DataFrame(outflowArray, index=discharges, columns=discharges)
    
    #Now, try to figure out storages in ac-ft for Modified Puls that works
    storageArray = np.zeros((numDischarges,numDischarges))
    timestepSeconds = SsarrReachObj.timestepHrs * 60 * 60
    #For each row (previous discharge), loop through each column to ramp up inflows
    for i, prevDischarge in enumerate(discharges):
        tempDischarges = [discharges[0]]
        tempStorages = [0.] #Initialize to 0 ac-ft at minimum discharge
        O1 = prevDischarge
        #Establish starting storage in ac-ft
        if i == 0: #first row
            S1 = 0
        else: #interpolate from the table that we are building
            S1 = dischargeStorTbl(O1)
        for j, meanInflow in enumerate(discharges):
            if i == 0 and j == 0: continue #no need to do anything, we already initialized it
            O2 = dfSSARR.iloc[i,j] #final outflow from SSARR routing with same previous discharge and mean inflow. This is what we are trying to match
            if np.isnan(O2): continue #skip it
            S2 = timestepSeconds/AC_FT_TO_CF*(meanInflow + (S1*AC_FT_TO_CF/timestepSeconds - O1/2) - O2/2) #ending storage (from Modified Puls Equation)
            storageArray[i,j] = S2
            tempDischarges.append(O2)
            tempStorages.append(S2)
        #set up a temporary discharge-storage table for the previous release values
        #If any extrapolation is required, it would be very minor and on the extreme edges, so it's probably OK
        dischargeStorTbl = interpolate.interp1d(tempDischarges, tempStorages, bounds_error = False, fill_value='extrapolate')
        tempTbl = pd.DataFrame()
        tempTbl["Discharge"] = tempDischarges
        tempTbl["Storage"] = tempStorages
        #if i == 0: break #TODO

    #Now, set up the Modified Puls table
    dfModPulsTbl = pd.DataFrame()
    dfModPulsTbl["Discharge"] = discharges
    dfModPulsTbl["Storage"] = np.diag(storageArray)
    
    return dfModPulsTbl


#####  Main  ###################################################################
# MAIN CODE

if __name__ == "main" or __name__ == "__main__":
    # Code starts here when called as script.  Otherwise can be loaded as a module
    #Reach is the south santiam from waterloo to the confluence
    storageList = [0,
        38,
        93,
        149,
        316,
        455,
        593,
        871,
        1149,
        1705,
        2816,
        4205,
        5594,
        6983,
        8372,
        11150,
        13928,
        16706,
        19483,
        22261,
        27817,
        38929,
        47263,
        55597,
        70080,
        116800
        ]
    dischargeList = [0,
        65,
        160,
        255,
        541,
        779,
        1016,
        1492,
        1968,
        2919,
        4822,
        7200,
        9578,
        11957,
        14335,
        19092,
        23849,
        28606,
        33362,
        38119,
        47632,
        66660,
        80930,
        95200,
        120000,
        200000
        ]
    
    
    timestepHours = 1
    reachObj = SsarrReach(timestepHrs = timestepHours)
    reachObj.buildWithKTS(2, 0.1, 5)
    #reachObj.init_flows(1000)
    #qOut = reachObj.routeOneStep(1000,10000)
    
    #discharges = np.logspace(log10(100), log10(100000), num=20)
    
    #Convert the SSARR reach to Modified Puls
    #TODO this works perfectly if the number of subreaches in SSARR is 1. It doesn't work perfectly if number of subreaches is > 1
    #When number of subreaches > 1 in SSARR, hysteresis occurrs in the discharge-storage function. 
    discharges = np.round(np.linspace(100, 50000, num=20))
    discharges = [50,100,250,500, 1000,2000,3000,4000,4650,5500,6780,7500,9500,11500,13500,15500,17500,19300,25000, 180000]
    SsarrReachObj = reachObj
    dfModPulsTbl = convertSSARRtoModPuls(reachObj, discharges)
    #Data point 1: good match to SSARR with n=0.1, KTS=20, phases=2 by using 2 subreaches and multiplying storages by 0.12
    #Data point 1: good match to SSARR with n=0.1, KTS=5, phases=2 by using 2 subreaches and multiplying storages by 0.35
    
    #Route a test hydrograph through to validate it
    dfInflows = pd.read_csv("../Data/CandidateHydrographs.csv", parse_dates=True, index_col='DateTime')
    #Convert it to the time interval of interest for routing
    dfHydrographs = dfInflows.resample("%dH" %timestepHours).mean()
    dfHydrographs["SSARR"] = reachObj.routeHydrograph(dfHydrographs["Inflow"])
    '''
    reachObj.numSubreaches = 2
    dfHydrographs["SSARR_2_subreaches"] = reachObj.routeHydrograph(dfHydrographs["Inflow"])
    '''
    
    reachModPuls = ModPulsReach(dfModPulsTbl["Storage"], dfModPulsTbl["Discharge"], 1)
    reachModPuls = ModPulsReach(dfModPulsTbl["Storage"]*0.35, dfModPulsTbl["Discharge"], 2)
    #reachModPuls = ModPulsReach(dfModPulsTbl["Storage"]*0.12, dfModPulsTbl["Discharge"], 2) #TODO why does it work better with 2 subreaches and dividing storage by 2?
    #reachModPuls = ModPulsReach(dfModPulsTbl["Storage"], dfModPulsTbl["Discharge"], 5)
    reachModPuls.setTimeStepSeconds(timestepHours*60*60)
    dfHydrographs["ModPuls"] = reachModPuls.routeHydrograph(dfHydrographs["Inflow"])
    
    '''
    for ratio in np.arange(0.000045,0.000075,0.000005):
        for nSubreaches in range(1,4):
            reachModPuls = ModPulsReach(dfModPulsTbl["Storage"]*ratio, dfModPulsTbl["Discharge"], nSubreaches) #TODO why does it work better with 2 subreaches and dividing storage by 2?
            reachModPuls.setTimeStepSeconds(timestepHours*60*60)
            dfHydrographs["ModPuls-%0.6f, %s subreaches" %(ratio, nSubreaches)] = reachModPuls.routeHydrograph(dfHydrographs["Inflow"])
    '''
    #dfHydrographs.plot()
    
    
    ###############################################################################
    #Check to see how different the routing results are when using daily average data vs hourly data
    #Question: should we reduce the number of subreaches when we move from hourly to daily?
    #Answer: probably not. We are way off using daily flows anyways, why fool ourselves into thinking we are precise?
    plt.close('all')
    fig, ax = plt.subplots(1)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(DateFormatter("%b-%d"))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    
    dfDaily = dfInflows.resample("1D").mean()
    reachDaily = ModPulsReach(storageList, dischargeList, 5)
    reachDaily.setTimeStepSeconds(60*60*24)
    dfDaily["Routed"] = reachDaily.routeHydrograph(dfDaily["Inflow"])
    reachDaily = ModPulsReach(storageList, dischargeList, 1)
    reachDaily.setTimeStepSeconds(60*60*24)
    dfDaily["Routed-1Sub"] = reachDaily.routeHydrograph(dfDaily["Inflow"])
    
    dfHourly = dfInflows.resample("1H").mean()
    reachHourly = ModPulsReach(storageList, dischargeList, 5)
    reachHourly.setTimeStepSeconds(60*60)
    dfHourly["Routed"] = reachHourly.routeHydrograph(dfHourly["Inflow"])
    
    df3Hr = dfInflows.resample("3H").mean()
    reach3Hr = ModPulsReach(storageList, dischargeList, 5)
    reach3Hr.setTimeStepSeconds(60*60*3)
    df3Hr["Routed"] = reach3Hr.routeHydrograph(df3Hr["Inflow"])
    reach3Hr = ModPulsReach(storageList, dischargeList, 2)
    reach3Hr.setTimeStepSeconds(60*60*3)
    df3Hr["Routed-2Subs"] = reach3Hr.routeHydrograph(df3Hr["Inflow"])

    #dfHourly.plot(ax=ax)
    
    #dfDaily.plot(ax=ax)
    
    ax.plot(dfHourly.index, dfHourly["Inflow"], label="Hourly Inflow", color="green")
    ax.step(dfDaily.index, dfDaily["Inflow"], where="post", label="Daily Inflow", color="green", ls="--")
    ax.plot(dfHourly.index, dfHourly["Routed"], label="Hourly Routed (5 subreaches)", color="orange")
    ax.step(dfDaily.index, dfDaily["Routed"], where="post", label="Daily Routed (5 subreaches)", color="blue")
    ax.step(dfDaily.index, dfDaily["Routed-1Sub"], where="post", label="Daily Routed (1 subreaches)", color="tab:gray")
    ax.plot(df3Hr.index, df3Hr["Routed"], label="3-hour Routed (5 subreaches)", color="tab:purple")
    ax.plot(df3Hr.index, df3Hr["Routed-2Subs"], label="3-hour Routed (2 subreaches)", color="tab:brown")
    dfHrtoDaily = dfHourly.resample("1D").mean()
    ax.step(dfHrtoDaily.index, dfHrtoDaily["Routed"], where="post", label="Hourly Routed, daily average (5 subreaches)", color="orange", ls="--")
    ax.set_ylabel("Flow (cfs)")
    ax.legend()
    