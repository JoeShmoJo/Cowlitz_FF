"""
Module to deal with DSS file reading/writing.
This module can use Python 3 classes (e.g. pandas dataframes)
The objects and methods in this module are a mirror
of those in the HEC Java library.
This is done to make sure the code that calls this module
can be run in either Jython (using the HEC libraries) or Python 3 (without HEC libraries)

Typical usage:

#Open file
dssFile = HecDss.open("SomeFile.dss")
#Read data
df = dssFile.readDF("//DWORSHAK/FLOW-OUT//1DAY/TEST/")
unitStr = df.units
typeStr = df.type
#Write data
dssFile.writeDF(df, "//DWORSHAK/FLOW-OUT//1DAY/COPIED/", unitStr, typeStr)
#Close file
dssFile.close()
"""

################################################################################
# IMPORTS

import sys
import datetime
#import copy
import os
import re
from datetime import timedelta
from utilsTime import toPyDatetime
from utilsIO import ensure_dir

if sys.version_info[0] == 3:
    import numpy as np
    import pandas as pd
    try:
        from pydsstools.heclib.dss import HecDss as DSS_External
        from pydsstools.core import TimeSeriesContainer as TSC_External
    except ImportError:
        raise AssertionError("Need to install the pydsstools library")
else:
    #In Jython, just keep it arbitrary to make sure the functions in this module can still be imported (but not the classes)
    DSS_External = object
    TSC_External = object

################################################################################
# STATIC INPUT

INTERVAL_DICT = {"15MIN" :   15,
                 "1HOUR" :   60, 
                 "3HOUR" :  180,
                 "4HOUR" :  240,
                 "6HOUR" :  360,
                 "12HOUR":  720,
                 "1DAY"  : 1440
                 }

VALID_TYPE_STR = ["PER-AVER", "INST-VAL"]

#DSS times are usually stored as minutes past this fixed base date
BASE_DATETIME_DSS = datetime.datetime(1899,12,31,0,0) #This is a ._julian_base_date['julianBaseDate'] of 0
BASE_DATE_DSS = datetime.date(1899,12,31) 

################################################################################
# CLASS DEFINITIONS

class HecDss():
    """
    Class for reading/writing data from a DSS file. 
    This object represents a single DSS file.
    This is a mirror of hec.heclib.dss.HecDss
    which is a subclass of hec.hecmath.DSSFile
    
    Typical usage:
        dssFile = HecDss.open("myDSSFile.dss")
        df = dssFile.readDF("//DWORSHAK/FLOW-OUT//1DAY/TEST/")
        unitStr = df.units
        typeStr = df.type
        #Write data
        dssFile.writeDF(df, "//DWORSHAK/FLOW-OUT//1DAY/COPIED/", unitStr, typeStr)
        #Close file
        dssFile.close()
    """
    def __init__(self):
        self.fileObj = None #pydsstools.heclib.dss.HecDss object
        self.fileName = "" #The full filename string of the file
        self.startTime = None #datetime.datetime (optional)
        self.endTime = None #datetime.datetime (optional)
    
    #def open(self, fileName, startTime=None, endTime=None):
    @classmethod
    def open(cls, fileName, startTime=None, endTime=None):
        """
        Creates an instance of the class by opening a DSS file. 
        If it does not exist, it will create a blank one and open it.
        Optional arguments: startTime and endTime (datetime.datetime objects)
        
        This method is a mirror of theJava method

        #TODO handle overloaded .open methods in Java (timewindow option)
        """
        clsObj = cls() #Calls __init__ and returns an instance of the current class
        clsObj.openFile(fileName, startTime, endTime)
        return clsObj

    def openFile(self, fileName, startTime=None, endTime=None):
        """
        Opens the DSS file. 
        If it does not exist, it will create a blank one and open it.
        Optional arguments: startTime and endTime (datetime.datetime objects)
        """
        if not self.fileObj is None:
            #The DSS file is already open
            if fileName == self.fileName:
                #It's the same one
                return self
            else:
                #If we want to open a new file, have to close the existing file
                self.close()
        self.fileName = fileName
        if os.path.exists(fileName):
            #Take whatever version exists
            self.fileObj = DSS_External.Open(fileName) #pydsstools.heclib.dss.HecDss
        else:
            #When creating new files, do it in version 6 for stability
            ensure_dir(fileName)
            self.fileObj = DSS_External.Open(fileName, version=6) #pydsstools.heclib.dss.HecDss
        #You can't actually open the file with a time window, but save the start/end times
        #so they will be used as default when reading time series data
        if not startTime is None:
            self.startTime = toPyDatetime(startTime)
            self.endTime = toPyDatetime(endTime)
        return self
    
    def readTSC(self, pathnameStr, timeStart=None, timeEnd=None):
        """
        Reads the given pathname from the DSS file
        Returns a TimeSeriesContainer object
        """
        df = self.readDF(pathnameStr,timeStart, timeEnd)
        unitStr = df.units
        typeStr = df.type
        tsc = createTSC(df.index, df["value"].values, unitStr, typeStr, pathname=pathnameStr)
        return tsc
    
    def readDF(self, pathname, timeStart=None, timeEnd=None, dropMissing=False):
        """
        Reads the given pathname from the DSS file

        :param str pathname: The pathname
        :param DateTime timeStart: the start time to begin retrieving.
        :param DateTime timeEnd: end time of the extract
        :param boolean dropMissing: If true, missing values will be dropped. If false, they will be kept as na
        
        Returns : DataFrame with index as the TimeStamp, and column named "value" as the data
            If the input pathname is totally invalid, raise an error
            If the pathname is fine, but there is just no data, return an empty DataFrame
            There's also two attributes attached to this dataframe
            Be careful, the moment you run any other functions (e.g. dropna), these will get wiped away
                1. units (e.g. "cfs")
                2. type (e.g. "PER-AVER")
            Best practice is to call this method and then immediately get the units and type if you need them
        """
        dfEmpty = pd.DataFrame(columns=["value"])
        dfEmpty.units = "N/A"
        dfEmpty.type = "N/A"
        pathname = pathname.upper()
        if pd.isna(pathname) or not DSSPathString().isValid(pathname): #not valid
            raise Exception("Read failed: %s" %(pathname))
        pathname = removeDPart(pathname)
        #deal with time windows
        startTime = self.startTime #default
        endTime = self.endTime #default
        if not timeStart is None: #use the input arguments if supplied
            startTime = toPyDatetime(timeStart)
            endTime = toPyDatetime(timeEnd)
        if startTime is None:
            ts = self.fileObj.read_ts(pathname)
        else:
            ts = self.fileObj.read_ts(pathname, window=(startTime, endTime))

        if ts is None or ts.values is None:# or (-901 in ts.values and "IR-CENTURY" in pathname): 
            print("Failed to read pathname: %s" %pathname)
            return dfEmpty
        elif min(ts.values) >=-902. and max(ts.values) <= -901.:
            #-901/902 is missing to DSS. 
            # All values are missing. 
            print("Failed to read pathname: %s" %pathname)
            return dfEmpty
        #dtIndex = ts.pytimes #This works, but it takes forever #pydsstools/src/time_series.pyx

        #BACKWARDS COMPATABILITY: The formatting of the pydsstools changed subsantially. A few checks here will determine how the object will be handled.
        
        if isinstance(ts.times, list): #This is the old format, where the times attribute is a list of integers
            startPyDate = ts.startPyDateTime
            timeZero = ts.times[0] #Sometimes this object is a dictionary, we need it in a list in order to pull the first value
            if ts.granularity == 60: #Minute granularity
                if ts.interval <= 0: #Irregular data, uses a different system
                    if ts._julian_base_date["julianBaseDate"] != 0: 
                        raise AssertionError("Julian Base date isn't 0, haven't encountered this before")
                    #ts.startPyDateTime is not accurate with irregular data, it just returns the input start time of the time window supplied rather than the actual first time. 
                    dtIndex = [BASE_DATETIME_DSS + timedelta(minutes=(tt)) for tt in ts.times]
                else:
                    #pydsstools/src/time_series.pyx says to use .startJulianDate, which is the offset number of days since 31Dec1899 00:00
                    #But that's the same thing as the startPyDate, so no need to worry about that
                    dtIndex = [startPyDate + timedelta(minutes=(tt-timeZero)) for tt in ts.times]
            else:
                raise AssertionError("Invalid granularity: %" %ts.granularity)
            dfOut = pd.DataFrame(index = dtIndex, columns=["value"])
            dfOut["value"] = ts.values
        elif isinstance(ts.times, dict): #This is the latest format, where the times attribute is a dictionary. This works for regular data intervals.
            startPyDate = ts.startPyDateTime
            if ts.granularity == 60: #Minute granularity
                values_col = ts.values
                times_col = pd.date_range(start=startPyDate, periods=len(values_col), freq='{}s'.format(ts.interval))
                dfOut = pd.DataFrame(values_col, columns=['value'])
                dfOut = dfOut.set_index(times_col)
            else:
                raise AssertionError("Invalid granularity: %" %ts.granularity)
        else:
            raise AssertionError("pydsstools formatting invalid.")
        if dropMissing:
            dfOut = dfOut[dfOut.value != -901] #drop the missing values
            dfOut = dfOut[dfOut.value != -902] #drop the missing values
        else:
            dfOut = dfOut.replace(-901, np.nan)
            dfOut = dfOut.replace(-902, np.nan)
        dfOut = dfOut.astype(np.float64, copy=False)
        dfOut.units = ts.units
        dfOut.type = ts.type
        return dfOut
    
    def hecTime_to_datetime(self, hecTime_object):
        '''
        DEPRECATED
        TODO, function to be removed
        Quick function to convert the hectime object to datetime. The hectime object is from the pydsstools library.
        #Writing a dictionary that takes the abreviated month as a key and the integer as the value
        '''
        month_dict = {
        "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
        "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12}

        #Parsing through the date (dateString method)
        formatted_date = hecTime_object.dateString()
        year = formatted_date[-4:]
        month_string = formatted_date[2:-4] #returns the month as a string(apr, aug, etc)
        day = formatted_date[:2]
        month = month_dict[month_string.lower()] #Gets the month integer value from the dictionary

        #Parsing through the time (timeString method)
        formatted_time = hecTime_object.timeString()
        hour, minute, second = formatted_time.split(':')
        
        #Creating a datetime object
        new_dt = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))

        return new_dt

    def _write(self, tsc):
        """
        Writes a TimeSeriesContainer to dss
        This tsc should already be well formed
        This is an internal function because we want calling functions to use .writeTSC or .writeTS or .write
        """
        if self.fileObj is None:
            if self.fileName != "":
                #The user must have forgotten to call .open or might have closed the dss file already
                #Go ahead and open it
                self.openFile(self.fileName)
            else:
                raise AssertionError("DSS file is not open, can't write data")
        #Run checks on tsc to see if it is ready to write
        tsc.checkReadyToWrite()
        #Write it out
        self.fileObj.put_ts(tsc)
        return True
    
    def writeTSC(self, tsc):
        """
        Writes a TimeSeriesContainer to DSS.
        This function can handle regular date, irregular data, and does fine if the data has gaps in it too
        Returns True if write was successful, false if not
        """

        if tsc.interval <= 0:
            #Irregular data, just write it all out
            self._write(tsc)
            return True
        #Regular data
        is_continuous = tsc.isContinuous()
        if is_continuous:
            #Write it all out as one big continous chunk
            self._write(tsc)
        else:
            #May be regular data (e.g. 1HOUR), but it's discontinous and full on missing some time periods.
            #Change output time interval, broken up indexes, and frequency depoending on 
            pathname = tsc.pathname
            if '1HOUR' in pathname:
                idx = "MS"          #this is used to break up the data into monthly chunks of data to write to DSS for hourly data
                frequency = "H"     #used to create continuous time stamps to merge data into continuous data series
            elif '1DAY' in pathname:
                idx = "YS"          #this is used to break up the data into yearly chunks of data to write to DSS for daily data
                frequency = "D"     #used to create continuous time stamps to merge data into continuous data series
            else:
                raise AssertionError("Unsupported timestep for pathname: %s" %pathname)
            startDate = getPyDateTimeFromDSSTimeMinutes(tsc.times[0])
            endDate = getPyDateTimeFromDSSTimeMinutes(tsc.times[-1])

            #Create a dataframe with continous, regular interval times without any gaps. Then, fill in the values from the actual data
            df = pd.DataFrame(index=tsc.times)
            df["value"] = tsc.values
            dfContinuous = pd.DataFrame(index=pd.date_range(startDate, endDate, freq = frequency))
            dfContinuous = pd.merge(dfContinuous, df, how='left', left_index=True, right_index=True)
            
            #we want to break datastream up in the following only if we have more than one month of hourly data, or more than 2 years of daily data
            #thus the thresholds
            one_timestep_threshold = 31 if "HOUR" in pathname else  2  
            divisor = 24 if "HOUR" in pathname else 365 
            if len(tsc.times)/divisor > one_timestep_threshold:
                #Break the continuous frame into consecutive periods:
                #[start -> 1st interior period start), [period start -> next), ..., [last period start -> end]
                #(The previous version of this loop dropped the final stretch, wrote degenerate
                #middle chunks, never updated tsc.times, and called _write on the wrong object.)
                periodStartDates = pd.date_range(startDate, endDate, freq = idx)
                boundaries = [startDate]
                for p in periodStartDates:
                    if startDate < p < endDate:
                        boundaries.append(p)
                boundaries.append(endDate)
                for j in range(len(boundaries) - 1):
                    dfPeriod = dfContinuous.loc[boundaries[j]:boundaries[j + 1]]
                    if j < len(boundaries) - 2 and len(dfPeriod) > 1:
                        #drop the shared boundary timestep so consecutive chunks don't overlap
                        dfPeriod = dfPeriod.iloc[:-1]
                    if dfPeriod.isna().sum().values[0] == len(dfPeriod):
                        #Every single value is missing, skip it
                        continue
                    tsc.startDateTime = dfPeriod.index[0].strftime("%d%b%Y %H:%M:00") # "15JUL2019 19:00:00"
                    tsc.numberValues = len(dfPeriod)
                    tsc.times = list(dfPeriod.index.to_pydatetime())
                    tsc.values = dfPeriod[dfPeriod.columns[0]].values
                    self._write(tsc)
            else:
                #Still need to merge to continuous data for shorter periods, this is important for writing shorter indexed peak time windows
                #Idx = pd.date_range(startDate, endDate)
                dfPeriod = pd.DataFrame(index=pd.date_range(startDate, endDate, freq = frequency))
                dfPeriod = pd.merge(dfPeriod, dfContinuous, how='left', left_index=True, right_index=True)
                if dfPeriod.isna().sum().values[0] == len(dfPeriod):
                    #Every single value is missing, skip it
                    return False
                tsc.startDateTime = dfPeriod.index[0].strftime("%d%b%Y %H:%M:00") # "15JUL2019 19:00:00"
                tsc.numberValues = len(dfPeriod)
                tsc.times = list(dfPeriod.index.to_pydatetime())
                tsc.values = dfPeriod[dfPeriod.columns[0]].values
                self._write(tsc)
        return True

    def writeDF(self, df, pathname, unitStr, typeStr):
        """
        Writes a DataFrame with time series data to DSS.
        Basically the same thing as .writeSeries
        :param  pd.DataFrame    df:  time-series data to be written out. Index = datetime, only one column with values
        :return True if written, False if not
        """
        if isinstance(df, pd.Series):
            #Just in case the user accidentally gave a Series instead of a DataFrame
            return self.writeSeries(df, pathname, unitStr, typeStr)
        if len(df.columns) > 1: raise AssertionError("Only dataframes with a single column allowed")
        return self.writeSeries(df[df.columns[0]], pathname, unitStr, typeStr)

    def writeSeries(self, s, pathname, unitStr, typeStr):
        """
        Writes a Series object with time series data to DSS.
        :param  pd.Series    s:  time-series data to be written out. Index = datetime, only one column with values
        :param  str       pathname:   DSS path string to save 
        :param  str        unitStr: e.g. "cfs"
        :param  str        typeStr: e.g. "INST-VAL"
        :return True if written, False if not
        """
        if s.isna().sum() == len(s):
            #Every single value is missing, skip it
            return False
        times = s.index
        values = s.values
        tsc = createTSC(times, values, unitStr, typeStr, pathname = pathname)
        self.writeTSC(tsc)
        return True
        

    def write(self, tsc):
        """
        Mirror of method in hec.hecmath.DSSFile
        """
        self.writeTSC(tsc)

    def put(self, tsc):
        """
        Mirror of method in hec.heclib.dss.HecDss
        """
        self.writeTSC(tsc)

    def get(self, pathname, *args):
        """
        Mirror of overloaded method in hec.heclib.dss.HecDss.

        Returns a TimeSeriesContainer object

        in Java, the various calls are:
        get(java.lang.String pathname)
        get(java.lang.String pathname, boolean readEntireSet)
        get(java.lang.String pathname, java.lang.String startTime, java.lang.String endTime)
        """
        numOptionalArgs = len(args)
        if numOptionalArgs == 0:
            return self.readTSC(pathname)
        elif numOptionalArgs == 1:
            if args[0] is True:
                #Read the whole timewindow
                pathname = pathname.removeDPart()
                return self.readTSC(pathname, timeStart=None, timeEnd=None)
            else:
                #Normal
                return self.readTSC(pathname)
        elif numOptionalArgs == 2:
            #startTime, endTime specified
            return self.readTSC(pathname, timeStart=args[0], timeEnd=args[1])
        else:
            raise AssertionError("Unrecognized input arguments: %s" %args)
    
    def read(self, pathname, *args):
        """
        Same thing as .get
        """
        return self.get(pathname, args)
    
    def close(self):
        """Closes the dss file"""
        if self.fileObj is None: return None #The file isn't open yet
        self.fileObj.close()
        self.fileObj = None
        self.fileName = ""

    def done(self):
        """Same thing as .close"""
        self.close()

class TimeSeriesContainer(TSC_External):
    """
    Single time series with times and corresponding values. 
    Only real metadata is units, type, and DSS pathname string.
    This has to be a subclass of the pydsstools.core.TimeSeriesContainer
    to work properly with writing data using the pydsstools packages.
    """
    def __init__(self):
        """
        The attributes here match up with pydsstools.core.TimeSeriesContainer
        They are slighly different than the hec java class hec.io.TimeSeriesContainer
        """
        #.pathname here is identical to .fullName in the HEC java class
        #Need to duplicate the attribute so it will work in Jython as well as Python
        #The ._syncAttributes method does the work, it is to be called before writing the data
        self.fullName = "" #e.g. "//Mica/Elev//1DAY/TEMP/"
        self.pathname = "" #e.g. "//Mica/Elev//1DAY/TEMP/"
        #times are a list of datetime.datetime objects. 
        #This is different than hec.io.TimeSeriesContainer, which uses times in minutes past the BASE_DATETIME_DSS
        self.times = [] 
        self.values = [] #list of corresponding values
        self.startDateTime = None #e.g. "15JUL2019 19:00:00" (dateTimeObject.strftime("%d%b%Y %H:%M:00"))
        self.numberValues = 0
        self.units = "" # e.g. "cfs"
        self.type = "" # either "INST-VAL" or "PER-AVER"
        self.interval = -1 #If Regular Data (e.g. 1HOUR), interval in minutes. If irregular data (e.g. IR-CENTURY), -1

        #Other attributes in hec.io.TimeSeriesContainer that are not needed
        self.fileName = ""

    def setPathname(self, pathname):
        """
        Convenience method for setting the dss pathname (always uppercase)
        """
        self.pathname = pathname.upper()
        self.fullName = self.pathname
    
    def isContinuous(self):
        """
        Check whether the timeseries is continuous data (with no gaps) or not. 
        Check to see if the theoretical # of values matches the actual number of values
        returns True/False
        """
        if self.interval <= 0:
            return False
        numValuesExpected = (self.times[-1] - self.times[0]).total_seconds()/60./self.interval + 1
        is_continuous = True #True if there are no gaps in the data, False if gaps
        if numValuesExpected != len(self.values):
            is_continuous = False
        return is_continuous
    
    def _syncAttributes(self):
        """
        Make sure .fullName and .pathname are in sync, one of them may have been messed with by a calling script
        """
        if self.fullName != self.pathname:
            #Use whichever is longer
            fullNameLen = len(self.fullName)
            pathnameLen = len(self.pathname)
            if fullNameLen > pathnameLen:
                self.pathname = self.fullName
            else:
                self.fullName = self.pathname

    def checkReadyToWrite(self):
        """
        returns True if the object is formatted correctly and ready to write.
        raises an error if something is wrong
        """
        self._syncAttributes()
        if len(self.values) == 0: 
            raise AssertionError("No values in the time series")
        if not DSSPathString().isValid(self.pathname): 
            raise AssertionError("Invalid dss pathname: %s" %self.pathname)
        if len(self.values) != len(self.times): 
            raise AssertionError("Number of values: %s doesn't match number of times: %s" %(len(self.values), len(self.times)))
        if self.type not in VALID_TYPE_STR:
            raise AssertionError("Type string not valid: %s" %self.type)
        if not isinstance(self.times[0], datetime.datetime):
            raise AssertionError(".times must be datetime.datetime, they are currently: %s" %self.times[0].__class__)
        return True
    
    def clone(self):
        '''Returns a copy of the timeseriescontainer'''
        # copy.deepcopy doesn't work, have to do it the hard way
        #return copy.deepcopy(self)
        tscOut = TimeSeriesContainer()
        tscOut.fullName = self.fullName
        tscOut.pathname = self.pathname
        tscOut.times = self.times
        tscOut.values = self.values
        tscOut.startDateTime = self.startDateTime
        tscOut.numberValues = self.numberValues
        tscOut.units = self.units
        tscOut.type = self.type
        tscOut.interval = self.interval
        tscOut.fileName = self.fileName
        return tscOut

class DSSPathString():
    """Simple class that provides methods for manipulating DSS pathnames.
    A valid DSS pathname is something like:
    "//MICA/ELEV/01JAN2010/1DAY/TEST/"
    
    Mirror of hec.lang.DSSPathString

    #TODO need to mirror hec.heclib.dss.DSSPathname methods too
    """
    def __init__(self, pathname=""):
        self.pathname = "" 
        if pathname != "":
            self.setPathname(pathname)
    def checkValid(self, pathname):
        """Raises an error if the pathname isn't valid"""
        if not self.isValid(pathname):
            raise AssertionError("Something is wrong with the format of the dss path: %s" %pathname)
    def isValid(self, pathname):
        """Checks that pathname is a string and that it fits the proper DSS format"""
        if not isinstance(pathname, str):
            #raise Exception('Please make sure that the DSS path is a string. Input:%s' %pathname)
            pathname = str(pathname)
        pattern =  r"^/([^/]*)/([^/]*)/([^/]*)/([^/]*)/([^/]*)/([^/]*)/$"
        return  bool(re.match(pattern + r'\Z', pathname))
    def getPathname(self):
        """Returns the pathname (uppercase)"""
        return self.pathname.upper()
    def setPathname(self, pathname):
        self.checkValid(pathname)
        self.pathname = pathname
    def getAPart(self):
        return self.pathname.split("/")[1]
    def getBPart(self):
        return self.pathname.split("/")[2]
    def getCPart(self):
        return self.pathname.split("/")[3]
    def getDPart(self):
        return self.pathname.split("/")[4]
    def getEPart(self):
        return self.pathname.split("/")[5]
    def getFPart(self):
        return self.pathname.split("/")[6]
    def setAPart(self, newStr):
        listParts = self.pathname.split("/")
        listParts[1] = newStr
        self.pathname = "/".join(listParts)
    def setBPart(self, newStr):
        listParts = self.pathname.split("/")
        listParts[2] = newStr
        self.pathname = "/".join(listParts)
    def setCPart(self, newStr):
        listParts = self.pathname.split("/")
        listParts[3] = newStr
        self.pathname = "/".join(listParts)
    def setDPart(self, newStr):
        listParts = self.pathname.split("/")
        listParts[4] = newStr
        self.pathname = "/".join(listParts)
    def setEPart(self, newStr):
        listParts = self.pathname.split("/")
        listParts[5] = newStr
        self.pathname = "/".join(listParts)
    def setFPart(self, newStr):
        listParts = self.pathname.split("/")
        listParts[6] = newStr
        self.pathname = "/".join(listParts)
    def toString(self):
        return self.getPathname()

class DSSPathname(DSSPathString):
    """
    hec.heclib.dss.DSSPathname mirror. Everything here should be in DSSPathString
    """

# END OF CLASS DEFINITIONS
################################################################################
# FUNCTION DEFINITIONS

def removeDPart(pathname):
    """Take out the d part (date) from a dss pathname"""
    dssPath = DSSPathString(pathname)
    dssPath.setDPart("")
    return dssPath.getPathname()

def replaceFPart(pathname, newFPart):
    """Replaces the F-part from a dss pathname"""
    dssPath = DSSPathString(pathname)
    dssPath.setFPart(newFPart)
    return dssPath.getPathname()

def getDSSEPartFromIntervalInMinutes(intervalMinutes):
    """
    The E-part of DSS files is in the format "1DAY", "30MIN", etc.
    But it's hard to get this information if all you have is a number that represents
    the timestep interval in minutes. 
    This function does not cover all possible time intervals, just the most
    commonly used ones in modeling
    
    :param int intervalMinutes: The time step in minutes
    :return str: A string representing the corresponding DSS E part (e.g. "1DAY")
                 Returns None if no match was found
    """
    #First, reverse the dictionary
    revDict = {}
    for k in INTERVAL_DICT.keys():
        v = INTERVAL_DICT[k]
        revDict[v] = k
    #Now do the lookup
    if intervalMinutes in revDict:
        dssEPart = revDict[intervalMinutes]
    else:
        return None
    return dssEPart

def getIntervalInMinutesFromDSSEPart(dssEPart):
    """
    The reverse of :func:`getDSSEPartFromIntervalInMinutes`
    
    :param str dssEPart: A string representing the corresponding DSS E part (e.g. "1DAY")
    :return: The time step in minutes
             Returns None if no match was found
    """
    if dssEPart.upper() in INTERVAL_DICT:
        intervalMinutes = INTERVAL_DICT[dssEPart.upper()]
    else:
        return None
    return intervalMinutes

def getPyDateTimeFromDSSTimeMinutes(dssTimeMinutes):
    """
    Returns a python datetime.datetime object corresponding to the input DSS time in minutes
    e.g. input 0, returns 01Jan1900 00:00:00
    e.g. input 1440, returns 02Jan1900 00:00:00
    """
    if isinstance(dssTimeMinutes, float) or isinstance(dssTimeMinutes, int):
        outDT = BASE_DATETIME_DSS + timedelta(minutes=dssTimeMinutes)
    elif isinstance(dssTimeMinutes, datetime.datetime):
        outDT = dssTimeMinutes
    else:
        raise AssertionError("Input %s is unrecognized class: %" %(dssTimeMinutes, dssTimeMinutes.__class__))
    return outDT

def getDSSTimeMinutesFromPyDateTime(pyDateTime):
    """
    Inverse of .getPyDateTimeFromDSSTimeMinutes. 
    Shouldn't really be used much, we should always be operating in the python datetime world,
    not the DSS time world.
    But if we need to put something back to a true java TimeSeriesContainer, we need this
    """
    if hasattr(pyDateTime, "minute"): #a true datetime object
        return int((pyDateTime - BASE_DATETIME_DSS).total_seconds()/60.)
    else: #datetime.date object
        #Assume defined at 2400 hours, not 0000 hours
        return int((pyDateTime - BASE_DATE_DSS).total_seconds()/60. + 1440)

def createTSC(times, vals, unitStr, typeStr, bpart=None, cpart=None, epart=None, fpart=None, pathname=None):
    """
    Create a TimeSeriesContainer object for output, populating all
    necessary fields in the object. 
    
    :param list times: A list of datetime.datetime or datetime.date objects with times to be stored
                       If datetime.date provided, they are converted to datetime and assumed end of day values
    :param list vals: A list of values to be stored
    :param str unitStr: e.g. ft or cfs
    :param str typeStr: e.g. INST-VAL or PER-AVER
    Optional arguments. The user needs to supply either the pathname argument, or else
      all of the bpart, cpart, epart, fpart arguments
    :param str pathname: The dss pathname (e.g. "//MICA/ELEV//1DAY/TEST/")
    :param str bpart: The b-part ("location")
    :param str cpart: The c-part ("parameter")
    :param str epart: The time interval (e.g. "1DAY")
    :param str fpart: The f-part ("version")

    Returns a TimeSeriesContainer object
    """

    tsc = TimeSeriesContainer()
    tsc.numberValues = len(vals)
    tsc.values = vals
    if not isinstance(times, list): 
        #PyDssTools is picky about it being a list, not a numpy Array or pandas DateTimeIndex
        times = list(times)
    if len(times) > 0:
        #if isintance(times[0], datetime.date): #This statement doesn't always work, sometimes is true if the object is datetime.datetime
        if hasattr(times[0], "real"): #integer or float
            times = [getPyDateTimeFromDSSTimeMinutes(t) for t in times]
        if not hasattr(times[0], "minute"): #If only date supplied with no time, assume it is at end of day (2400 hours)
            times = [datetime.datetime(t.year, t.month, t.day, 0, 0) + timedelta(days=1) for t in times]
        startDateTime = times[0]
        tsc.startDateTime = startDateTime.strftime("%d%b%Y %H:%M:00") # "15JUL2019 19:00:00"
    else:
        tsc.startDateTime = BASE_DATETIME_DSS.strftime("%d%b%Y %H:%M:00") #default value
    tsc.times = times
    tsc.units = unitStr
    tsc.type = typeStr
    #Deal with the pathname
    if not pathname is None:
        #pathname supplied, use it
        pathFinal = pathname
    else:
        #no pathname supplied, must have supplied pathname parts
        if cpart is None:
            raise AssertionError("Need to supply either the pathname argument or the pathname parts to create a valid time series container")
        pathFinal = ("//%s/%s//%s/%s/" % (bpart, cpart, epart, fpart))
    tsc.setPathname(pathFinal)
    epart = DSSPathString(pathFinal).getEPart()
    #Figure out the interval (defined in minutes)
    if "IR" in epart.upper(): 
        #irregular data
        tsc.interval = -1
    else:
        #regular data
        tsInterval = getIntervalInMinutesFromDSSEPart(epart)
        tsc.interval = tsInterval
    return tsc

def prepareTSCont(vals, tims, bpart, cpart, epart, fpart, unitStr, typeStr):
    """
    Mirror of the function in cTsUtils.py
    """
    return createTSC(tims, vals, unitStr, typeStr, bpart, cpart, epart, fpart)

# END OF FUNCTION DEFINITIONS
################################################################################
#TESTING/DEBUGGING
if __name__ == "__main__":
    
    dssFile = "testData/CRT_Observed.dss"
    pathIn = "//MICA-POOL/ELEV//1DAY/OBSERVED/"
    pathOut = "//TEMP/ELEV//1DAY/TEMP/"
    startDate = datetime.date(2015, 3, 25)
    endDate = datetime.date(2015, 5, 28)
    #Test
    dssObj = HecDss.open(dssFile)
    df1 = dssObj.readDF(pathIn, startDate, endDate)
    unitStr = df1.units
    typeStr = df1.type
    print(df1)
    '''
    dssTimes = df1.index
    tsc = createTSC(dssTimes, df1[df1.columns[0]].values, unitStr, typeStr, pathname=pathOut)
    dssObj.writeTSC(tsc)
    dssObj.close()
    #Repeat for IR-Century
    dssFile = "testData/CRT_Forecasts_2020L.dss"
    pathIn = "//DCDB/VOLUME[APR-AUG] FCST/01JAN1900/IR-CENTURY/2020-LEVEL/"
    pathOut = "//TEMP/ELEV//IR-CENTURY/TEMP/"
    dssObj.open(dssFile)
    df1 = dssObj.readDF(pathIn, startDate, endDate)
    tscIR = createTSC(df1.index, df1[df1.columns[0]].values, unitStr, typeStr, pathname=pathOut)
    dssObj.writeTSC(tscIR)
    dssObj.close()
    '''
    