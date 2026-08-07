"""
This module calculates an hourly hydrograph from daily average flows 
And instantaneous peak records.

Generally, a user would just call the "createHourlyUsingSpline" function

Ryan Cahill
May 2022

"""

import pandas as pd
import numpy as np
import datetime
from datetime import timedelta

###############################################################################
#CONFIGURATION

################################################################################
# CLASS DEFINITIONS


################################################################################
# FUNCTION DEFINITIONS

def prepDFIndexes(dfDaily, dfPeaks):
    """
    Function to prepare the dataframe objects for the hydrograph shaping process
    """
    #DSS thinks daily average defined at 2400 hours, python thinks it's defined at 0000 hours
    #Reset index to just be the date of daily average data without the time (period-average)
    if not isinstance(dfDaily.index, pd.PeriodIndex):
        if len(dfDaily) > 0:
            dfDaily = dfDaily.set_index(dfDaily.index.to_period("D") - timedelta(days=1))
    #Reset index to just be the date of daily average data without the time (instantaneous)
    #For recent years, the peak comes along with a time, but for older years, just a day is given (2400 hours)
    #Subtract a second to get it in the proper date
    if isinstance(dfPeaks.index, pd.DatetimeIndex):
        if len(dfPeaks) > 0:
            dfPeaks = dfPeaks.set_index((dfPeaks.index - timedelta(seconds=1)).date)
    return dfDaily, dfPeaks
    
def expandByOneDay(df):
    """
    When reading daily data, we often want to go one day before
    
    The beginning/end of the timewindow usually aren't interesting, we just want to make sure we have continuous data.
    Expand the dataset on the front end by a day
    Just copy the first value
    
    Input argument can be a Series or DataFrame
    """
    startDate = df.index[0] + datetime.timedelta(days=-1)
    valueStart = df.iloc[0]
    df.loc[startDate] = valueStart
    df = df.sort_index()
    return df

def create3PtDF(dfDaily):
    '''
    Create dataframe with 3 points per day defined (beginning of day, noon, end)
    Extent is from 0000 on the first day to 2400 on the last day
    '''
    df3Pts = dfDaily.resample("12h").ffill()
    df3Pts = df3Pts.set_index(df3Pts.index)
    #The last point ends at 1200 hours, want to have it go all the way to 2400 hours
    df3Pts.loc[df3Pts.index[-1] + timedelta(hours=12)] = np.nan
    df3Pts['value'] =  np.nan
    return df3Pts

def fillAroundPeak(df3Pts, dfDaily, dfPeaks):
    '''
    For all the instantaneous maximum peak records, fill data for the day of peak flow
    '''
    for peakDate, peakRow in dfPeaks.iterrows():
        peakDate = peakDate
        qPeak = peakRow['value']
        
        #Day 1 is the day before peak, day 2 is the peak day, and day 3 is the day after peak
        if not peakDate in dfDaily.index: continue
        q2 = dfDaily.loc[peakDate]['value']
        if peakDate - timedelta(days=2) not in dfDaily.index: #Going beyond start of data
            q0 = q1 = q2
        else:
            q0 = dfDaily.loc[peakDate - timedelta(days=2)]['value']
            q1 = dfDaily.loc[peakDate - timedelta(days=1)]['value']
        if peakDate + timedelta(days=2) not in dfDaily.index: #Going beyond end of data
            q3 = q4 = q2
        else:
            q3 = dfDaily.loc[peakDate + timedelta(days=1)]['value']
            q4 = dfDaily.loc[peakDate + timedelta(days=2)]['value']
        if pd.isnull(q2): continue #No daily flow, just annual peak. Skip
        '''
        #Diagnostic check to see if daily max value is concurrent with the peak max
        dfDailyWindow = dfDaily.loc[((dfDaily.index > (peakDate - timedelta(days=15)))&(dfDaily.index < (peakDate + timedelta(days=15))))]
        dfDailyMaxDate = dfDailyWindow.idxmax().values[0].to_timestamp()
        dfDailyMaxDate = datetime.date(dfDailyMaxDate.year, dfDailyMaxDate.month, dfDailyMaxDate.day)
        daysDiff = abs((peakDate - dfDailyMaxDate).days)
        if daysDiff>1: print(" Peak date: %s, Daily date: %s" %(peakDate, dfDailyMaxDate))
        '''
        #Get the residual flow (above the minimum value) for the 3 days around peak flow
        qMin = min(q1, q2, q3)
        r1 = q1 - qMin
        r2 = q2 - qMin
        r3 = q3 - qMin
        #Calculate the time of day that we want the instantaneous peak to land
        #0<=t<=1, with t=0 at beginning of day and t=1 at end of day
        if q2 < q1: t = 0
        elif q2 < q3: t = 1
        else: t = max(0,min(1,0.4 + (r3-r1)/(r1+r2+r3)))
        #Round it to the nearest hour (still between 0 and 1)
        t = round(t*24,0)/24.
        tHour = max(1,min(23, 24*t))
        #Now, deal with the data point at beginning of day for the peak day (day 2)
        #We'll call that point A
        #Point B is the end of the day for the peak day
        #Calculate minimum and maximum possible values for qA
        #Assume that on day 1, the slope of the hydrograph would never go negative
        #Case 1 would be to hold the flow constant for the first 12 hours, then pick end of day value such that daily volume is conserved
        #Case 2 would be to hold the flow constant for the last 12 hours, then pick the mid-day value such that daily volume is conserved
        #Do some basic algebra with trapezoidal areas (1st half and 2nd half) to get these equations
        v1 = q1*24 #cfs-hours
        vPeakDay = q2*24 #cfs-hours
        #Get The flow target at the beginning of day 1
        dtDay1Begin = datetime.datetime(peakDate.year, peakDate.month, peakDate.day,0) - timedelta(days=1)
        if dtDay1Begin in df3Pts.index and not pd.isnull(df3Pts.loc[dtDay1Begin]['value']):
            #If the value has already been defined, use it
            q01 = df3Pts.loc[dtDay1Begin]['value']
        else:
            #Use a simple assumption of halfway between day 0 flow and day 1 flow
            q01 = (q0+q1)/2 
        qA_Case1 = 2*(v1 - q01*12)/12. - q01
        qA_Case2 = (v1 - 12*q01/2.)/18.
        qAMax = min(qPeak, max(q1, qA_Case1))
        #The minimum value for qA is that which would cause qB to be exactly equal to qP
        #If qA was lower than this, then qB would have to be higher than qP to balance the daily volume
        qA_Min_For_qB = 2*(vPeakDay - (24-tHour)*qPeak)/tHour - qPeak
        qAMin = max(qA_Case2, qA_Min_For_qB)
        
        #Now deal with point B (end of the day for the peak day)
        #Get the flow target at the end of day 3
        v3 = q3*24 #cfs-hours
        dtDay3End = datetime.datetime(peakDate.year, peakDate.month, peakDate.day,0) + timedelta(days=2)
        if dtDay3End in df3Pts.index and not pd.isnull(df3Pts.loc[dtDay3End]['value']):
            #If the value has already been defined, use it
            q34 = df3Pts.loc[dtDay3End]['value']
        else:
            #Use a simple assumption of halfway between day 0 flow and day 1 flow
            q34 = (q3+q4)/2 
        qB_Case1 = 2*(v3 - q34*12)/12. - q34
        qB_Case2 = (v3 - 12*q34/2.)/18.
        qBMax = min(qPeak, max(q3, qB_Case1))
        #The minimum value for qB is that which would cause qA to be exactly equal to qP
        #If qB was lower than this, then qA would have to be higher than qP to balance the daily volume
        qB_Min_For_qA = 2*(vPeakDay - tHour*qPeak)/(24-tHour) - qPeak
        qBMin = max(qB_Case2, qB_Min_For_qA)
        
        #Drop the original mid-day point for the day of peak
        midDayDateTime = datetime.datetime(peakDate.year,peakDate.month,peakDate.day, 12)
        if midDayDateTime in df3Pts.index:
            df3Pts = df3Pts.drop(midDayDateTime)
        
        #Now, actually set qA and qB
        #If t = 0, we want to use as high of value as possible for qA
        #If t = 1, we want to use as low of a value as possible for qA
        #First, make sure that a solution using a triangular shape is actually possible. 
        #We could run into 
        #The minimum possible volume would occur if the beginning/end of the peak day were as low as possible
        vMin = tHour*(qAMin+qPeak)/2. + (24-tHour)*(qPeak + qBMin)/2.
        #The maximum possible volume would occur if the beginning/end of the peak day were as high as possible
        vMax = tHour*(qAMax+qPeak)/2. + (24-tHour)*(qPeak + qBMax)/2.
        if vPeakDay < vMin or vPeakDay > vMax:
            #Uh oh, we are in a funky situtation where it is impossible to satisfy the daily avg volume using a triangular shape
            #print("Triangular solution failed for peak on %s, using another intermediate point" %peakDate)
            #print("\tqTarget=%0.0f cfs\n\tqMin=%0.0f cfs\n\tqMax=%0.0f cfs\n\tHour=%0.0f" %(q2,vMin/24., vMax/24., tHour))
            #Need to add an intermediate point
            tHrsNewPt = 12
            qA = (1-t)*qAMax + t*qAMin
            qB = t*qBMax + (1-t)*qBMin
            if t > 0.5: #peak occurs late
                t1 = tHrsNewPt
                t2 = tHour - tHrsNewPt
                t3 = 24 - tHour
                qNew = ((vPeakDay - t3*(qPeak + qB)/2. - t1/2*qA - t2/2*qPeak)/(t1/2 + t2/2))
            else: #peak occurs early
                t1 = tHour
                t2 = tHrsNewPt - tHour
                t3 = 24 - tHrsNewPt
                qNew = ((vPeakDay - t1*(qPeak + qA)/2. - t3/2*qB - t2/2*qPeak)/(t3/2 + t2/2))
            dtNewPt = datetime.datetime(peakDate.year, peakDate.month, peakDate.day, tHrsNewPt)
            df3Pts.loc[dtNewPt] = qNew
        else:
            #Triangular shape works just fine
            qA = (1-t)*qAMax + t*qAMin
            #We know qA and qPeak, now we can find qB (flow at the end of the day of peak)
            qB = 2*(vPeakDay - tHour*(qA+qPeak)/2.)/(24-tHour) - qPeak
            #Double-check qB is still reasonable, and adjust if needed
            qB = min(qBMax, max(qBMin, qB))
            #Double-back on qA
            qA = 2*(vPeakDay - (24-tHour)*(qB+qPeak)/2.)/tHour - qPeak
            #Make sure we did it right
            vCheck = tHour*(qA+qPeak)/2. + (24-tHour)*(qPeak+qB)/2.
            qCheck = vCheck/24.
            if round(qCheck,0) != round(q2,0):
                print("Solution failed to converge on daily average flow on day of peak: %s" %peakDate)
        
        #Set the 3 points
        dtPeak = datetime.datetime(peakDate.year,peakDate.month,peakDate.day, int(tHour))
        dtA = datetime.datetime(peakDate.year,peakDate.month,peakDate.day, 0)
        dtB = dtA + timedelta(days=1)
        df3Pts.loc[dtPeak] = qPeak
        df3Pts.loc[dtA] = qA
        df3Pts.loc[dtB] = qB
        df3Pts = df3Pts.sort_index()
    return df3Pts

def fillTheRest(df3Pts, dfDaily):
    """
    After we've filled in the data around the peak, need to deal with the rest of the data
    In general, the value at midnight will be the average of the flow of the 2 days
    The point at noon will be figured such that the daily average volume is conserved
    
    This function is really slow, could be majorly improved. But not worth the effort. 
    """
    
    series3Pts = df3Pts['value']
    #dfMidnight = df3Pts.loc[df3Pts.index.hour == 0]
    i = 0
    for dt, row in dfDaily.iterrows():
        #Get flow for previous, current, next day
        qCurrentDay = row['value']
        if i == 0:
            qPrevDay = qCurrentDay
        else:
            qPrevDay = dfDaily.loc[dt - timedelta(days=1)]['value']
        if i == (len(dfDaily)-1):
            qNextDay = qCurrentDay
        else:
            qNextDay = dfDaily.loc[dt + timedelta(days=1)]['value']
        i += 1
        #If no valid daily data, skip
        if pd.isnull(qCurrentDay): 
            continue
        dayStart = datetime.datetime(dt.year, dt.month, dt.day, 0)  
        dayEnd =  dayStart + timedelta(days=1)
        #Do the beginning of day
        if pd.isnull(series3Pts.loc[dayStart]):
            q0 = (qPrevDay + qCurrentDay)/2.
            series3Pts.loc[dayStart] = q0
        else:
            q0 = series3Pts.loc[dayStart]
        #Do the end of day
        if i==len(dfDaily) or pd.isnull(series3Pts.loc[dayEnd]):
            q1 = (qNextDay + qCurrentDay)/2.
            series3Pts.loc[dayEnd] = q1
        else:
            q1 = series3Pts.loc[dayEnd]
        #Do the middle of the day
        #OK, need to do the noon point
        dtNoon = dayStart + timedelta(hours=12)
        if dtNoon in series3Pts.index:
            if pd.isnull(series3Pts.loc[dtNoon]):
                vTarget = qCurrentDay*24
                qNoon = (vTarget - 6*q0 - 6*q1)/12
                series3Pts.loc[dtNoon] = qNoon
    df3Pts['value'] = series3Pts
    return df3Pts

def fillNoons(df3Pts, dfDaily):
    """
    After we've filled in the data around the peak and the midnight values, 
    need to fill in the noon values. 
    The point at noon will be figured such that the daily average volume is conserved
    
    This function is really slow, could be majorly improved. But not worth the effort. 
    """
    
    #series3PtsNoon = df3Pts.loc[df3Pts.index.hour==12]['value']
    series3Pts = df3Pts['value']
    #dfMidnight = df3Pts.loc[df3Pts.index.hour == 0]
    for dt, row in dfDaily.iterrows():
        #If no valid daily data, skip
        qCurrentDay = row['value']
        if pd.isnull(qCurrentDay): 
            continue
        dayStart = datetime.datetime(dt.year, dt.month, dt.day, 0)  
        dayEnd =  dayStart + timedelta(days=1)
        if dayStart in series3Pts: 
            q0 = series3Pts.loc[dayStart]
        else: 
            q0 = series3Pts.iloc[0]
        if dayEnd in series3Pts: 
            q1 = series3Pts.loc[dayEnd]
        else: 
            q1 = series3Pts.iloc[-1]
        #Do the middle of the day (noon)
        dtNoon = dayStart + timedelta(hours=12)
        if dtNoon in series3Pts.index:
            if pd.isnull(series3Pts.loc[dtNoon]):
                vTarget = qCurrentDay*24
                qNoon = (vTarget - 6*q0 - 6*q1)/12
                series3Pts.loc[dtNoon] = qNoon
    df3Pts['value'] = series3Pts
    return df3Pts

def applyFloorToNoons(df3Pts, minFlow):
    '''
    Sometimes, the 3 point method produces flow lower than a reasonable minimum flow value at noon
    This step is usually intended to be the final step (after estimating noon values and values around the peak)
    It's assumed that the midnight values are already OK
    
    Moves the noon point upward to the minimum value, and then reduces the larger of the two midnight point to compensate
    Also adjusts the noon points at the adjacent days to make sure that daily volume is conserved
    '''
    series3Pts = df3Pts['value']
    for dtNoon,row in df3Pts.loc[((df3Pts.index.hour==12) & (df3Pts["value"]<minFlow))].iterrows():
        #print("Noon value was negative on: %s" %dtNoon)
        qNoon = row['value']
        qDiff = minFlow - qNoon
        #Find whether beginning of day or end of day is higher, and adjust it
        dtDayStart = dtNoon - timedelta(hours=12)
        dtDayEnd = dtNoon + timedelta(hours=12)
        qDayStart = series3Pts[dtDayStart]
        qDayEnd = series3Pts[dtDayEnd]
        if qDayEnd > qDayStart:
            #Adjust the end of day
            qNew = qDayEnd - 2*qDiff
            #Adjust noon of next day
            dtNoonNext = dtNoon + timedelta(days=1)
            if dtNoonNext in series3Pts.index:
                qNewNoon = series3Pts[dtNoonNext] + qDiff
                series3Pts.loc[dtDayEnd] = qNew
                series3Pts.loc[dtNoonNext] = qNewNoon
        else:
            #Adjust the beginning of day
            qNew = qDayStart - 2*qDiff
            #Adjust noon of next day
            dtNoonPrev = dtNoon - timedelta(days=1)
            if dtNoonPrev in series3Pts.index:
                qNewNoon = series3Pts.values[0] + qDiff
                series3Pts.loc[dtDayStart] = qNew
                series3Pts.loc[dtNoonPrev] = qNewNoon
        series3Pts.loc[dtNoon] = minFlow
    df3Pts['value'] = series3Pts
    return df3Pts

def createHourlyUsingAvg(dfDaily, dfPeaks):
    '''
    Outputs a dataframe of hourly data, using daily average data and instantaneous peak as input
    midnight values are the average of the previous and next days. Noon points adjusted to make volumes match
    Create the dataframe that will have 3 points per day (beginning, middle (somewhere), and end)
    
    DEPRECATED
    
    This method is generally not preferred--the spline method is more robust
    '''
    #Create the dataframe that will have 3 points per day (beginning, middle (somewhere), and end)
    df3Pts = create3PtDF(dfDaily)
    df3Pts = fillAroundPeak(df3Pts, dfDaily, dfPeaks)
    df3Pts = fillTheRest(df3Pts, dfDaily)
    dfHourly = df3Pts.resample("1h").mean()
    dfHourly = dfHourly.interpolate(method='linear', limit=24)
    return dfHourly
    
def createHourlyUsingSpline(dfDaily, dfPeaks, minFlow=0):
    '''
    Outputs a dataframe of hourly data, using daily average data and instantaneous peak as input
    Midnight values assigned using cubic spline, which was generated using daily average values at 1200 hours
    Then, erase the values at 1200 hours and replace with whatever is needed to balance the daily average volume.
    '''
    dfDaily, dfPeaks = prepDFIndexes(dfDaily, dfPeaks)
    dfDaily.index = dfDaily.index.to_timestamp()
    dfPeaks.index = pd.to_datetime(dfPeaks.index)
    dfSpline = create3PtDF(dfDaily)
    #Set the daily average value for the day at 1200 hours
    dfDailyAtNoon = dfDaily.set_index(dfDaily.index + timedelta(hours=12))
    dfSpline['value'] = dfDailyAtNoon['value']
    #Cubic spline interpolation to fill in all values
    dfSpline = dfSpline.interpolate(method="cubic", limit=1)
    #Now we've got the values at 0000 and 2400 for the day that will be pretty smooth
    #Make sure we don't go below realistic minimum
    
    dfSpline = dfSpline.clip(lower=minFlow)
    
    #Get rid of the values at 1200 and then we'll recalculate to conserve the daily avg volume
    dfSpline.loc[dfSpline.index.hour==12] = np.nan
    #Fill in first and last value as the first valid value and last valid value
    dfSpline.iat[0, 0] = dfSpline.head().bfill().iloc[0,0]
    dfSpline.iat[-1, 0] = dfSpline.tail().ffill().iloc[-1,0]
    #Fill in data around peak, then the values at 1200
    dfSpline = fillAroundPeak(dfSpline, dfDaily, dfPeaks) 
    dfSpline = fillNoons(dfSpline, dfDaily)
    dfSpline = applyFloorToNoons(dfSpline, minFlow)
    dfSplineHourly = dfSpline.resample("1h").mean()
    dfSplineHourly = dfSplineHourly.interpolate(method='linear', limit=24)
    return dfSplineHourly

#####  Main  ###################################################################
# MAIN CODE
if __name__ == "main" or __name__ == "__main__":
    # Code starts here when called as script.  Otherwise can be loaded as a module
    pass
    
    