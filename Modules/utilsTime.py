"""
Module with utility functions dealing with dates and times
"""
################################################################################
# IMPORTS
import datetime
from datetime import timedelta
import re
################################################################################
# STATIC INPUT

#All lists of 14-period information in this script should be in this order
#The 1st value (index 0) is AG1 and the last value (index 13) is JUL
#AG2 is called AUG, and AP2 is called APR.
PERIOD_NAMES = ["AG1","AUG","SEP","OCT","NOV","DEC","JAN","FEB","MAR","AP1","APR","MAY","JUN","JUL"]
DATE_STR_14PER = ["31JAN","28FEB","31MAR","15APR","30APR","31MAY","30JUN","31JUL","15AUG","31AUG","30SEP","31OCT","30NOV","31DEC"] #Feb will be overwritten during leap years
DAYS_IN_14PERIOD = [15, 16, 30, 31, 30, 31, 31, 28, 31, 15, 15, 31, 30, 31] #Feb will be overwritten during leap years
daysIn14PeriodDict = dict(zip(PERIOD_NAMES, DAYS_IN_14PERIOD))

MONTHS_DICT = {"JAN": 1,"FEB": 2,"MAR": 3,"APR": 4,"MAY": 5,"JUN": 6,"JUL": 7,"AUG": 8,"SEP": 9,"OCT": 10,"NOV": 11,"DEC": 12}

MIDNIGHT_TIME = datetime.time(0,0)
################################################################################
# CLASS DEFINITIONS

# END OF CLASS DEFINITIONS
################################################################################
# FUNCTION DEFINITIONS
def isLeap(year):
    """Returns true if a leap year, false if not"""
    remainder = year % 4
    if remainder == 0:
        return True
    else:
        return False

def getPowerYear(periodName, calYear):
    """Returns the power year for a given dateStr and calendar year.
    The Power Year starts on August 1st. Sep 30, 1929 is power year 1930.
    e.g. if periodName="DEC" and calYear=1930, output is 1931"""
    powYear = calYear
    if periodName in ["AUG","AU1","AU2", "AG1","AG2","SEP","OCT","NOV","DEC"]:
        powYear = calYear + 1
    elif not periodName in PERIOD_NAMES:
        raise AssertionError("Invalid Period Name: %s" %periodName)
    return powYear

def get14PerIdx(periodName):
    """Input period name (e.g. "SEP"), get the index (e.g. 2), with AG1 being 0"""
    idx = PERIOD_NAMES.index(periodName)
    return idx

def get14PeriodsAsDates(year):
    """
    Returns a list of 14-period dates (datetime.date) for a given calendar year
    Accounts for leap years
    """
    periods = list(DATE_STR_14PER)
    if isLeap(year):
        periods[1] = "29FEB"
    else:
        periods[1] = "28FEB"
    periodDates = [datetime.datetime.strptime(p+str(year)+" 0001", "%d%b%Y %H%M") for p in periods]
    periodDates = [datetime.date(d.year, d.month, d.day) for d in periodDates]
    return periodDates

def getDDMMM(dateInput):
    """
    Input a datetime.date or datetime.datetime object, returns in format DDMMM, e.g. "05MAR"
    """
    return dateInput.strftime("%d%b").upper()

def toPyDatetime(dateInput, dateOnly=False):
    '''
    Converts an input into a Python datetime.datetime object.
    If needed, the inputs will be assumed to be at 2400 hours
    Inputs can be:
        string("YYYY-MM-DD")
        string("YYYY-Mon-DD")
        string("DDMMMYYYY HHMM")
        datetime.date (assumed at 2400 hours of the input)
        datetime.datetime (will do nothing)
        
    So if a datetime.date is supplied for January 10th, 
        return value will be a datetime at 0000 hours on January 11th
    '''
    if isinstance(dateInput, datetime.datetime):
        dateOutput = dateInput
    elif isinstance(dateInput, str):
        dateInput = dateInput.upper()
        date_string_pattern = r"^\d{4}-\d{1,2}-\d{1,2}$"
        date_string_pattern2 = r"^\d{4}-[A-Z]{3}-\d{1,2}$"
        date_string_pattern3 = r"^\d{1,2}[A-Z]{3}\d{4} \d{2}\d{2}$"
        string_validity = bool(re.match(date_string_pattern + r'\Z', dateInput))
        string_validity2 = bool(re.match(date_string_pattern2 + r'\Z', dateInput))
        string_validity3 = bool(re.match(date_string_pattern3 + r'\Z', dateInput))
        if string_validity or string_validity2:
            datestring_split = dateInput.split('-')
            if string_validity2:
                datestring_split[1] = MONTHS_DICT[datestring_split[1]]
            else:
                if int(datestring_split[1]) > 12 and int(datestring_split[1]) < 1:
                    raise Exception('Month is out of range. Please make sure it is 1-12.')
            if int(datestring_split[2]) > 31 and int(datestring_split[2]) < 1:
                 raise Exception('Date is out of range. Please make sure that it is 1-31.')
            datestring_dt = datetime.datetime(int(datestring_split[0]), int(datestring_split[1]), int(datestring_split[2])) + timedelta(days=1)
            dateOutput =  datestring_dt
        elif string_validity3:
            if int(dateInput[10:12]) == 24:
                hour_dt = 0
                add_day = True
            else:
                hour_dt = int(dateInput[10:12])
                add_day = False
            date_split = [int(dateInput[5:9]), MONTHS_DICT[dateInput[2:5]], int(dateInput[:2]), hour_dt, int(dateInput[12:])] #Splitting the string to year, month, day, hour, minute   
            datestring_dt = datetime.datetime(date_split[0], date_split[1], date_split[2], date_split[3], date_split[4])
            if add_day:
                datestring_dt = datestring_dt + timedelta(days=1) #So the value will be 0000 at the next day
            dateOutput = datestring_dt
        else: 
            raise Exception('If the date input is a string, the following formats are accepted: \nYYYY-MM-DD\nYYYY-Mon-DD\nDDMonYYYY HHMM\nYou have input: '+ dateInput)
    elif isinstance(dateInput, datetime.date):
        dateOutput = datetime.datetime(dateInput.year, dateInput.month, dateInput.day) + timedelta(days=1)
    else:
        raise AssertionError("Input not supported: %s, class: %s" %(dateInput, dateInput.__class__))
    
    if dateOnly:
        date_only_output = datetime.date(dateOutput.year, dateOutput.month, dateOutput.day)
        if dateOutput.hour == 0 and dateOutput.minute == 0:
            date_only_output -= timedelta(days=1)
        return date_only_output
    return dateOutput

def toPyDate(dateInput):
    '''
    Converts an input into a Python datetime.date object, assumed at 2400 hours to be similar to DSS.
    Inputs can be:
        string("YYYY-MM-DD")
        datetime.date (will do nothing)
        datetime.datetime (if defined at 0000 hours on 05Jan, will return 04Jan (assumed 2400)
        HecTime
    '''
    if isinstance(dateInput, datetime.date):
        if hasattr(dateInput, "to_pydatetime") or hasattr(dateInput, "minute"):
            #This isn't a pure python datetime.date, 
            #but a pandas._libs.tslibs.timestamps.Timestamp or datetime.datetime object
            dateOutput = datetimeToDate(dateInput)
        else:
            #A pure python datetime.date object
            dateOutput = dateInput
    elif isinstance(dateInput, str):
        dateParts = dateInput.split("-")
        if len(dateParts) != 3: raise AssertionError("Invalid input provided: %s" %dateInput)
        dateOutput = datetime.date(int(dateParts[0]), int(dateParts[1]), int(dateParts[2]))
    elif isinstance(dateInput, datetime.datetime):
        dateOutput = datetimeToDate(dateInput)
    elif hasattr(dateInput, "addDays"):
        #Java HecTime object, defined at 2400 hours
        dateOutput = datetime.date(dateInput.year(), dateInput.month(), dateInput.day())
    else:
        raise AssertionError("Input not supported: %s, class: %s" %(dateInput, dateInput.__class__))
    return dateOutput

def datetimeToDate(dtIn):
    '''
    Converts a python datetime.datetime object to a python datetime.date object
    if defined at 0000 hours on 05Jan, will return 04Jan (assumed 2400)
    if defined at 0001 hours on 05Jan, will return 05Jan
    '''
    dateOutput = datetime.date(dtIn.year, dtIn.month, dtIn.day)
    if dtIn.time() == MIDNIGHT_TIME:
        #If defined at 05Jan at 0000, we want to return a date at 04Jan (assumed 2400 hours)
        dateOutput = dateOutput - timedelta(days=1)
    return dateOutput

# END OF FUNCTION DEFINITIONS
################################################################################