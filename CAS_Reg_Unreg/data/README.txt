CAS_Reg_Unreg / data
====================

This folder holds only the small tabular inputs that belong to this task.
There is no observed-data DSS here on purpose.

ObsData_RegUnreg.dss was removed on 7 Aug 2026 because it duplicated records
already held in the CAS_Unreg_FF study. The canonical observed-data store for
BOTH tasks is now:

    ../../CAS_Unreg_FF/data/obsData.dss

That file is DSS version 6. Scripts reach it with a relative path from their own
location, so nothing needs a machine-specific path.

Local files in this folder
--------------------------
  Events_Detailed.csv        Analyst-defined event start/end dates
  MOS_Inflow_ESRD.csv        Mossyrock inflow ESRD special-curve table
  MOS_stor_rating_curve.csv  Mossyrock elevation-storage rating curve
  MOS_96.csv                 Mossyrock CDB inflow source table


Records this task uses from ../../CAS_Unreg_FF/data/obsData.dss
---------------------------------------------------------------
R = read      W = written back into the shared store

#Reservoir_Inflow_And_Basin_Peaks.py
  R  //MOS/FLOW-IN//1HOUR/CWMS/
  R  //MOS/FLOW-OUT//1HOUR/CWMS/
  R  //MOS/STOR//1HOUR/CWMS/
  R  //MAY/FLOW-IN//1HOUR/CWMS/
  R  //MAY/FLOW-OUT_PEAKCLEAN_2009_2026//1HOUR/CWMS/
  R  //MAY/STOR_PEAKCLEAN_2009_2026//1HOUR/CWMS/
  W  //<res>/FLOW-IN-CALC-RAW//1HOUR/CWMS/
  W  //<res>/FLOW-IN-CALC-RAW-PEAKS//1HOUR/CWMS/
  W  /LOCAL/CASTLE_ROCK_LOCAL/FLOW//1HOUR/DERIVED/

#Inflow_Volume_Correction.py
  R  //MAY/FLOW-LOCAL//1HOUR/CWMS/
  R  //<res>/FLOW-IN-CALC-RAW//1HOUR/CWMS/
  W  //<res>/FLOW-IN-CALC-CLEANED-VOLCOR//1HOUR/CWMS/
  W  //MAY/FLOW-LOCAL-SHAPED//1HOUR/CWMS/

#Create_Ensembles.py
  R  //MOS/ELEV//1HOUR/CWMS/
  R  //MOS/FLOW-IN-CALC-CLEANED-VOLCOR//1HOUR/CWMS/
  R  //MOS/FLOW-OUT//1HOUR/CWMS/
  R  //MAY/ELEV//1HOUR/CWMS/
  R  //MAY/FLOW-LOCAL-SHAPED//1HOUR/CWMS/
  R  //MAY/FLOW-OUT_PEAKCLEAN_2009_2026//1HOUR/CWMS/
     (writes ensemble members out to the ResSim watershed shared/ensemble.dss,
      not back into obsData.dss)

#ExtractResSimEnsembleResults.py
  R  //MOS/ELEV//1HOUR/CWMS/
  R  //MOS/ELEV-RULECURVE//1HOUR/CENWP-CALC/
  R  //MOS/FLOW-OUT//1HOUR/CWMS/
  R  //MAY/ELEV//1HOUR/CWMS/
  R  //MAY/FLOW-IN-CALC-RAW-PEAKS//1HOUR/CWMS/
  R  //MAY/FLOW-OUT_PEAKCLEAN_2009_2026//1HOUR/CWMS/
     (writes to ../data/results.dss)

src/ResSim_Inflows/Cowlitz_Inflows_for_ResSim.py
  R  //MOS/ELEV//1HOUR/CWMS/
  R  //MOS/FLOW-IN//1HOUR/CWMS/
  R  //MOS/FLOW-OUT//1HOUR/CWMS/
  R  /COWLITZ RIVER AT RANDLE/14231000/FLOW//1DAY/USGS/
  R  /COWLITZ RIVER AT RANDLE/14231000/FLOW//1HOUR/USGS/
  R  /COWLITZ RIVER AT RANDLE/14231000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/
  R  /COWLITZ RIVER AT MOSSYROCK/14235000/FLOW//1DAY/USGS/
  R  /COWLITZ RIVER AT MOSSYROCK/14235000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/
  R  /COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW//1DAY/USGS/
  R  /COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/
  R  //MAYFIELD OUTFLOW/FLOW//1HOUR/USGS 14238000/
  R  /MOS6A/MOS/FLOW-IN//1DAY/2020_Modified_Flow_Report/
     (writes ../../output/ResSimInflows.dss:
      //MOSSYROCK/FLOW-IN//1HOUR/FOR_RESSIM/
      //CASTLE ROCK/FLOW-LOCAL//1HOUR/FOR_RESSIM/)

#MOS_Special_Release_MinFloodPool.py
  PENDING -- still points at the removed ObsData_RegUnreg.dss. It will be
  repointed once the ResSim inputs are final, because the volume-corrected
  inflows and the shaped Mayfield local it uses are being replaced by what is
  in ResSimInflows.dss. Records it needs today:
     //MOS/ELEV-RULECURVE//1Hour/CENWP-CALC/
     //MOS/ELEV-USGS//1Day/USGS/
     //MOS/FLOW-IN-CALC-CLEANED-VOLCOR//1Hour/CWMS/
     //MAY/FLOW-LOCAL-SHAPED//1Hour/CWMS/
     /COWLITZ RIVER BELOW MAYFIELD DAM, WA/14238000/FLOW//1Hour/USGS/

Critical_Duration_Correlation.py reads no DSS. It takes CSVs from the other
task: CAS_Unreg_FF/output/unreg_durations_massbalance.csv,
CAS_Unreg_FF/output/wy_peak_records.csv, and
CAS_Unreg_FF/data/CastleRock_USGS_peaks.csv.


Caution
-------
Scripts marked W write into the shared store that CAS_Unreg_FF also depends on.
Back up obsData.dss before running them, and do not run them casually.
