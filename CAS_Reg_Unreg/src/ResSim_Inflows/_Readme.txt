This folder contains code to create hourly inflows for Mossyrock and locals at Castle Rock for use in ResSim modeling.
These time series are not precise, usually developed by shaping daily average data to hourly.
But they can be used to run ResSim and show how different peak flows at Castle Rock would be under different operational strategies.

The observed data that this script relies on is not in this folder to avoid duplication.
It is here: "<repo root>/CAS_Unreg_FF/data/obsData.dss"
(the script now resolves this itself from its own location, so it can be run from anywhere)

I did add records to that file, notably the 2020 modified flows and daily/hourly/peak records at USGS locations on the Cowlitz River. 
As well as MOS FLow-IN and Flow-OUT from CWMS (had been dropped from this file, but I need it)

Main script is "Cowlitz_Inflows_for_ResSim.py"
The shared modules it uses (utilsDSS, HydrologicRouting, HourlyHydrographShaping)
now live in "<repo root>/Modules" so they are not duplicated per project.

The output is in "../../output/ResSimInflows.dss"

Ryan Cahill
6Aug2026