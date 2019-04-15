#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1979-01-01/to/2017-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "area": "50/0/-50/180",	
    "levelist": "100/150/200/250/350/450/550/650/750/800/850/900/950/1000",
    "levtype": "pl",
    "param": "60.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00",
    "type": "an",
    "format": "netcdf",
    "target": "EI_AN0012_1979_2017_PV60128_3D.nc",
})
