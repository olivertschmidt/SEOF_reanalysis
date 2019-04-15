#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1979-01-01/to/2017-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levelist": "10",
    "levtype": "pl",
    "param": "131.128/132.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/12:00:00",
    "type": "an",
    "format": "netcdf",
    "target": "EI_AN0012_1979_2017_UV131128_10hPa.nc",
})
