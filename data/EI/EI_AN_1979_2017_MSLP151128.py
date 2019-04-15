#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class"     : "ei",
    "dataset"   : "interim",
    "date"      : "1979-01-01/to/2017-12-31",
    "expver"    : "1",
    "grid"      : "0.75/0.75",
    "levtype"   : "sfc",
    "param"     : "151.128",
    "step"      : "0",
    "stream"    : "oper",
    "time"      : "00:00:00/06:00:00/12:00:00/18:00:00",
    "type"      : "an",
    "format"    : "netcdf",
    "target"    : "EI_AN_1979_2017_MSLP151128.nc",
})
