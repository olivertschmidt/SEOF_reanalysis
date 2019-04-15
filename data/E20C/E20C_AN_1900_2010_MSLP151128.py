#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class"     : "e2",
    "dataset"   : "era20c",
    "date"      : "1900-01-01/to/2010-12-31",
    "expver"    : "1",
    "levtype"   : "sfc",
    "grid"      : "0.75/0.75",
    "param"     : "151.128",
    "stream"    : "oper",
    "time"      : "00:00:00/12:00:00",
    "type"      : "an",
    "format"    : "netcdf",
    "target"    : "E20C_AN_1900_2010_MSLP151128.nc",
})
