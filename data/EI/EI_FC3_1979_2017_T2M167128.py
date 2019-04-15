#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
    
server = ECMWFDataServer()
    
server.retrieve({
    'class'     : "ei",
    'dataset'   : "interim",
    'type'      : "fc",
    'stream'    : "oper",
    'expver'    : "1",
    'levtype'   : "sfc",
    'param'     : "167.128",
    'step'      : "3",
    'grid'      : "0.75/0.75",
    'time'      : "00:00:00/12:00:00",
    'date'      : "1979-01-01/to/2017-12-31",
    'format'    : "netcdf",
    'target'    : "EI_1979_2017_T167128.nc"
})
