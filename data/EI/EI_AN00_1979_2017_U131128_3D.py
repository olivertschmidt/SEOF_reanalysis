#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
 
def retrieve_era20c_mnth():
    """      
       A function to demonstrate how to iterate efficiently over all months,
       for a list of years of the same decade (eg from 2000 to 2009) for an ERA-20C synoptic monthly means request.
       You can extend the number of years to adapt the iteration to your needs.
       You can use the variable 'target' to organise the requested data in files as you wish.
    """
    yearStart = 2000
    yearEnd = 2001
    monthStart = 1
    monthEnd = 12
    requestMonthList = []
    for year in list(range(yearStart, yearEnd + 1)):
        for month in list(range(monthStart, monthEnd + 1)):
            requestMonthList.append('%04d-%02d-01' % (year, month))
    requestMonths = "/".join(requestMonthList)
    target_pl = "era20c_from_%s_to_%s_pl.grb" % (yearStart, yearEnd)
    era20c_mnth_pl_request(requestMonths, target_pl)
 
def era20c_mnth_pl_request(requestMonths, target):
    """      
        An ERA era20c request for analysis, pl data.
        You can change the keywords below to adapt it to your needs.
        (eg add or remove levels, parameters, times etc)
    """
    server.retrieve({
        "class": "e2",
        "stream": "mnth",
        "type": "an",
        "dataset": "era20c",
        "date": requestMonths,
        "expver": "1",
        "levtype": "pl",
        "levelist": "1/5/10/50/100/150/200/250/350/450/550/650/750/800/850/900/950/1000",
        "param": "60.128",
        "target": target,
        "grid" : "1.5/1.5",
        "time": "00"
    })
if __name__ == '__main__':
    retrieve_era20c_mnth()


#!/usr/bin/env python
#from ecmwfapi import ECMWFDataServer
#server = ECMWFDataServer()
#server.retrieve({
#    "class": "ei",
#    "dataset": "interim",
#    "date": "1979-01-01/to/2017-12-31",
#    "expver": "1",
#    "grid": "0.75/0.75",
#    "levelist": "100/150/200/250/350/450/550/650/750/800/850/900/950/1000",
#    "levtype": "pl",
#    "param": "131.128",
#    "step": "0",
#    "stream": "oper",
#    "time": "00:00:00",
#    "type": "an",
#    "format": "netcdf",
#    "target": "EI_AN00_1979_2017_U131128_3D.nc",
#})

# E20C_MONTHLYMEAN00_1900_2010_U131128_3D