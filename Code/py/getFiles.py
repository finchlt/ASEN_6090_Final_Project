#!/usr/bin/python

import gpsToolBox
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--station", dest="station", default=None, 
                  help="PBO Station name", metavar="FILE")
parser.add_option("-d", "--doy",dest="doy",
                  help="Day of Year")
parser.add_option("-y", "--year",dest="year",
                  help="Year")


(options, args) = parser.parse_args()
#print(options)

if options.station is not None:
    station=options.station
    year=options.year
    doyI=int(options.doy)
    doy="%03d"%doyI

    gpsToolBox.getRinexData(station,year,doy)

    cmd='mv '+station+doy+'0.mat ./../../Data/mat/'+station+doy+'0_'+year+'.mat'
    print(cmd)
    os.system(cmd)
