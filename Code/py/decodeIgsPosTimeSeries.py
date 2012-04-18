#!/usr/bin/python

import gpsToolBox
import os
import csv
import numpy as np
import scipy.io

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--station", dest="station", default=None, 
                  help="PBO Station name", metavar="FILE")

(options, args) = parser.parse_args()
print(options)

if options.station is not None:
    station=options.station.upper()
    fName=station+'.pbo.igs05.csv'
    cmd='wget -nc ftp://data-out.unavco.org/pub/products/position/'+station+'/'+fName

    print(cmd)
    os.system(cmd)

    # find no of readings in CSV
    ifile = open(fName,'r')
    reader = csv.reader(ifile)
    
    irow = 0
    for row in reader:
        irow=irow+1
        
    ifile.close()
    nr = irow+1
    
    # read in CSV file
    result={}
    result['year']=np.zeros((1,nr),int)
    result['month']=np.zeros((1,nr),int)
    result['day']=np.zeros((1,nr),int)
    result['soln']=np.zeros((3,nr),float) # row1to3 -> E,N,V
    result['sig']=np.zeros((3,nr),float)  # row1to3 -> E,N,V
    
    ifile = open(fName,'r')
    reader = csv.reader(ifile)
    
    irow = 0
    for row in reader:
        irow=irow+1
        
        if irow<10:
            continue
        
        print(irow)

        icol=0
        for col in row:
            if icol == 0:
                date=col.split('-')
                result['year'][0,irow]=int(date[0])
                result['month'][0,irow]=int(date[1])
                result['day'][0,irow]=int(date[2])
            elif icol == 1:
                result['soln'][1,irow]=float(col)
            elif icol == 2:
                result['soln'][0,irow]=float(col)
            elif icol == 3:
                result['soln'][2,irow]=float(col)
            elif icol == 4:
                result['sig'][1,irow]=float(col)
            elif icol == 5:
                result['sig'][0,irow]=float(col)
            elif icol == 6:
                result['sig'][2,irow]=float(col)

            icol=icol+1
        print('%d/%d/%d'%(result['year'][0,irow], result['month'][0,irow], result['day'][0,irow]))
        print('soln E:%10.6f|N:%10.6f|V:%10.6f' % (result['soln'][0,irow], result['soln'][1,irow], result['soln'][2,irow]) )
        print('sig  E:%10.6f|N:%10.6f|V:%10.6f' % (result['sig'][0,irow], result['sig'][1,irow], result['sig'][2,irow]) )
        
    ifile.close()

    scipy.io.savemat(station+'_timeSeries.mat',result)
    cmd='mv '+station+'_timeSeries.mat ./../../Data/'+station+'_timeSeries.mat'
    print(cmd)
    os.system(cmd)
