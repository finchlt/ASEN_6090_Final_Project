def getRinexData(station, year, doy):
    import gpsToolBox
    import ftplib
    import os
    import scipy.io

    print('Retrieving RINEX for '+station+' for '+doy+'-'+year)
    ftp=ftplib.FTP('garner.ucsd.edu')
    ftp.login()
    
    print('Connected to FTP')
    
    file=station+doy+'0.'+year[2:4]+'o.Z'
    print('Retrieving: '+file)
    
    fid=open(file,'wb')
    filename='/pub/rinex/'+year+'/'+doy+'/'+file
    print('\tLocated at '+filename)
    
    ftp.retrbinary("RETR " + filename, fid.write)
    
    print('Closing FTP connection')
    ftp.close()
    
    os.system('gunzip -f '+file)
    
    file=station+doy+'0.'+year[2:4]+'o'
    data=gpsToolBox.readRinexObs(file)
    
    print('Saving to '+station+doy+'0.mat')
    scipy.io.savemat(station+doy+'0.mat',data)
    print('Done')

def readRinexObsHeader(filename):
    """Reads in the Rinex Header of an Observation File
    Usage: readRinexObsHeader(filename)
    Output:
    dictionary with following fields
      RecvXYZ      -> Approx Position of Recv from Header
      nObservables -> Number of Observables
      observables  -> List of Observables
      year         -> Year of first observation
      month        -> Month of first observation
      day          -> day of first observation
    Assumes Rinex v2.11
    """
    import numpy as np
    
    fid = open(filename,'r')
    result={}
    for line in fid:
        if line.find('RINEX VERSION') != -1: #Some Error checking
            if line.split()[0] != '2.11':
                print('Warning: File Rinex Version isn\'t 2.11')
            if line.split()[1] != 'OBSERVATION':
                print('Warning: Rinex File isn\'t Observation FIle')
        if line.find('APPROX POSITION XYZ') != -1:
            result['recvXYZ']=np.zeros((1,3),float)
            i=0
            for tmp in line.split()[0:3]:
                result['recvXYZ'][0,i]=float(tmp)
                i=i+1            
        if line.find('TYPES OF OBSERV') != -1:
            result['nObservables']=int(line.split()[0])
            result['observables']=line.split()[1:result['nObservables']+1]
        if line.find('TIME OF FIRST OBS') != -1:
            result['year']=int(line.split()[0])
            result['month']=int(line.split()[1])
            result['day']=int(line.split()[2])
        # print comment from header, jic
        if line.find('COMMENT') != -1:
            print('COMMENT: '+line[0:-8])
        if line.find('END OF HEADER') != -1:
            break
        
    fid.close()
    
    return result
    
def readRinexObs(filename):
    """Reads in the Rinex Observation File
    Usage: readRinexObs(filename)
    Output:
    dictionary with following fields
      Header       -> Dictionary with the fields from readRinexObsHeader()
      Observable   -> An observable key for each observable found in file (32xnEpoch)
    Assumes Rinex v2.11
    """
    import numpy as np
    import math
    
    fid = open(filename,'r')
    header={}

    #Read in Header
    while 1:
        # for line in fid
        # cannot use for loop here, since for will call next() which
        # makes the fid object unseekable, which is needed later on.
        line=fid.readline()
        if line.find('RINEX VERSION') != -1: #Some Error checking
            if line.split()[0] != '2.11':
                print('Warning: File Rinex Version isn\'t 2.11')
            if line.split()[1] != 'OBSERVATION':
                print('Warning: Rinex File isn\'t Observation FIle')
        if line.find('APPROX POSITION XYZ') != -1:
            header['recvXYZ']=np.zeros((1,3),float)
            i=0
            for tmp in line.split()[0:3]:
                header['recvXYZ'][0,i]=float(tmp)
                i=i+1            
        if line.find('TYPES OF OBSERV') != -1:
            header['nObservables']=int(line.split()[0])
            header['observables']=line.split()[1:header['nObservables']+1]
        if line.find('TIME OF FIRST OBS') != -1:
            header['year']=int(line.split()[0])
            header['month']=int(line.split()[1])
            header['day']=int(line.split()[2])
        # print comment from header, jic
        if line.find('COMMENT') != -1:
            print('COMMENT: '+line[0:-8])
        if line.find('END OF HEADER') != -1:
            break
    
    eohByte=fid.tell() # End of Header byte
    
    result={}
    result['header']=header

    # read in entire file to find out number of epochs
    # this is stupid, need to find a more efficient way
    # hint: find growable numpy entity
    nEpochs=0
    for line in fid:
        nSV=int(line[30:-1].split('G')[0])
        nEpochs=nEpochs+1
        for i in range(nSV):
            for o in range(math.ceil(header['nObservables']/5)): #read in req no of lines
                line=fid.readline()

    print('Number of Epochs: '+repr(nEpochs))

    # seek back to begining of observations (after header)
    fid.seek(eohByte)

    # setup the observable arrays
    for o in header['observables']:
        result[o]=np.zeros((32,nEpochs))

    result['hour']=np.zeros((1,nEpochs))
    result['min']=np.zeros((1,nEpochs))
    result['sec']=np.zeros((1,nEpochs))
    
    # Start reading in Observations
    epoch=0
    for line in fid:
        result['hour'][0,epoch]=line.split()[3]
        result['min'][0,epoch]=line.split()[4]
        result['sec'][0,epoch]=line.split()[5]
        nSV=int(line[30:-1].split('G')[0])
        prnInView=line[30:-1].split('G')[1:]
        
        #print('epoch: '+repr(epoch)+' | nSV: '+repr(nSV))
        #print(repr(result['hour'][0,epoch])+' : '+repr(result['min'][0,epoch])+' : '+repr(result['sec'][0,epoch]))
        #print('-----')
        
        # read in observables for each PRN
        for i in range(nSV):
            #print('PRN: '+prnInView[i])
            prn=int(prnInView[i])
            line=fid.readline()
            obs=line[:-1].split()
            
            # first line of observations
            obsLine=1
            nObs=0
            nObsLine = math.ceil(header['nObservables']/5)
            
            for o in range(header['nObservables']):
                #print(repr(o)+' : '+header['observables'][o]+' : '+line[16*nObs+1:16*nObs+16])
                try:
                    result[header['observables'][o]][prn-1,epoch]=float(line[16*nObs+1:16*nObs+16])
                except ValueError: # handle blank space in Rinex file
                    result[header['observables'][o]][prn-1,epoch]=0                                
                nObs = nObs + 1
                if (o+1)%5 == 0: #time to read in new line
                    obsLine = obsLine + 1
                    line=fid.readline()
                    obs=line[:-1].split()
                    nObs=0
        #print('----')
        epoch=epoch+1
        
        #break
    fid.close()
    
    return result

def readRinexNav(filename):
    """Reads in the Rinex Nav File
    Usage: readRinexNav(filename)
    Output:
    dictionary with following fields
      Header       -> Dictionary with the fields from readRinexObsHeader()
    Assumes Rinex v2.11    
    """
def rot(angle, axis):
    """ R = rot(angle, axis)
    Input:
      angle in degrees
      axis = {1,2,3}
    Output:
      R -> 3x3 numpy ndarray
    """
    import numpy as np
    R = np.eye(3,dtype=float)
    cang = np.cos(angle*np.pi/180)
    sang = np.sin(angle*np.pi/180)

    if (axis == 1):
        R[1,1]=cang
        R[2,2]=cang
        R[1,2]=sang
        R[2,1]=-sang
    if (axis == 2):
        R[0,0]=cang
        R[2,2]=cang
        R[0,2]=-sang
        R[2,0]=sang
    if (axis == 3):
        R[0,0]=cang
        R[1,1]=cang
        R[1,0]=-sang
        R[0,1]=sang

    return R

def wgsxyz2lla(xyz):
    """ lla = wgsxyz2lla(xyz)
    Input:
      numpy array (1x3) of ECEF coordinates
    Output:
      Dictionary with following elements
        'lat'
        'lon'
        'alt'
    """

    import numpy as np
    
    # check if input is numpy array
    if type(xyz) != np.ndarray:
        print("Input is not a numpy array, blurg!")
        return -1

    if xyz.shape != (1,3):
        print("Input is not of required dimensions (1x3)")
        return -1

    # define some constants
    A_EARTH = 6378137;
    flattening = 1/298.257223563;
    NAV_E2 = (2-flattening)*flattening; # also e^2
    rad2deg = 180/np.pi;

    x = xyz[0,0]
    y = xyz[0,1]
    z = xyz[0,2]
    
    result={}

    if ((x==0.) and (y==0.) and (z==0.)):
        print("ERROR: XYZ at center of earth")
        return 0
    
    if ((x==0.) and (y==0.)):
        result['lon']=0.0
    else:
        result['lon']=np.arctan2(y,x)*rad2deg
        
    # make prelim guesses
    rho = np.sqrt(x*x + y*y)
    templat = np.arctan2(z,rho)
    tempalt = np.sqrt(x*x + y*y + z*z) - A_EARTH
    rhoerror = 1000.0
    zerror = 1000.0

    while ( (np.abs(rhoerror)>1e-6) or (np.abs(zerror)>1e-6) ):
        slat=np.sin(templat)
        clat=np.cos(templat)
        q = 1 - NAV_E2*slat*slat;
        r_n = A_EARTH/np.sqrt(q);
        drdl = r_n*NAV_E2*slat*clat/q; # d(r_n)/d(latitutde)
        
        rhoerror = (r_n + tempalt)*clat - rho;
        zerror   = (r_n*(1 - NAV_E2) + tempalt)*slat - z;
        
        aa = drdl*clat - (r_n + tempalt)*slat;
        bb = clat;
        cc = (1 - NAV_E2)*(drdl*slat + r_n*clat);
        dd = slat;
        
        invdet = 1.0/(aa*dd - bb*cc);
        templat = templat - invdet*(+dd*rhoerror -bb*zerror);
        tempalt = tempalt - invdet*(-cc*rhoerror +aa*zerror);

    result['lat']=templat*rad2deg
    result['alt']=tempalt
    
    return result

def wgslla2xyz(lla):
    """ xyz = wgslla2xyz (lla)
    Input:
      lla -> Dictionary with following elements
              'lat'
              'lon'
              'alt'
    Output:
      xyz -> numpy array (1x3) of ECEF coordinates
    """
    import numpy as np

    if ( (('lat' in lla.keys())==False) and (('lon' in lla.keys())==False) and (('alt' in lla.keys())==False) ):
        print("required keys in lla dictionary doesnt exists!")
        return -1    
        
    if ( (lla['lat']<-90) or (lla['lat']>90) or (lla['lon']<-180) or (lla['lon']>360) ):
        print("Invalid lat/lon values")
        return -1
    
    # some constants
    A_EARTH = 6378137;
    flattening = 1/298.257223563;
    NAV_E2 = (2-flattening)*flattening; # also e^2
    deg2rad = np.pi/180;

    slat = np.sin(lla['lat']*deg2rad);
    clat = np.cos(lla['lat']*deg2rad);
    r_n = A_EARTH/np.sqrt(1 - NAV_E2*slat*slat);

    xyz = np.zeros((1,3))
    xyz[0,0]=(r_n + lla['alt'])*clat*np.cos(lla['lon']*deg2rad)
    xyz[0,1]=(r_n + lla['alt'])*clat*np.sin(lla['lon']*deg2rad)
    xyz[0,2]=(r_n*(1 - NAV_E2) + lla['alt'])*slat

    return xyz

def wgsxyz2enu(xyz, lla):
    """ enu=wgsxyz2enu(xyz, lla)
    Input:
      xyz -> numpy array (1x3) of ECEF coordinates
      lla -> Dictionary with following elements
              'lat'
              'lon'
              'alt'    
    Output:
      numpy array (1x3) with ENU 
    """
    import numpy as np
    
    # check if input is numpy array
    if type(xyz) != np.ndarray:
        print("Input is not a numpy array, blurg!")
        return -1

    if xyz.shape != (1,3):
        print("Input is not of required dimensions (1x3)")
        return -1    

    refxyz = wgslla2xyz(lla)
    if (refxyz.any() == -1):
        return -1

    diffxyz = xyz - refxyz
    
    R1 = rot(90+lla['lon'],3)
    R2 = rot(90-lla['lat'],1)
    R = np.dot(R2,R1)

    enu = np.dot(R,np.transpose(diffxyz))

    return np.transpose(enu)
