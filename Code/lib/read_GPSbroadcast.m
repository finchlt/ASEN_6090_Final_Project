function [gps_ephem,ionoparams] = read_GPSbroadcast(navfilename)

%==========================================================================
%==========================================================================
% [gps_ephem,ionoparams] = read_GPSbroadcast(navfilename)
%
% Read in an IGS Broadcast Ephemeris for all GPS satellites and construct
%   matrix of all ephemeris values.
%
%
% Author: Ben K. Bradley
% date: 07/19/2009
%
%
% INPUT:             Description                                     Units
%
%  navfilename    - name of IGS broadcast ephemeris file to read in   tring
%
%
% OUTPUT:       
%    
%  gps_ephem      - matrix of gps satellite orbit parameters
%
%     = [prn M0 delta_n e sqrt_a Loa i perigee ra_rate i_rate Cuc Cus Crc...
%            Crs Cic Cis Toe IODE GPS_week Toc Af0 Af1 Af2 0 health]
%  
%  ionoparams     - parameters for the Klobuchar  [A0 A1 A2 A3 B0 B1 B2 B3]
%                    ionospheric model                        
%          
%
% Coupling:
%
%  none
%
% References:
% 
%  [1] Interface Control Document: IS-GPS-200D
%        < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%  [2] RINEX GPS Format, Version 2, (Table A4)
%        < http://www.ngs.noaa.gov/CORS/instructions2/ >
%
%==========================================================================
%==========================================================================


% Open desired IGS Ephemeris File
% =========================================================================
%fid = fopen('brdc1960.09n','r');
if (exist(navfilename,'file') == 2)

    fid = fopen(navfilename,'r');
else
    error(sprintf('Unable to find broadcast file: %s',navfilename), 'ERROR!');
end


% Step through the header of the file and pull out iono parameters
% =========================================================================
headerend   = [];
headeralpha = [];  ALPHA = [];
headerbeta  = [];  BETA  = [];

while (isempty(headerend) == 1)
   tline     = fgetl(fid); 
   headerend = findstr(tline,'END OF HEADER');
   
   headeralpha = findstr(tline,'ION ALPHA');
   if (isempty(headeralpha) == 0)
       
      [A0, remain] = strtok(tline);
      [A1, remain] = strtok(remain);
      [A2, remain] = strtok(remain);
      [A3]         = strtok(remain);
      
      ALPHA = [str2num(A0) str2num(A1) str2num(A2) str2num(A3)];
      
   end
   
   headerbeta = findstr(tline,'ION BETA');
   if (isempty(headerbeta) == 0)
       
      [B0, remain] = strtok(tline);
      [B1, remain] = strtok(remain);
      [B2, remain] = strtok(remain);
      [B3]         = strtok(remain);
      
      BETA = [str2num(B0) str2num(B1) str2num(B2) str2num(B3)];
      
   end
   
end

ionoparams = [ALPHA BETA];





j = 1;

% Enter main loop to read the rest of the ephemeris file
%==========================================================================
%==========================================================================
while 1
    
    % Load next line in ephemeris file
    tline = fgetl(fid);
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(tline), break, end
   
    

        %-----------------------------------------------------------------
        % Read in variables of the FIRST line this satellite's ephemeris
        %-----------------------------------------------------------------
        prn = str2num(tline(1:2));   % prn number of satellite
       
        Af0 = str2num(tline(23:41)); % clock bias
        Af1 = str2num(tline(42:60)); % clock drift
        Af2 = str2num(tline(61:79)); % clock drift rate
        
        
        
        %-----------------------------------------------------------------
        % SECOND LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in second line of satellite ephemeris
        
        IODE = str2num(tline(4:22));    %Issue of Data (Ephemeris)  
        
        Crs  = str2num(tline(23:41));   %Amplitude of the Sine Harmonic Correction 
                                        %  Term to the orbit radius
      
        delta_n= str2num(tline(42:60)); %Mean Motion Difference from Computed Value, rad/s
        
        M0   = str2num(tline(61:79));   %Mean Anomaly at Reference Time, rad
        
        %-----------------------------------------------------------------
        % THIRD LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in third line of satellite ephemeris
        
        Cuc = str2num(tline(4:22));     %Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Argument of Latitude
       
        ecc = str2num(tline(23:41));    %Eccentricity
        
        Cus = str2num(tline(42:60));    %Amplitude of the Sine Harmonic Correction
                                        %  Term to the Argument of Latitude
       
        sqrt_a = str2num(tline(61:79)); %Square Root of the Semi-Major Axis, m^1/2
        
        %-----------------------------------------------------------------
        % FOURTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in fourth line of satellite ephemeris
        
        Toe = str2num(tline(4:22));     %Reference Time Ephemeris (sec into GPS week)
        
        Cic = str2num(tline(23:41));    %Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Angle of Inclination
        
        Loa = str2num(tline(42:60));    %Longitude of Ascending Node of Orbit Plane
                                        %  at Weekly Epoch, rad
      
        Cis = str2num(tline(61:79));    %Amplitude of the Sine Harmonic Correction
                                        %  Term to the Angle of Inclination
         
        %-----------------------------------------------------------------                    
        % FIFTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in fifth line of satellite ephemeris
        
        incl = str2num(tline(4:22));    %Inclination Angle at Reference Time, rad
        
        Crc  = str2num(tline(23:41));   %Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Orbit Radius
     
        perigee = str2num(tline(42:60));%Argument of Perigee, rad
        
        ra_rate = str2num(tline(61:79));%Rate of Change of Right Ascension, rad
        
        %-----------------------------------------------------------------
        % SIXTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in sixth line of satellite ephemeris
        
        i_rate = str2num(tline(4:22));   %Rate of change of inclination angle, rad/s
        
        %str = tline(23:41);             %codes on L2 channel (unecessary)
       
        GPS_week = str2num(tline(42:60));%GPS Week Number (to go with Toe)
        
        %str   = tline(61:79);           %L2 flag
        
        %-----------------------------------------------------------------
        % SEVENTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %includes: SV accuracy, SV health, TGD, IODC
        
        health = str2num(tline(23:41)); %Satellite health (0.00 = usable)
        
        
        %-----------------------------------------------------------------
        % EIGHTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in eighth line of satellite ephemeris
        
        Toc = 0; %Time of clock
        
        
        
        
        gps_ephem(j,:) = [prn M0 delta_n ecc sqrt_a Loa incl perigee ra_rate i_rate Cuc Cus Crc Crs Cic Cis Toe IODE GPS_week Toc Af0 Af1 Af2 0 health];
        
        j = j + 1;   
        
end


fclose(fid);




