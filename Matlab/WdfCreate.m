% WdfCreate creates a valid WDF file with given spectral data.
%
% Many different properties need to be supplied to make a valid WDF file.
% If it is desired to store x, y spectral data from MATLAB in a wdf file
% this script will automate the process.
%
% wdf = WdfCreate(xlist,spectra,filename) creates a new WdfCreate object
% with all the necessary properties to create a valid wdf file. The xlist
% and spectral data are necessary properties and must be supplied. The
% filename is an optional property that automatically saves the WDF file,
% otherwise properties may be adjusted first and the file saved at a later
% time.
%
% wdf.WriteWdf(filename) is a method that saves a new WDF file to the given
% filename with the properties of the given WdfCreate object (in this case
% 'wdf').
%
% Accessible Properties
% =====================
%
% After a WdfCreate object is created there are a number of accessible
% properties that can be changed.
%
% Metadata
%
% wdf.TimeStart - start time of the measurement in MATLAB datetime
% wdf.TimeEnd - end time of the measurement in MATLAB datetime
% wdf.SpectralUnits - unit type of the supplied spectral data (see 
% WiREDataUnit)
% wdf.LaserWaveNum - wavelength of the laser used to take spectral data
% wdf.User - user's name
% wdf.Title - measurement's title
%
% Data properties
%
% wdf.Spectra - array of all the spectral data points supplied as a row 
% vector in the case of a single spectra or a 2d matrix in the case of a
% time series style measurement (points in spectra across rows, each
% seperate spectra across columns)
%
% Xlist properties
%
% wdf.XList - array of all the x axis data points supplied as a 1D column
% or row vector
% wdf.XListDataType - data type of the supplied x axis data (see 
% WiREDataType)
% XListUnitType - unit type of the supplied x axis data (see 
% WiREunitType)
%
% ylist properties
%
% Ylist - ylist data
% YListDataType - data type of the ylist data (see WiREDataType)
% YListUnitType - unit type of the ylist data (see WiREUnitType)

% Copyright 2022 Renishaw plc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

classdef WdfCreate < handle
    %% Private WDF properties
    properties (Access = protected)
        % file handle
        Handle = -1;
        % metadata
        Signature = 'WDF1';
        Version = 1;
        Size = 512;
        Flags = 0;
        Uuid = [0,0,0,0];
        NTracks = 0;
        Status = 0;
        NPoints = 0;
        NSpectra = 0;
        NCollected = 0;
        NAccum = 1;
        YListCount = 1;
        XListCount = 0;
        OriginCount = 3;
        AppName = 'WiRE'
        AppVersion = [5,5,0,0];
        ScanType = WiREScanBasicType.Unspecified;
        MeasurementType = WiREMeasurementType.Single;
        % origin properties
        OriginDataType = [WiREDataType.Time,WiREDataType.Flags, ...
            WiREDataType.Checksum];
        OriginUnitType = [WiREDataUnit.FileTime,WiREDataUnit.Arbitrary, ...
            WiREDataUnit.Arbitrary];
        OriginName = {'Time','Flags','Checksum'}
        OriginValues = 0;
        OriginPrimaryDataType = [0,0,0];
        
    end
    
    %% Public WDF properties
    properties (GetAccess = public, SetAccess = public)
        % metadata
        TimeStart = datetime('now');
        TimeEnd = datetime('now');
        SpectralUnits = WiREDataUnit.Counts;
        LaserWaveNum = 0;
        User = getenv('USERNAME');
        Title = 'MATLAB generated file';
        % data properties
        Spectra = 0
        % xlist properties
        XList = 0;
        XListDataType = WiREDataType.Spectral;
        XListUnitType = WiREDataUnit.RamanShift;
        % ylist properties
        Ylist = 1;
        YListDataType = WiREDataType.SpatialY;
        YListUnitType = WiREDataUnit.Pixels;
        
    end
    
    methods
        % Constructor: creates a WdfCreate object storing all of the
        % necessary properties to create a valid wdf file. A single
        % spectrum or time series can be passed in then if a filename is
        % provided the wdf file will be created.
        function this = WdfCreate(xList,spectra,varargin)
            if nargin > 3
                throw(MException('Renishaw:SPD:WiRE:InputArguments', ...
                    'WdfCreate can accept a maximum of 3 inputs'));
            end
            
            this.XList = xList;
            this.Spectra = spectra;
            
            if isempty(this.Spectra) || length(size(this.Spectra)) > 2
                throw(MException('Renishaw:SPD:WiRE:InputArguments', ...
                    'Spectra must be a 1D vector for single scans and a 2D matrix for time series'));
            end
            if isempty(this.XList) || length(size(this.XList)) > 2 || ...
                    (size(this.XList,1) > 1 && size(this.XList,2) > 1)
                throw(MException('Renishaw:SPD:WiRE:InputArguments', ...
                    'X list must be a 1D vector'));
            elseif size(this.XList,1) > 1
                this.XList = this.XList';
            end
            
            this.NPoints = size(this.Spectra,2);
            this.NSpectra = size(this.Spectra,1);
            this.NCollected = size(this.Spectra,1);
            this.XListCount = length(this.XList);
            
            if this.XListCount ~= this.NPoints
                throw(MException('Renishaw:SPD:WiRE:InputArguments', ...
                    'X list points and spectra points must be the same length'));
            end
            
            if this.XListCount < 2
                throw(MException('Renishaw:SPD:WiRE:InputArguments', ...
                    'At least 2 data points are required to make a functional wdf file'));
            end
            if this.NSpectra > 1
                this.MeasurementType = WiREMeasurementType.Series;
                this.OriginPrimaryDataType = [1,0,0];
            end
            
            this.OriginValues = [1:this.NSpectra; ...
                zeros(1,this.NSpectra);zeros(1,this.NSpectra)];
            
            if nargin ==3
                this.WriteWdf(varargin{1})
            end
        end
    end
    
    %% Protected helper methods
    methods (Access = protected)
        % Converts MATLAB datetime into Windows file time (number of 100
        % nanosecond intervals since 1/1/1601)
        function wftTime = WftConvert(this,time)
            windowsFileTimeStart = datetime(1601,1,1);
            wftTime = seconds(time-windowsFileTimeStart).*10^7;
            
        end
    end
    
    %% Public methods
    methods (Access = public)
        % writes a valid wdf file to the given filepath using a WdfCreate's
        % properties
        function WriteWdf(this,fileName)
            this.Handle = fopen(fileName,'wb');
            
            %%% WRITE WDF HEADER %%%
            %   0; uint32 signature; Magic number to check that this is a WDF file (WDF_BLOCKID_FILE)
            fwrite(this.Handle,this.Signature,'char*1','l');
            %   4; uint32 version; The version of this specification used by this file.
            fwrite(this.Handle,this.Version,'uint32','l');
            %   8; uint64 size; The size of this block (512 bytes)
            fwrite(this.Handle,this.Size,'uint64','l');
            %  16; uint64 flags; Flags from the WdfFlags enumeration
            fwrite(this.Handle,this.Flags,'uint64','l');
            %  24; uint32 uuid[4]; a file unique identifier - never changed once allocated
            fwrite(this.Handle,this.Uuid,'uint32','l');
            %  40; uint64 unused0;
            fwrite(this.Handle,0,'uint64','l');
            %  48; uint32 unused1;
            fwrite(this.Handle,0,'uint32','l');
            %  52; uint32 ntracks; if WdfXYXY flag is set - contains the number of tracks used
            fwrite(this.Handle,this.NTracks,'uint32','l');
            %  56; uint32 status; file status word (error code)
            fwrite(this.Handle,this.Status,'uint32','l');
            %  60; uint32 npoints; number of points per spectrum
            fwrite(this.Handle,this.NPoints,'uint32','l');
            %  64; uint64 nspectra; number of actual spectra (capacity)
            fwrite(this.Handle,this.NSpectra,'uint64','l');
            %  72; uint64 ncollected; number of spectra written into the file (count)
            fwrite(this.Handle,this.NCollected,'uint64','l');
            %  80; uint32 naccum; number of accumulations per spectrum
            fwrite(this.Handle,this.NAccum,'uint32','l');
            %  84; uint32 ylistcount; number of elements in the y-list (>1 for image)
            fwrite(this.Handle,this.YListCount,'uint32','l');
            %  88; uint32 xlistcount; number of elements for the x-list
            fwrite(this.Handle,this.XListCount,'uint32','l');
            %  92; uint32 origincount; number of data origin lists
            fwrite(this.Handle,this.OriginCount,'uint32','l');
            %  96; utf8   appname[24]; application name (utf-8 encoded)
            fwrite(this.Handle,this.AppName,'char*1','l');
            for n = 1:24-length(this.AppName)
                fwrite(this.Handle,0,'uint8','l');
            end
            % 120; uint16 appversion[4]; application version (major,minor,patch,build)
            fwrite(this.Handle,this.AppVersion,'uint16','l');
            % 128; uint32 scantype; scan type - WdfScanType enum
            fwrite(this.Handle,uint32(this.ScanType),'uint32','l');
            % 132; uint32 type; measurement type - WdfType enum
            fwrite(this.Handle,uint32(this.MeasurementType),'uint32','l');
            % 136; uint64 time_start; collection start time as FILETIME
            fwrite(this.Handle,this.WftConvert(this.TimeStart),'uint64','l');
            % 144; uint64 time_end; collection end time as FILETIME
            fwrite(this.Handle,this.WftConvert(this.TimeEnd),'uint64','l');
            % 152; uint32 units; spectral data units (one of WdfDataUnits)
            fwrite(this.Handle,uint32(this.SpectralUnits),'uint32','l');
            % 156; float  laserwavenum; laser wavenumber
            fwrite(this.Handle,this.LaserWaveNum,'single','l');
            % 160; uint64 spare[6];
            for n = 1:6
                fwrite(this.Handle,0,'uint64','l');
            end
            % 208; utf8   user[32]; utf-8 encoded user name
            fwrite(this.Handle,this.User,'char*1','l');
            for n = 1:32-length(this.User)
                fwrite(this.Handle,0,'uint8','l');
            end
            % 240; utf8   title[160]; utf-8 encoded title
            fwrite(this.Handle,this.Title,'char*1','l');
            for n = 1:160-length(this.Title)
                fwrite(this.Handle,0,'uint8','l');
            end
            % 400; uint64 padding[6]; padded to 512 bytes
            for n = 1:6
                fwrite(this.Handle,0,'uint64','l');
            end
            % 448; uint64 free[4]; available for third party use
            for n = 1:4
                fwrite(this.Handle,0,'uint64','l');
            end
            % 480; uint64 reserved[4]; reserved for internal use by WiRE
            for n = 1:4
                fwrite(this.Handle,0,'uint64','l');
            end
            
            %%% WRITE DATA BLOCK %%%
            % Block Header
            %   0; uint32 type_id
            fwrite(this.Handle,'DATA','char*1','l');
            %   4; uint32 unique_id
            fwrite(this.Handle,0,'uint32','l');
            %   8; uint64 block_length
            fwrite(this.Handle,16+(4.*this.NSpectra.*this.NPoints),'uint64','l');
            
            % Block Content
            %  16; float  spectra[nspectra]
            fwrite(this.Handle,this.Spectra','single','l');
            
            %%% WRITE XLST BLOCK %%%
            % Block Header
            % 0; uint32 type_id
            fwrite(this.Handle,'XLST','char*1','l');
            %   4; uint32 unique_id
            fwrite(this.Handle,0,'uint32','l');
            %   8; uint64 block_length
            fwrite(this.Handle,24+(4.*length(this.XList)),'uint64','l');
            
            % Block Content
            %  16; uint32 xlist_data_type
            fwrite(this.Handle,uint32(this.XListDataType),'uint32','l');
            %  20; uint32 xlist_unit_type
            fwrite(this.Handle,uint32(this.XListUnitType),'uint32','l');
            %  24; float  xlist[nxlist]
            fwrite(this.Handle,this.XList,'single','l');
            
            %%% WRITE YLST BLOCK %%%
            % Block Header
            % 0; uint32 type_id
            fwrite(this.Handle,'YLST','char*1','l');
            %   4; uint32 unique_id
            fwrite(this.Handle,0,'uint32','l');
            %   8; uint64 block_length
            fwrite(this.Handle,24+(4.*this.YListCount),'uint64','l');
            
            % Block Content
            %  16; uint32 ylist_data_type
            fwrite(this.Handle,uint32(this.YListDataType),'uint32','l');
            %  20; uint32 ylist_unit_type
            fwrite(this.Handle,uint32(this.YListUnitType),'uint32','l');
            %  24; float  ylist
            fwrite(this.Handle,this.Ylist,'single','l');
            
            %%% WRITE ORGN BLOCK %%%
            % Block Header
            %   0; uint32 type_id
            fwrite(this.Handle,'ORGN','char*1','l');
            %   4; uint32 unique_id
            fwrite(this.Handle,0,'uint32','l');
            %   8; uint64 block_length
            fwrite(this.Handle,20+this.OriginCount.*(24+8.*this.NSpectra),'uint64','l');
            
            % Block Content
            %  16; uint32 number_origins
            fwrite(this.Handle,this.OriginCount,'uint32','l');
            for n = 1:this.OriginCount
                if this.OriginPrimaryDataType(n)
                    %  20; uint32 WdfDataType enum
                    fwrite(this.Handle,uint32(this.OriginDataType(n))+uint32(2^31),'uint32','l');
                else
                    fwrite(this.Handle,uint32(this.OriginDataType(n)),'uint32','l');
                end
                %  24; uint32 WdfUnitType enum
                fwrite(this.Handle,uint32(this.OriginUnitType(n)),'uint32','l');
                %  28; utf8   origin_name[16]
                fwrite(this.Handle,char(this.OriginName(n)),'char*1','l');
                for nn = 1:16-length(char(this.OriginName(n)))
                    fwrite(this.Handle,0,'uint8','l');
                end
                %  44; double values[nspectra]
                fwrite(this.Handle,this.OriginValues(n,:),'uint64','l');
            end
            
            % Close file handle
            fclose(this.Handle);
        end
    end
end