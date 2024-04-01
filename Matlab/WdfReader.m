% WdfReader  Provides access to selected data in WDF files.
%
% Renishaw's WDF file format is designed for storing Raman spectral data.
% The WdfReader class provides access to a subset of the data in
% a WDF file, including:
%   * The spectral data;
%   * Data origin lists (see below for more details); and
%   * Key file metadata (such as measurement title and username).
%   * Saturated spectra and cosmic ray removed spectra
%   * Loadings and explained variance from principal components analysis.
% It allows modified spectra and custom masks to be written back to the 
% file.
%
% wdf = WdfReader(FILENAME, fileAccess) creates a new WdfReader object,
% opening the WDF file specified by FILENAME with fileAccess = 'rb'
% (read-only) or 'r+b' (read & write). 
%
% wdf.Close() closes the file.
%
%
% Accessing spectral data
% =======================
%
% SPECTRA = wdf.GetSpectra(IndexStart, IndexEnd) reads the specified
% spectra from the file.  Indices are 1-based.
%
% wdf.WriteSpectra(IndexStart, Data) writes the spectrum intensity lists 
% in Data to the file, starting at the specified spectrum index. The 
% length of the new spectra must be the same as those that are being 
% replaced and the number of spectra in Data must be not more than the 
% number of spectra in the file with indices >= IndexStart. 
% 
% Spectral datasets are associated with an 'X-list' and 'Y-list'; for
% normal spectra the Y-list contains only one element and can be ignored.
% The meaning of the data in these lists can be determined by the relevant
% properties of the WdfReader object.
%
% XLIST = wdf.GetXList() retrieves the X-list data.
% YLIST = wdf.GetYList() retrieves the Y-list data.
%
% For large files that cannot fit into memory, data can be processed in
% manageable consecutive chunks using the following set of functions:
%    wdf.StartChunkwiseReading();
%    [X, indices] = wdf.GetNextChunk(numberOfSpectra);
%    trueOrFlase  = wdf.AreMoreChunks();
%    fraction     = wdf.GetChunkwiseProgress();
% (note that numberOfSpectra is optional for GetNextChunk()).
%
% wdf.WriteXList(XList) overwrites the entire wavenumber list.  This is 
% only permissible if the new XList is the same size as the list that 
% it is replacing.
%
%
% WDF file metadata
% =================
%
% WdfReader exposes the following properties:
%    Title               The measurement title.
%    Username            The name of the user who acquired the data.
%    MeasurementType     The measurement type of the data in the file.
%    ScanType            The scan type used to acquire the data.
%    LaserWavenumber     The wavenumber of the laser used.
%    Count               The actual number of spectra that are stored in
%                        the file (the number of spectra collected).
%    SpectralUnits       A WiREDataUnit value indicating the units of
%                        measurement associated with the spectral data.
%    XListType           A WiREDataType value indicating the type (or
%                        meaning) of the values stored in the X-list.
%    XListUnits          A WiREDataUnit value indicating the units of
%                        measurement associated with the X-list data.
%    YListType           A WiREDataType value indicating the type (or
%                        meaning) of the values stored in the Y-list.
%    YListUnits          A WiREDataUnit value indicating the units of
%                        measurement associated with the Y-list data.
%    PointsPerSpectrum   The number of elements in each spectrum / dataset.
%                        Equal to (XListLength * YListLength).
%    DataOriginCount     The number of Data Origin Lists in the file.
%    Capacity            The maximum number of spectra that can be stored
%                        in the file.
%    ApplicationName     The name of the software used to acquire the data.
%    ApplicationVersion  The software version used to acquire the data (as
%                        a 4-element vector: [major, minor, patch, build]).
%    XListLength         The number of elements in the X-list.
%    YListLength         The number of elements in the Y-list.
%    AccumulationCount   The number of accumulations co-added for each
%                        spectrum.
%
% Accessing Data Origin List values
% =================================
%
% WDF files may contain one or more Data Origin Lists.  Each list stores a
% unique type of data (such as time, or position in the X-axis) for each
% spectrum / dataset in the file.  Each Data Origin List is identified by
% the WiREDataType of the values it stores.  Lists are either 'primary'
% (important) or 'alternate' (less important).  The data in each list can
% either be stored as floating-point or integer values, and the caller must
% use the correct method to access the data (otherwise meaningless values
% will be returned); this can typically be inferred from the WiREDataType
% value associated with the list.
% 
% INFO = wdf.GetOriginListInfo() returns an Nx4 cell-array where each row
% describes a Data Origin List in the file.  The columns contain:
%   (1) A LOGICAL value indicating if the list is a 'primary' origin list;
%   (2) The list data type (as a WiREDataType value);
%   (3) The units of measurement (as a WiREDataUnit value); and
%   (4) The list name.
%
% V = wdf.GetOriginListValues(WiREDataType, IndexStart, IndexEnd) reads the
% specified values from the Data Origin List with the specified data type,
% treating the binary data in the origin list as double-precision floating
% point values.  Indices are 1-based.
%
% V = wdf.GetOriginListInfo(WiREDataType, IndexStart, IndexEnd) is
% equivalent to GetOriginListValues(), but treats the binary data in the
% origin list as 64-bit integer values.
%
% Times = wdf.GetOriginTimes(IndexStart,IndexEnd) Reads the timestamps that 
% for spectra in the specified index range. 
%
% wdf.WriteOriginTimes(Timestamps,IndexStart,IndexEnd) overwrites the
% consecutive spectrum time stamps between IndexStart and IndexEnd. 
%
% [Xcoord, Ycoord] = wdf.GetOriginCoords(IndexStart,IndexEnd) obtains lists 
% of the X and Y coordinates at which data were collected for the specified 
% spectrum index range. 
% 
% [Zcoord] = wdf.GetOriginZCoords(IndexStart,IndexEnd) obtains a list of 
% the Z coordinates at which data were collected for the specified spectrum 
% index range.
%
%
% Accessing Origin Flags
% ======================
%
% The origin flags array is a list of 64-bit values having an entry for each 
% spectrum in the wdf file. The bits correspond to various flags set for the 
% corresponding spectrum. Those bits important for use outside Renishaw are 
% detailed below. 
%   Bit     Flag meaning
%   1       Saturated
%   2       Error
%   3       Cosmic Ray removed
%   49      Current mask
%   50 - 64 Custom mask
%
% [Saturated CosmicRay] = wdf.GetOriginFlags(IndexStart,IndexEnd) will return
% two lists. The first is a list of all spectra that have been flagged as 
% suffering from saturation. The second is a list of all spectra that have
% had cosmic rays removed. INDEXSTART and INDEXEND define the range of
% spectra to consider, where INDEXSTART is 1 or greater.
%
% SATURATED and COSMICRAY return vectors of numbers that correspond to the
% Nth spectrum. For instance an appearance of the number 59 means that the
% 59th spectrum has been affected by the relevant flag.
%
% [FlaggedIndices] = GetOriginFlagsBit(FlagBits,IndexStart,IndexEnd) 
% gets an array of spectrum indices for spectra which have the flag 
% specified by FlagBits set. E.g. CurrentMask = 
% GetOriginFlagsBit(49,IndexStart,IndexEnd) will get a list of all 
% spectra in range [IndexStart,IndexEnd] that are currently masked.  
%
% [Success] = AddMaskToFlags(Mask, MaskBit) is intended for adding custom 
% masks into the file origin flags. Mask is a column vector of spectrum 
% indices corresponding to spectra to be masked. MaskBit is a scalar or 
% row vector spcifying the bits to mask. MaskBit should be in the range  
% [50,64], corresponding to the custom mask area of the origin flags.
%
% [Success] = RemoveFlagBits(MaskBit) is intended for clearing custom mask 
% bits. MaskBit is a row vector of integers between 50 and 64 inclusive 
% corresponding to the custom masks area of the origin flags. E.g. 
% RemoveFlagBits(50:64) will clear all custom flags. 
%
%
% Accessing the Loadings Data from Principal Component Analysis
% =============================================================
%
% [Loadings VarianceExplained] = wdf.getLoadingsData(PCANumber) will
% return both the loadings for all principal components computed in the
% WiRE software and the corresponding percentage of the overall variance
% that each principal component represents.
% 
% %LOADINGS is a 2D array in which the first column corresponds to the
% Raman shift in wavenumbers, and each subsequent column the corresponding
% loadings for each principal component at each Raman shift. The second
% column is PC1, the third PC2 etc, listed in order of decreasing
% percentage of variance explained. 
%
% VARIANCEXPLANINED is a vector listing each of the percentage variance 
% explained values. THe first entry corresponds to PC1, the second to PC2
% etc.
%
% PCANUMBER corresponds to the specific PCA loadings data that are required.
% If one has performed PCA three times, there will be three sets of PCA data
% stored. If the user wanted to extract the loadings from the first
% analysis, PCANUMBER would be equal to 1.
%
% Accessing all the Map Data
% ==========================
%
% [MapData, MapOutput] = wdf.getMapData() will return all the map data in
% the specified wdf file, along with the p-set description of each map.
%
% MAPDATA is a 2D array in row order. Each row contains the map data, and
% there are k rows, for k maps performed.
%
% MAPOUTPUT is a cell array containing the p-set of each map found in the
% file. This allows selection of specific maps, depending on the content of
% the p-set. A cell array is used because the p-set is not a standard
% format between different map types.
%
% Accessing the WMAP Block
% ========================
%
% WMAP = wdf.GetWMAPBlock() returns a struct containing the following
% fields: Location, StepSize, numPoints, linefocus_size and flag. These are
% accessed using either WMAP.Location (if it is already run), or
% wdf.GetWMAPBlock().Location.
%
% LOCATION is the X, Y, Z position of the map.
%
% STEPSIZE is the step size used in X, Y, Z. If the motion is not performed
% in the Z-axis, this will default to 1.
%
% NUMPOINTS is the number of points in the scan (ie. pixels) in X, Y and Z.
% If no steps are taken in the Z-axis, this will default to 1.
%
% LINEFOCUS_SIZE is the size of the line focus used, in X Y and Z. The axes
% not used default to 0.
%
% FLAG is a flag indicating the type of scan performed, usually given by
% rastor or snake, followed by rowmajor or columnmajor, though more
% complicated setups are possible.
%
% Accessing the White Light Image
% ===============================
%
% [WhiteLightImage,XCoords,YCoords]=wdf.GetWLImage() returns the white
% light image found in the wdf file.
%
% WHITELIGHTIMAGE is a 3d array, of size (NumberOfPixelsInX * NumberOfPixelsInY * 3).
% This represents an RGB image the size of the map.
%
% XCOORDS is the sequence of x-values for each pixel of the white light
% image.
%
% YCOORDS is the sequence of y-values for each pixel of the white light
% image.
%
% Note, alternatively the image is length(XCoords) * length(YCoords) * 3.
%
% Unsupported methods
% ===================
%
% The following methods are exposed by WdfReader but are intended for use
% only by Renishaw:
%    GetBlockID(),  GetFullBlockID(),  GetNextAvailableBlockUID(),
%    GetPrimaryOriginLists(), GetCurrentMaskInfo().

% Copyright 2012 - 2022 Renishaw plc.
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

classdef WdfReader < handle
    %% Private fields
    properties (Access = protected)
        m_NextSpectrumIndex = 1;
        Handle = -1;
        XListOffset = int64(-1);
        YListOffset = int64(-1);
        OriginsOffset = int64(-1);
        DataOffset = int64(-1);
        OriginListInfo = cell(0, 4);
        PCAOffset = int64(-1);
    end;
    
    %% Public properties relating to WDF metadata
    properties (GetAccess = public, SetAccess = protected)
        Title = '';
        Username = '';
        MeasurementType = WiREMeasurementType.Unspecified;
        ScanType = WiREScanType(0);
        LaserWavenumber = nan;
        Count = nan;
        SpectralUnits = WiREDataUnit.Arbitrary;
        XListType = WiREDataType.Arbitrary;
        XListUnits = WiREDataUnit.Arbitrary;
        YListType = WiREDataType.Arbitrary;
        YListUnits = WiREDataUnit.Arbitrary;
        PointsPerSpectrum = nan;
        DataOriginCount = nan;
        Capacity = nan;
        ApplicationName = '';
        ApplicationVersion = nan(1, 4);
        XListLength = nan;
        YListLength = nan;
        AccumulationCount = nan;
        OPut = cell(1,2);
        TimeStamp = nan;
        FileAccess = 'rb';
    end;
    
    %% Constructor, destructor and Close
    methods
        % Constructor: opens an existing WDF file using the protected
        % constructor.  The passed-in fileaccess is used if supplied ('rb',
        % 'r+b' supported), defaults to 'rb' if not supplied.  
        function this = WdfReader(fileName, fileAccess)
            if (nargin > 1)
                if ~(strcmp(fileAccess,'rb') || strcmp(fileAccess,'r+b'))
                    throw(WdfError('WdfReader does not support %s file access.  Use rb for read access or r+b to write to the file.', fileAccess));
                end;
            else
                fileAccess = 'rb';
            end;
            this.ProtectedConstructor(fileName, fileAccess);
        end;
        
        % Closes the file.
        function Close(this)
            if (this.Handle ~= -1)
                fclose(this.Handle);
                this.Handle = -1;
            end;
        end;
        
        % Destructor: just calls the Close() method.
        function delete(this)
            this.Close();
        end;
    end;
    
    %% Protected helper methods
    methods (Access = protected)
        % Opens a WDF file with the specified access mode, reads the WDF
        % header, and validates the DATA, XLST, YLST and ORGN blocks.
        function ProtectedConstructor(this, fileName, fileAccess)
            % Open the file
            this.Handle = fopen(fileName, fileAccess);
            if (this.Handle == -1)
                throw(WdfError('Cannot open file "%s".', fileName));
            end;
            this.FileAccess = fileAccess;
            
            % Parse the header, locate the data block and check its size is
            % consistent with expected value.
            this.ReadWdfHeader();
            [this.DataOffset, ~, dataSize] = LocateBlock(this, 'DATA');
            if (this.DataOffset == -1)
                throw(WdfError('Cannot locate Spectral Data block.'));
            end;
            if (dataSize < (16 + (4 * this.PointsPerSpectrum * this.Capacity)))
                throw(WdfError('Spectral Data block size inconsistent with WDF header.'));
            end;
            
            % Locate the X-list and Y-list blocks, and validate their sizes.
            [this.XListOffset, ~, dataSize] = LocateBlock(this, 'XLST');
            if (this.XListOffset == -1)
                throw(WdfError('Cannot locate X-list block.'));
            end;
            if (dataSize < (24 + (4 * this.XListLength)))
                throw(WdfError('X-list block size inconsistent with WDF header.'));
            end;
            
            [this.YListOffset, ~, dataSize] = LocateBlock(this, 'YLST');
            if (this.YListOffset == -1)
                throw(WdfError('Cannot locate Y-list block.'));
            end;
            if (dataSize < (24 + (4 * this.YListLength)))
                throw(WdfError('Y-list block size inconsistent with WDF header.'));
            end;
            
            % Read the X and Y list types / units
            this.XListType = GetXListType(this);
            this.YListType = GetYListType(this);
            this.XListUnits = GetXListUnits(this);
            this.YListUnits = GetYListUnits(this);
            
            % If the Data Origin List count is non-zero, attempt to locate
            % then size-check the ORGN block.  Finally, read in the data
            % origin list info.
            if (this.DataOriginCount ~= 0)
                [this.OriginsOffset, ~, dataSize] = LocateBlock(this, 'ORGN');
                if (this.OriginsOffset == -1)
                    throw(WdfError('Cannot locate Data Origin List block.'));
                end;
                if (dataSize < (20 + (this.DataOriginCount * (24 + (this.Capacity * 8)))))
                    throw(WdfError('Data Origin List block size inconsistent with WDF header.'));
                end;
            end;
            this.ReadOriginListInfo();

        end;
           
        % Reads the WDF header block, extracting key header fields and
        % storing them in class properties.
        function ReadWdfHeader(this)
            % Confirm that the signature, version and size fields
            % contain the expected values.
            this.AssertSeek(0, 'bof');
            blockID = this.AssertRead(1, 'uint32');
            blockUID = this.AssertRead(1, 'uint32');
            blockLength = this.AssertRead(1, 'uint64');
            if (~isequal(blockID, this.GetBlockID('WDF1')) || ...
                    ~(isequal(blockUID, 0) || isequal(blockUID, 1)) || ...
                    ~isequal(blockLength, 512))
                throw(WdfError('File does not use a recognised WDF format / version.'));
            end;
            
            % Read key fields from the main file header
            this.AssertSeek(60, 'bof');
            this.PointsPerSpectrum = this.AssertRead(1, 'uint32');
            this.Capacity = this.AssertRead(1, 'uint64');
            this.Count = this.AssertRead(1, 'uint64');
            this.AccumulationCount = this.AssertRead(1, 'uint32');
            this.YListLength = this.AssertRead(1, 'uint32');
            this.XListLength = this.AssertRead(1, 'uint32');
            this.DataOriginCount = this.AssertRead(1, 'uint32');
            this.ApplicationName = this.ReadUtf8String(24);
            this.ApplicationVersion = this.AssertRead([1 4], 'uint16');
            this.ScanType = WiREScanType(this.AssertRead(1, '*uint32'));
            this.MeasurementType = WiREMeasurementType(this.AssertRead(1, '*uint32'));
            this.AssertSeek(152, 'bof');
            this.SpectralUnits = WiREDataUnit(this.AssertRead(1, '*uint32'));
            this.LaserWavenumber = this.AssertRead(1, 'float32');
            this.AssertSeek(208, 'bof');
            this.Username = this.ReadUtf8String(32);
            this.Title = this.ReadUtf8String(160);
            this.AssertSeek(136,'bof');
            this.TimeStamp = this.AssertRead(1, 'uint64');
        end;
        

        
        % Searches for a specific block within the WDF file, by ID and
        % optionally UID.  Returns the location (offset), UID and length
        % (in bytes) of the first matching block found, or -1 if no
        % matching block was located.
        function [Location, BlockUID, BlockSize] = LocateBlock(this, TargetID, TargetUID)
            % Set default return values if block is not found
            Location = -1;
            BlockUID = [];
            BlockSize = uint64(0);
            
            % Initialise search
            found = false;
            offset = int64(512);
            TargetID = this.GetBlockID(TargetID);
            
            % Iteratively walk over all data-blocks in the file
            while ((~found) && (feof(this.Handle) == 0))
                % Attemp to jump to next block
                if (fseek(this.Handle, offset, 'bof') ~= 0)
                    return;
                end;
                
                % Read block header, and abort search if header fields not
                % read successfully
                blockID = fread(this.Handle, 1, 'uint32');
                uid = fread(this.Handle, 1, 'uint32');
                bSize = fread(this.Handle, 1, 'uint64');
                if (isempty(blockID) || isempty(uid) || (numel(bSize) ~= 1) || (bSize < 16))
                    return;
                end;
                
                % Check if we have found the requested block
                if (nargin == 3)
                    found = (isequal(TargetID, blockID) && isequal(TargetUID, uid));
                else
                    found = isequal(TargetID, blockID);
                end;
                
                % Store block info if match found, otherwise prepare to
                % look for the next block
                if (found)
                    Location = offset;
                    BlockUID = uid;
                    BlockSize = bSize;
                else
                    offset = offset + bSize;
                end;
            end;
        end;
        
        % Searches for all blocks within the WDF file by ID.
        % It returns the set of locations, block UIDs associated with the given ID, 
        % block sizes (as arrays) and the number found (as a number).
        function [Locations, BlockUIDs, BlockSizes, NumberFound] = LocateAllBlock(this, TargetID)
            % Set default return values if block is not found
            Locations = [-1];
            BlockUIDs = [];
            BlockSizes = [0];
            NumberFound = 0;
            
            % Initialise search
            offset = int64(512);
            TargetID = this.GetBlockID(TargetID);
            
            % Iteratively walk over all data-blocks in the file
            while (feof(this.Handle) == 0)
                % Attemp to jump to next block
                if (fseek(this.Handle, offset, 'bof') ~= 0)
                    return;
                end;
                
                % Read block header, and abort search if header fields not
                % read successfully
                blockID = fread(this.Handle, 1, 'uint32');
                uid = fread(this.Handle, 1, 'uint32');
                bSize = fread(this.Handle, 1, 'uint64');
                if (isempty(blockID) || isempty(uid) || (numel(bSize) ~= 1) || (bSize < 16))
                    return;
                end;
                
                
                % Check if we have found the requested block
                found = isequal(TargetID, blockID);
                if(found)
                   NumberFound = NumberFound + 1; 
                end
                
                % Store block info if match found, otherwise prepare to
                % look for the next block
                if (found)
                    Locations(NumberFound) = offset;
                    BlockUIDs(NumberFound) = uid;
                    BlockSizes(NumberFound) = bSize;
                    offset = offset + bSize;
                else
                    offset = offset + bSize;
                end;
            end;
        end;
        
        % Reads the specified number of bytes and converts them from UTF8
        % to a Matlab string, with trailing NULLs and whitespace trimmed.
        function s = ReadUtf8String(this, Length)
            s = this.AssertRead([1 Length], 'uint8=>char');
            s = deblank(char(unicode2native(s, 'UTF-8')));
        end;
        
       
        % Seeks to the requested position in the WDF file, raising an error
        % if the seek operation is unsuccessful.
        function AssertSeek(this, Position, Origin)
            if (fseek(this.Handle, Position, Origin) ~= 0)
                throwAsCaller(WdfError('Failed to seek to requested position within file.'));
            end;
        end;
        
        % Reads the requested data from the WDF file, raising an error if
        % the actual number of elements read is fewer than requested.
        function [Data] = AssertRead(this, RDims, Precision)
            [Data, readCount] = fread(this.Handle, RDims, Precision);
            if (readCount ~= prod(RDims))
                throwAsCaller(WdfError('Failed to read requested data from file.'));
            end;
        end;
        
        % Writes the specified data to the WDF file, in column order, and
        % raising an error if the actual number of elements written is
        % fewer than requested.
        function AssertWrite(this, Data, Precision)
            [writeCount] = fwrite(this.Handle, Data, Precision);
            if (writeCount ~= numel(Data))
                errorStr = 'Failed to write requested data to File.';
                if ~(strcmp(this.FileAccess,'r+b'))
                    errorStr = strcat(errorStr, '\nFile has not been opened for write access.');
                end;
                throwAsCaller(WdfError(errorStr));
            end;
        end;
        
        % Helper function used to read data origin list values.
        function [Data] = ReadOriginListData(this, ListType, IndexStart, IndexEnd, Precision)
            % Validate inputs.
            if (nargin ~= 5)
                throwAsCaller(WdfError('ListType, IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throwAsCaller(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throwAsCaller(WdfError('IndexEnd is out-of-range.'));
            end;
            listIndex = find(cellfun(@(x) isequal(x, ListType), this.OriginListInfo(:, 2)), 1);
            if (isempty(listIndex))
                throwAsCaller(WdfError('WDF file does not contain a Data Origin List with the specified type.'));
            end;
            
            % Read and return the data.
            count = (IndexEnd - IndexStart) + 1;
            listSize = 24 + (8 * this.Capacity);
            offset = this.OriginsOffset + 20 + ((listIndex - 1) * listSize) + 24 + ((IndexStart - 1) * 8);
            this.AssertSeek(offset, 'bof');
            Data = this.AssertRead(count, Precision);
        end;
        
        % Helper function used to read data origin list values.
        function WriteOriginListData(this, Data, ListType, IndexStart, IndexEnd, Precision)
            % Validate inputs.
            if (nargin ~= 6)
                throwAsCaller(WdfError('Data, ListType, IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throwAsCaller(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throwAsCaller(WdfError('IndexEnd is out-of-range.'));
            end;
            listIndex = find(cellfun(@(x) isequal(x, ListType), this.OriginListInfo(:, 2)), 1);
            if (isempty(listIndex))
                throwAsCaller(WdfError('WDF file does not contain a Data Origin List with the specified type.'));
            end;
            
            % write the new data.
            listSize = 24 + (8 * this.Capacity);
            offset = this.OriginsOffset + 20 + ((listIndex - 1) * listSize) + 24 + ((IndexStart - 1) * 8);
            this.AssertSeek(offset, 'bof');
            this.AssertWrite(Data, Precision);
        end;
        
        %Check that the list origin actually exists. Currently used to
        %check the origin list for the presence of (x,y) coordinates.
        function Check = CheckOriginListData(this, ListType, IndexStart, IndexEnd, Precision)
            % Validate inputs.
            if (nargin ~= 5)
                throwAsCaller(WdfError('ListType, IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throwAsCaller(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throwAsCaller(WdfError('IndexEnd is out-of-range.'));
            end;
            listIndex = find(cellfun(@(x) isequal(x, ListType), this.OriginListInfo(:, 2)), 1);
            if (isempty(listIndex))
                Check = false;
                return
            end;
            Check = true;
        end;
    end;
    

   
    
    %% Public methods that access key WDF data
    methods (Access = public)
        % Get the X-list type
        function [Result] = GetXListType(this)
            this.AssertSeek(this.XListOffset + 16, 'bof');
            Result = WiREDataType(this.AssertRead(1, '*uint32'));
        end;
        
        % Get the X-list units
        function [Result] = GetXListUnits(this)
            this.AssertSeek(this.XListOffset + 20, 'bof');
            Result = WiREDataUnit(this.AssertRead(1, '*uint32'));
        end;
        
        % Get the Y-list type
        function [Result] = GetYListType(this)
            this.AssertSeek(this.YListOffset + 16, 'bof');
            Result = WiREDataType(this.AssertRead(1, '*uint32'));
        end;
        
        % Get the Y-list units
        function [Result] = GetYListUnits(this)
            this.AssertSeek(this.YListOffset + 20, 'bof');
            Result = WiREDataUnit(this.AssertRead(1, '*uint32'));
        end;
        
        % Creates and stores a cell-array containing information about the
        % Data Origin Lists in the WDF file (one row per data origin list).
        function ReadOriginListInfo(this)
            ListInfo = cell(this.DataOriginCount, 4);
            listSize = 24 + (8 * this.Capacity);
            for n = 1:this.DataOriginCount
                this.AssertSeek(this.OriginsOffset + 20 + ((n - 1) * listSize), 'bof');
                listInfo = this.AssertRead([1 2], '*uint32');
                ListInfo{n, 1} = bitget(listInfo(1), 32) ~= 0;
                if (bitget(listInfo(1), 31) ~= 0)
                    ListInfo{n, 2} = WiREDataType(99); % custom
                else
                    ListInfo{n, 2} = WiREDataType(bitset(listInfo(1), 32, 0));
                end;
                ListInfo{n, 3} = WiREDataUnit(listInfo(2));
                ListInfo{n, 4} = this.ReadUtf8String(16);
            end;
            this.OriginListInfo = ListInfo;
        end;
   end; 
   

    
    %% Access to spectra in the file, plus X- and Y-list data
    methods
        % Read the X-list
        function [XList] = GetXList(this)
            this.AssertSeek(this.XListOffset + 24, 'bof');
            XList = double(this.AssertRead([1 this.XListLength], 'single'));
        end;
        
        % Read the Y-list
        function [YList] = GetYList(this)
            this.AssertSeek(this.YListOffset + 24, 'bof');
            YList = double(this.AssertRead([1 this.YListLength], 'single'));
        end;
        
        % Import a range of spectra, specified by start & end indices.
        function [Data] = GetSpectra(this, IndexStart, IndexEnd,type)
            % Validate inputs.
            if (nargin < 3)
                throw(WdfError('IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throw(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throw(WdfError('IndexEnd is out-of-range.'));
            end;
            
            % Read the data, and confirm the expected quantity was read
            offset = (this.DataOffset + 16) + ((IndexStart - 1) * 4 * this.PointsPerSpectrum);
            nRowsToRead = (IndexEnd - IndexStart) + 1;
            this.AssertSeek(offset, 'bof');
            Data = this.AssertRead(nRowsToRead * this.PointsPerSpectrum, 'single=>single');
            if nargin < 4 || strcmp(type,'double') 
                % Re-shape the data and convert to double-type
                Data = double(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
            elseif strcmp(type,'single')
                 Data = single(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
            else
                 throw(WdfError('Type must be single or double'))
            end
        end;
        
        % Write spectra back to the file, starting at the specified start
        % spectrum index.
        function WriteSpectra(this, IndexStart, Data)
            if (nargin ~= 3)
                throw(WdfError('IndexStart and Data must be specified.'));
            end;
            if (size(Data, 2) ~= this.PointsPerSpectrum)
                throw(WdfError('Number of columns in Data must match PointsPerSpectrum.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throw(WdfError('IndexStart is out-of-range.'));
            end;
            if (IndexStart + size(Data, 1) - 1 > this.Count)
                throw(WdfError('IndexStart and number of rows in Data would result in writing beyond end of file.'));
            end;
            
            % Write the data.
            offset = (this.DataOffset + 16) + ((IndexStart - 1) * 4 * this.PointsPerSpectrum);
            this.AssertSeek(offset, 'bof');
            this.AssertWrite(Data', 'float32');
        end;
        
        % Write X (frequency) list back to the file.  The length of the
        % original X list must be maintained. 
        function WriteXList(this, XList)
            % Argument checking
            if (nargin ~= 2)
                throw(WdfError('XList to write to file must be specified.'));
            end;
            if (size(XList,1) ~= 1 || size(XList,2) ~= this.XListLength)
                throw(WdfError('XList must be a row vector with the same number of points as the original X list.'));
            end;
            
            % Write the X list.
            this.AssertSeek(this.XListOffset + 24, 'bof');
            this.AssertWrite(XList, 'float32');
        end;
        
        
        % Initiates chunk-wise reading of the spectral data
        function StartChunkwiseReading(this)
            this.AssertSeek(this.DataOffset + 16, 'bof');
            this.m_NextSpectrumIndex = 1;
        end;
        
        % Reads the next chunk of spectra from the file; NumberOfSpectra is
        % optional.
        function [Data, Indices] = GetNextDataChunk(this, NumberOfSpectra)
            % Determine how many rows (spectra) to read into this chunk.
            if (nargin == 2)
                if (NumberOfSpectra < 1)
                    throw(WdfError('NumberOfSpectra must be greater-than-or-equal-to 1.'));
                end;
                nRowsToRead = NumberOfSpectra;
            else
                nRowsToRead = 4096;
            end;
            
            % If necessary, limit the number of rows to read to the number
            % of remaining available rows.
            if (this.m_NextSpectrumIndex + nRowsToRead > (this.Count + 1))
               nRowsToRead = (this.Count - this.m_NextSpectrumIndex) + 1;
            end;
            
            % Check there are *some* more spectra remaining.
            if (nRowsToRead < 1)
                throw(WdfError('All spectra have already been read in previous chunks.'));
            end;
            
            % Read the data.
            offset = (this.DataOffset + 16) + ((this.m_NextSpectrumIndex - 1) * this.PointsPerSpectrum * 4);
            this.AssertSeek(offset, 'bof');
            Data = this.AssertRead(nRowsToRead * this.PointsPerSpectrum, 'single');
            
            % Calculate the spectrum indices of the data we have just read.
            if (nargout > 1)
                Indices = this.m_NextSpectrumIndex + colon(0, nRowsToRead - 1);
            end;
                
            % Re-shape the data, convert to double type, and update the
            % next read-position marker.
            Data = double(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
            this.m_NextSpectrumIndex = this.m_NextSpectrumIndex + nRowsToRead;
        end;
        
        % Returns a boolean indicating if there are more chunks of spectra
        % to be read from the file.
        function [MoreChunks] = AreMoreChunks(this)
            MoreChunks = (this.m_NextSpectrumIndex <= this.Count);
        end;
        
        % Returns a value in the range [0-1] that indicates the fraction of
        % spectra already imported via chunk-wise reading.
        function [Fraction] = GetChunkwiseProgress(this)
            Fraction = (this.m_NextSpectrumIndex - 1) / this.Count;
        end;
    end;
    
    %% Access to Data Origin Lists
    methods
        % Returns a cell-array containing information about the Data Origin
        % Lists in the WDF file (one row per data origin list).
        function [ListInfo] = GetOriginListInfo(this)
            ListInfo = this.OriginListInfo;
        end;
        
        % Gets a cell-array containing a list of primary (non-alternate)
        % data origin list data types.
        function [PrimaryOriginLists] = GetPrimaryOriginLists(this)
            PrimaryOriginLists = cell(sum(cell2mat(this.OriginListInfo(:, 1))), 2);
            
            % Cycle over the data origin lists, and gather all primary
            % lists into the result.
            resultIndex = 1;
            for n = 1:this.DataOriginCount
                if (this.OriginListInfo{n, 1})
                    PrimaryOriginLists{resultIndex, 1} = sprintf('DataList%d', resultIndex - 1);
                    PrimaryOriginLists{resultIndex, 2} = int32(this.OriginListInfo{n, 2});
                    resultIndex = resultIndex + 1;
                end;
            end;
        end;
        
        % Gets a range of values from a Data Origin List, specified by list
        % data type, treating the binary data as double-precision floating
        % point values.
        function [Data] = GetOriginListValues(this, ListType, IndexStart, IndexEnd)
            Data = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
        end;
        
        % Gets a range of values from a Data Origin List, specified by list
        % data type, treating the binary data as 64-bit integer values.
        function [Data] = GetOriginListValuesInt(this, ListType, IndexStart, IndexEnd)
            Data = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'int64=>double');
        end;
        
        %Finds out which spectra have a saturated spectrum or cosmic ray 
        %removal flag associated with them. Returns the indices of the
        %flagged datasets corresponding to the (1-based) data array 
        %returned by the class constructor.
        function [Saturated, CosmicRay] = GetOriginFlags(this,IndexStart,IndexEnd)
            ListType = WiREDataType.Flags;
            Data = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'int64=>double');
            S = 7*ones(length(Data),1); Data = bitand(Data,S);
            %Find any byte with bit 1 set - DataSaturated flag
            A = find(Data==1); B = find(Data==3); C = find(Data==5); D = find(Data==7); 
            Saturated = sort([(A); (B); (C); (D)]); clear A B C D
            %Find any byte with bit 4 set - Cosmic Ray Removed flag
            A = find(Data==4); B = find(Data==5); C = find(Data==6); D = find(Data==7);
            CosmicRay = sort([(A); (B); (C); (D)]); clear A B C D
        end;
        
        %Extract all of the timestamps that are associated with the file.
        %Each spectrum has a timestamp associated with it.
        function Timestamps = GetOriginTimes(this,IndexStart,IndexEnd)
            ListType = WiREDataType.Time;
            Timestamps = this.ReadOriginListData(ListType, IndexStart,IndexEnd,'int64');
            unix_time = (Timestamps-116444736000000000)/10000000;
            Timestamps = datestr(((unix_time+4)/86400 + datenum(1970,1,1)),'dd/mm/yyyy HH:MM:SS');
        end
        
        %Overwrite time stamps located in the origin list. Time stamps
        %should be input in the 'dd/mm/yyyy HH:MM:SS' format. They will be
        %converted to the Windows date format.
        function WriteOriginTimes(this,Timestamps,IndexStart,IndexEnd)
            ListType = WiREDataType.Time;
            Timestamps = datenum(Timestamps,'dd/mm/yyyy HH:MM:SS') - datenum(1970,1,1);  
            unix_time = (86400*Timestamps) -4; 
            Timestamps = 116444736000000000 + unix_time*10000000;
            WriteOriginListData(this, Timestamps, ListType, IndexStart, IndexEnd, 'int64')
        end
        
        % Add a mask flag to the flags origin list.  Mask is a column 
        % vector of indices into the Spectra array to specify which spectra    
        % are to be masked. MaskBit is a row vector spcifying the bits to 
        % mask in the range 50 to 64 inclusive, corresponding to the 
        % custom mask area of the origin flags.   
        function [Success] = AddMaskToFlags(this, Mask, MaskBit)
            Success = false;
            % Validate inputs.
            if (nargin < 3)
                throwAsCaller(WdfError('AddMaskToFlags: Mask and MaskBit must be specified.'));
            end;
            Success = SetOriginFlagsBit(this, MaskBit, 1, Mask);
        end;

        % Clear the masks for all spectra at bits specified by MaskBit. 
        % MaskBit is a row vector of integers between 50 and 64 inclusive 
        % corresponding to the custom masks area of the origin flags.  
        % MaskBit = 50:64 will clear all custom flags
        function [Success] = RemoveFlagBits(this, MaskBit)
            Success = false;
            % Validate inputs.
            if(nargin < 2)
                throwAsCaller(WdfError('RemoveFlagBits: MaskBit must be specified.'))
            end
            Success = SetOriginFlagsBit(this, MaskBit, 0);
        end;
 
        % Get an array of flagged spectrum indices. FlagBits is a scalar 
        % or a row vector of integers between 1 and 64 inclusive. 
        % Optional inputs IndexStart and IndexEnd are scalar integers 
        % specifying the range of spectra for which the flag value is 
        % required, which default to 1 and wdf.Count respectively. The 
        % FlaggedIndices output is a column vector corresponding to the 
        % indices of spectra in the specified range for which ANY of the 
        % specified bits are set. FlaggedIndices will be empty if none 
        % of the specified bits are set for any spectrum in the range, 
        % or if there is no origin data in the file.  
        function [FlaggedIndices] = GetOriginFlagsBit(this, FlagBits, IndexStart, IndexEnd)
            % IndexStart and IndexEnd are optional
            if(nargin < 4)
                IndexEnd = this.Count();
            end
            if(nargin < 3)
                IndexStart = 1;
            end
            
            % Input checking
            if(nargin < 2)
                throwAsCaller(WdfError('GetOriginFlagsBit: MaskBit must be specified.'))
            end
            if (size(FlagBits,1) ~= 1 || any(FlagBits - floor(FlagBits)))
                throwAsCaller(WdfError('GetOriginFlagsBit: MaskBit must be a scalar or row vector of integer values.'));
            end
            if (any(FlagBits < 1) || any(FlagBits > 64))
                throwAsCaller(WdfError('GetOriginFlagsBit: all elements of MaskBit must be between 1 and 64 inclusive.'));
            end
            if (~isscalar(IndexStart) || (IndexStart ~= floor(IndexStart)))
                throwAsCaller(WdfError('GetOriginFlagsBit: IndexStart must be a scalar integer.'));
            end
            if (~isscalar(IndexEnd) || (IndexEnd ~= floor(IndexEnd)))
                throwAsCaller(WdfError('GetOriginFlagsBit: IndexEnd must be a scalar integer.'));
            end
            if(IndexStart < 1 || IndexEnd > this.Count )
               throwAsCaller(WdfError(sprintf('GetOriginFlagsBit: IndexStart and IndexEnd must correspond to spectra in the file, so lie in the range [1,%d]',this.Count)))
            end
            if(IndexEnd < IndexStart)
               throwAsCaller(WdfError('GetOriginFlagsBit: IndexStart should be less than or equal to IndexEnd.'))
            end
            
            % Test origin data, specifically Flags, exists in the file.
            originListIndex = find(cellfun(@(x) isequal(x, WiREDataType.Flags), this.OriginListInfo(:, 2)), 1);
            if (isempty(originListIndex))
                FlaggedIndices = [];      % No origin / flags data.
            else
                % Get the origin flags.  
                flagsData = GetOriginListValuesInt(this, WiREDataType.Flags, IndexStart, IndexEnd);
                % Generate the mask from the data returned
                S = sum(bitset(0,FlagBits,1)) * ones(length(flagsData),1);
                FlaggedIndices = bitand(flagsData, S);
                FlaggedIndices = find(FlaggedIndices);
            end;
        end;
        
        %Extract all x and y coordinates. Note that streamline and
        %streamline HR have transposed (x,y) coordinates relative to each
        %other.
        function [Xcoord, Ycoord] = GetOriginCoords(this,IndexStart,IndexEnd)
            ListType = WiREDataType.SpatialX;
            Check = this.CheckOriginListData(ListType,IndexStart,IndexEnd,'double');
            if Check
                Xcoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
                ListType = WiREDataType.SpatialY;
                Ycoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
            else
                Xcoord = []; Ycoord = [];
            end
        end

        %Extract all z coordinates.
        function Zcoord = GetOriginZCoords(this,IndexStart,IndexEnd)
            ListType = WiREDataType.SpatialZ;
            Zcoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
        end
        
    end;

    methods(Access = private)
        % Sets custom mask flags in the flags origin list.  MaskBit is a 
        % row vector of integers between 50 and 64 inclusive.  Value is the 
        % value to which all bits are set, scalar = 0 of 1.  Optional input 
        % Mask may be supplied as a column vector of spectral indeces for 
        % which the mask will be (re)set to the specified value.  If not 
        % supplied the mask will be set for all spectra.  
        function [Success] = SetOriginFlagsBit(this, MaskBit, Value, Mask)
            Success = false;

            % Mask is optional - if not supplied, mask all spectra
            if(nargin == 3)
                Mask = (1:this.Count)';
            end
            
            % Input checking
            if(nargin < 3)
                throwAsCaller(WdfError('SetOriginFlagsBit: MaskBit and Value must be specified.'))
            end
            if (size(MaskBit,1) ~= 1 || any(MaskBit - floor(MaskBit)))
                throwAsCaller(WdfError('SetOriginFlagsBit: MaskBit must be a scalar or row vector of integer values.'));
            end
            if (any(MaskBit < 50) || any(MaskBit > 64))
                throwAsCaller(WdfError('SetOriginFlagsBit: all elements of MaskBit must be in the custom mask range of between 50 and 64 inclusive.'));
            end
            if (~isscalar(Value) || (Value ~= floor(Value)) || (Value < 0) || (Value > 1))
                throwAsCaller(WdfError('SetOriginFlagsBit: Value must be scalar and may only take the values 0 or 1.'));
            end
            if(size(Mask,2) ~= 1 || any(Mask - floor(Mask)))
               throwAsCaller(WdfError(sprintf('SetOriginFlagsBit: Mask must be a column vector of integer values.')))
            end
            if(any(Mask < 1) || any(Mask > this.Count))
               throwAsCaller(WdfError(sprintf('SetOriginFlagsBit: all elements of Mask must correspond to spectra in the file, so lie in the range [1,%d]',this.Count)))
            end
            
            listIndex = find(cellfun(@(x) isequal(x, WiREDataType.Flags), this.OriginListInfo(:, 2)), 1);
            if (isempty(listIndex))
                throwAsCaller(WdfError('SetOriginFlagsBit: WDF file does not contain a Flags type data origin list.'));
            end;
            
            % Read the flags data.
            listSize = 24 + (8 * this.Capacity);
            offset = this.OriginsOffset + 20 + ((listIndex - 1) * listSize) + 24;
            this.AssertSeek(offset, 'bof');
            flagsData = this.AssertRead(this.Count, 'int64=>double');

            % Construct a # spectra length column vector of values to mask 
            % (0 = do not mask, 1 = mask)
            m = zeros(this.Count,1);
            m(Mask) = 1;

            % Sets the bits as specified on a copy of flagsData for all  
            % elements (spectra), then selects whether or not to apply the 
            % modification or not based on the mask.   
            S = flagsData;
            for bit = MaskBit
                S = bitset(S,bit,Value);
            end
            S = S .* m;
            flagsData = int64(S .* m + flagsData .* (~m));
            
            %Write flags data back into the file
            this.AssertSeek(offset, 'bof');
            this.AssertWrite(flagsData, 'int64');
            Success = true;
        end;
    end    
    %% Access to loadings data stored in PCA maps.
    
    methods (Access = public)
        %Access the loadings from PCA analysis data. This set of methods
        %makes use of the new PSet parser. 
        function [Loadings, VarEx] = getLoadingsData(this,PCAno)
            qq = 2;
            [NM, ~] = this.getNumberofID('MAP ');
            for ii = 1:NM
                [blockOffset, ~, ~] = LocateBlock(this,'MAP ',ii);
                % Offset from start of block to PSET size element = parsing
                % start is 20 bytes: 'MAP ' + 12 bytes + 'PSET'                
                this.OPut = this.PSetParser(blockOffset + 20);
                CP = this.FindParent('PCAR',PCAno);   %Only needed for PCA maps
                Disc = this.ReadCellValue('overlayXListType');
                if CP == 1 && Disc~=-1
                    V(qq-1) = this.ReadCellValue('%VarianceExplained'); %#ok<AGROW>
                    L(:,1) = this.ReadCellValue('overlayXList');
                    L(:,qq) = this.ReadCellValue('overlaySpectrum'); %#ok<AGROW>
                    qq = qq + 1;
                end;
            end;
            if exist('V','var')
                VarEx = V; Loadings = L;
            else 
                throw(WdfError('Invalid PCANUMBER'))
            end
        end;
        
        function [MapData, MapOutput] = getMapData(this)
            [Anchor, ~, ~] = LocateBlock(this,'YLST');
            [NM, UIDs] = this.getNumberofID('MA');
            NM = uint32(NM);
            this.AssertSeek(Anchor,'bof');
            NS = uint32(this.Count); MapD = zeros(NM,NS); 
            for ii = 1:NM
                [blockOffset, ~, ~] = LocateBlock(this,'MA',UIDs(ii));
                this.OPut = this.PSetParser(blockOffset + 20);
                for jj = 1:size(this.OPut,1)
                    if(strcmp(this.OPut{jj,1}, 'Label'))
                        MapOutput{ii} = char(this.OPut{jj,2}); %#ok<AGROW>
                    end
                end
                %MapOutput{ii} = this.OPut;
                NSpec = this.AssertRead(1,'uint64');
                for NRead = 1:NSpec
                    MapD(ii,NRead) = this.AssertRead(1,'single');
                end
            end
            MapData = MapD;
        end
    end
        
    methods (Access = protected)
        % PSet parser is recursive without limit to PSet nesting.
        % It supports explicit and implicit key value names.  Ecplicit key
        % value names (Type 107 elements) may be defined before or after
        % their corresponding value.  If explicit value names are not found
        % in the PSet, implicit key names are looked-up from those defined
        % in WiREKeys. The key IDs are assumed to be unique.  
        % Returns a nested cell array reflecting the PSet structure.  
        function OP = PSetParser(this,PsetStart)

            OP = cell(1,2); % Output for each PSet is M x 2 cell array.
            keys = cell(1,2);
            this.AssertSeek(PsetStart,'bof');
            psetSize = this.AssertRead(1,'uint32'); % Total size of block.
            bytesToRead = psetSize;                 % Remaining bytes to read.
            nVal = 1;   nKey = 1;
            % by = bytes of the single element, whether key name or value. 
            
            while (bytesToRead > 0)
                [Type, F, k] = this.ReadFlagKey();
                switch Type
                    case {63, 99}       % char or bool
                        Txt = char('uint8');    by = 1;
                    case 115            % short int
                        Txt = char('uint16');   by = 2;
                    case 105            % int / long
                        Txt = char('uint32');   by = 4;
                    case 119            % ulong
                        Txt = char('uint64');   by = 8;
                    case 114            % float 
                        Txt = char('single');   by = 4;
                    case 113            % double 
                        Txt = char('double');   by = 8;
                    case 116            % time: not supported - discard the values
                        this.AssertRead(1,'uint32'); 
                        this.AssertRead(1,'uint32'); by = 8;
                    case 117            % string value
                        [A, by] = this.ReadPSetString();
                    case 107            % key name string
                        [A, by] = this.ReadPSetString();
                        keys{nKey,1} = k;
                        keys{nKey,2} = A;
                        nKey = nKey+1;
                    case 112            % PSet
                        % Size of PSet + 4 for the size element
                        by = this.AssertRead(1,'uint32') + 4;   
                        % Next PSet advances bytes read (= psetSize - bytesToRead) 
                        % + 4 bytes for the parent PSet size + 4 bytes for FlagKey. 
                        A = this.PSetParser(PsetStart + psetSize - bytesToRead + 8);
                    otherwise  
                        throw(WdfError('PSet data type not recognised'));
                end;
                
                % Read the value type elements that might be arrays, i.e.
                % all excl. strings, time and psets, which are read in the
                % above switch statement 
                if (Type==63)||(Type==99)||(Type==115)||(Type==105)||(Type)==119||(Type==114)||(Type==113) 
                    if F==0
                        A = this.AssertRead(1,Txt);
                    elseif F==128  
                        A = this.FillPSetArray(Txt); 
                        by = 4 + by * length(A);    % +4 bytes for array size.  
                    end;
                    % Cast numeric types from the default double to specific type.  
                    switch Type
                        case 63             % bool
                            A = logical(A);
                        case 99             % char
                            A = uint8(A);
                        case 115            % short int
                            A = uint16(A);
                        case 105            % int / long
                            A = uint32(A); 
                        case 119            % long long
                            A = uint64(A);
                        case 114            % float
                            A = single(A);
                    end;
                end;
                % Set the output for all value types.
                if (Type==63)||(Type==99)||(Type==115)||(Type==105)||(Type)==119|| ...
                        (Type==114)||(Type==113)||(Type==117)||(Type==112)
                    OP{nVal,1} = k;
                    OP{nVal,2} = A;
                    nVal=nVal+1;
                elseif (Type==116)      % time: not supported.
                    OP{nVal,1} = k;
                    OP{nVal,2} = 'not supported';
                    nVal=nVal+1;
                end;
                % Update remaining bytes to read, additional 4 bytes for the FlagKey element.  
                bytesToRead = bytesToRead - by - 4;
            end;
            
            % Finished read: transfer key values to output array matching
            % values if they exist in keys array and looking them up if not.
            if ~(nVal==1)       % something to match
                for ii = 1:size(OP,1)
                    if (~(nKey==1) && ~isempty(keys(:,1)) && ~isempty(OP{ii,1}))
                        rMatch = find([keys{:,1}] == OP{ii,1});
                    else
                        rMatch = [];
                    end
                    if isempty(rMatch)
                        try
                            OP{ii,1} = char(WiREKeys(OP{ii,1}));
                        catch
                        end
                    else
                        OP{ii,1} = keys{rMatch,2};
                    end
                end;
            end;
        end;
        
         %Function to read the flag and key each time. These occur for 
         %every PSet data type.
         function [Tq, Fq, kq] = ReadFlagKey(this)
             Tq = this.AssertRead(1,'uint8');
             Fq = this.AssertRead(1,'uint8');
             kq = this.AssertRead(1,'uint16');
         end;
         
         %Function to enable certain data types read in strings
         function [Str, by] = ReadPSetString(this)
             Len = this.AssertRead(1,'uint32');
             Str = this.ReadUtf8String(Len); by = Len+4;
         end;
         
         %Function to read in PSet data to an array. Reads the flags to
         %determine whether data is scalar or vector.
         function A = FillPSetArray(this,Tx)
              Leng = this.AssertRead(1,'uint32'); A = zeros(1,Leng);
              for ii = 1:Leng
                  A(ii) = this.AssertRead(1,Tx);
              end 
         end;
   
         %Function to determine whether or not the current MAP being 
         %parsed contains the correct parent PSet. Compares the parent ID
         %that we want with the actual PSet parent ID
         function  CorrectParent = FindParent(this,TargetID,TargetUID)
             PIndex = 0; CorrectParent = 0;
             TargetID = dec2hex(this.GetBlockID(TargetID));
             WantPID = uint64(str2double(strcat(num2str(TargetUID),num2str(TargetID)))); %In hex
             Strings = this.OPut(:,1);
             for x = 1:length(this.OPut(:,1))
                 St = Strings{x};
                 if strcmp(St,'parent') == 1; PIndex = x; end
             end;
             if PIndex ~=0
                 ActualPID = uint64(str2double(dec2hex(this.OPut{PIndex,2})));
                 CorrectParent = isequal(WantPID,ActualPID);
             end;
         end;
         
         %Helper function to read out the value of a cell, given the 
         %corresponding PSet Property. CellValue must be a string. If the 
         %CellValue field does not exist, then the function will return a 
         %value of -1.
         function  MapT = ReadCellValue(this,CellValue)
             Strings = this.OPut(:,1); PIndex = -1;
             for x = 1:length(this.OPut(:,1))
                 St = Strings{x};
                 if strcmp(St,CellValue) == 1; PIndex = x; end
             end;
             if PIndex == -1; MapT = -1;
             else MapT = (this.OPut{PIndex,2}); end;
         end
         
         %Function to determine the number of PSets that have a given ID in
         %the code
         function [NP, BlockUIDs] = getNumberofID(this,ID)
             [Anchor, ~, ~] = LocateBlock(this,'YLST');
             this.AssertSeek(Anchor,'bof');
                 [~, BlockUIDs, datasize] = LocateAllBlock(this,ID);                
             NP = length(datasize);
         end;
    end;

    %% Get the white light image from the file
    methods
        function [WL,x,y]=GetWLImage(this)
            this.AssertSeek(0, 'bof');                              %go to beginning of file
            [wl_offset, ~, dataSize] = LocateBlock(this, 'WHTL');   %find the white light location and get data size
            if (wl_offset == -1)                                    %if no data found throw error
                throw(WdfError('Cannot locate WL image.'));
            end;
            this.AssertSeek(wl_offset+16, 'bof');                   %go to the image position
            imgdata = fread(this.Handle, dataSize, '*uint8');       %read all the data as bytes
            %make use of imread function to avoid having to write a function to decompress jpgs
            fid = fopen('img_temp.jpg', 'w');                       %open a temporary jpg file with write access
            fwrite(fid, imgdata,'*uint8');                          %write the data as bytes
            fclose(fid);                                            %close the temporary file
            WL=imread('img_temp.jpg');                              %read the temporary file in using imread function
            % read the EXIF data and create x and y axes
            info=imfinfo('img_temp.jpg');                           
            FOV=info.UnknownTags(2).Value./info.UnknownTags(3).Value;
            x=linspace(info.UnknownTags(1).Value(1),info.UnknownTags(1).Value(1)+FOV(1),info.Width);
            y=linspace(info.UnknownTags(1).Value(2),info.UnknownTags(1).Value(2)+FOV(2),info.Height);
        end;        
    end;
    
    %% Read the WMAP block for orientation
    methods
        function Returned = GetWMapblock(this)
            this.AssertSeek(0, 'bof');                              %go to beginning of file
            [wmap_offset, ~, ~] = LocateBlock(this, 'WMAP');        %find the WMAPBlock location and get data size
            if (wmap_offset == -1)                                  %if no data found throw error
                throw(WdfError('Cannot locate WMap section.'));
            end;
            flagValue = this.AssertRead(1, 'uint32');
            ignoreValue = this.AssertRead(1, 'uint32');                  % Ignore next uint32
            Returned.Location = this.AssertRead(3, 'float32');
            Returned.StepSize = this.AssertRead(3, 'float32');
            Returned.numPoints = this.AssertRead(3, 'uint32');
            Returned.linefocus_size = this.AssertRead(1, 'uint32');
            flag = this.DecimalBaseConverter(flagValue,2);          % Convert the 32bit flag to a binary sequence
            Returned.flag = this.ConvertFlagToString(flag);         % Convert the binary sequence into concatenated string of flag descriptions
        end;
    end;

    %% Get the current mask info from the mask block. 
    % Method allows analysis operations to add the current mask 'mask'
    % element property to the history log item. This is rendered as 'Masks'
    % by the data tree, with additional information being mined from the
    % Masks block at the point of rendering.  This method parses the WMSK
    % block from the file and returns the current mask element.  If there
    % is no current mask the returned value will be empty.  In this case
    % the 'Masks' element should not be added to the analysis item history.  
    methods
        function Info = GetCurrentMaskInfo(this)
            Info = [];                  % Return empty array if no WMSK block.
            this.AssertSeek(0, 'bof');  % Go to beginning of file.
            % Find the WMSK block location and get data size.
            [blockOffset, ~, ~] = LocateBlock(this, 'WMSK');

            if (blockOffset ~= -1)
                % Offset from start of block to PSET size element = parsing
                % start is 20 bytes: 'WMSK' + 12 bytes + 'PSET'
                maskBlock = this.PSetParser(blockOffset + 20);

                % Find the the currentmask element and return it's contents.
                r = find(ismember(maskBlock(:,1), 'currentmask'));
                if isempty(r)
                    throw(WdfError('Current mask element not located in mask block.'));
                else
                    Info = maskBlock{r,2};
                end
            end
        end;
    end;
    
    %% Query for next available block UID
    methods
        % Get the UID of the next available block given its ID
        function [UID] = GetNextAvailableBlockUID(this, BlockID)
            % Setup search parameters
            offset = int64(512);
            highestFoundUID = 0;
            targetID = this.GetBlockID(BlockID);
            
            % Iteratively walk over all blocks in the file
            while (feof(this.Handle) == 0)
                % Attemp to jump to next block
                this.AssertSeek(offset, 'bof');
                
                % Read block header, and throw an error if block header
                % fields could not be read
                blockID = fread(this.Handle, 1, 'uint32');
                uid = fread(this.Handle, 1, 'uint32');
                bSize = fread(this.Handle, 1, 'uint64');
                if (isempty(blockID) || isempty(uid) || (numel(bSize) ~= 1) || (bSize < 16))
                    if (feof(this.Handle))
                        break;
                    else
                        throw(WdfError('Error whilst searching for next available block UID.'));
                    end;
                end;
                
                % Check if this block is of the target type ID, and if so,
                % record the UID if it is larger than any previously-found
                % block UID
                if (isequal(targetID, blockID))
                    if (uid > highestFoundUID)
                        highestFoundUID = uid;
                    end;
                end;
                
                % Increment the offset
                offset = offset + bSize;
            end;
            
            % We've searched through the whole file, so now we can safely
            % conclude that we have found the largest existing block UID,
            % and return the next available value.
            UID = int32(highestFoundUID + 1);
        end;
        
    end;
    
    %% Public static helper methods
    methods (Static, Access = public)
        % Converts a 4-character block-ID to its UINT32 value
        function [ID] = GetBlockID(String)
            try ID = uint32(sum(String .* (256 .^ (0:3))));
            catch
                ID = uint32(sum(cat(2,String,'P ') .* (256 .^ (0:3))));
            end;
        end;
        
        % Converts a 4-character block-ID and UINT32 UID to a 8-byte unique
        % block identifier
        function [FullID] = GetFullBlockID(BlockIDString, BlockUID)
            blockID = uint64(sum(BlockIDString .* (256 .^ (0:3))));
            FullID = int64(blockID + bitshift(uint64(BlockUID), 32));
        end;
    end
    
    %% Private static helper methods
    methods(Static, Access = private)
        % Returns a string consisting of concatenated meaningful terms
        % describing the maparea type corresponding to the WdfMapAreaFlags 
        % enum defined in wdf.h
        function Flag = ConvertFlagToString(flag)
            Flag = '';
            if flag(1) == true
                Flag = strcat(Flag,'RandomPoints ', ' ');
            end
            if flag(2) == true
                Flag = strcat(Flag, 'ColumnMajor ', ' ');
            else
                Flag = strcat(Flag , 'RowMajor ', ' ');
            end
            if flag(3) == true
                Flag = strcat(Flag , 'Snake ', ' ');
            else
                Flag = strcat(Flag , 'Raster ', ' ');
            end
            if flag(4) == true
               Flag = strcat(Flag , 'LineFocusMapping ');
            end
            if flag(7) == true
               Flag = strcat(Flag , 'SurfaceProfile ');
            end
            if flag(8) == true
               Flag = strcat(Flag , 'Line ');
            end
        end
        
        % Converts an integer decimal value to a vector representing its
        % value in the supplied Base.  The returned vector is organised 
        % least significant element first and padded to 8 elements.  
        function out = DecimalBaseConverter(Num,Base)
            out=[];
            i=1;
            while Num > 0
                out(i)=mod(Num,Base);
                Num = floor(Num/Base);
                i=i+1;
            end
            % Pad to length 8
            if(8-length(out) >= 0)
                out = [out,zeros(1,8-length(out))];
            else
                throw(WdfError('not enough elements to represent the data'));
            end
        end
    end

end