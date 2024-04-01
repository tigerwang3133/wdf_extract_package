% Attach the wdf object to the WdfUtils object, but only allow set
% access, to prevent accessing all other methods within the wdf object
% from here.
%
% To use, set up a wdf object using 
% wdf = WdfReader(WdfName);
% Then set up a WdfUtils object as 
% Utils = WdfUtils;
% Set the wdf in the utils object to the required one as
% Utils.wdf = wdf;
% We can then call the methods in here using (for example)
% MeanSpectrum = Utils.GetMeanSpectrum(false)
% Which will generate the mean spectrum, without a waitbar.

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

classdef WdfUtils < handle
    
    properties (SetAccess = public, GetAccess = private)
        wdf = '';
    end;
    methods (Access = public)
        
        % Calculates the mean spectrum of the data. Use the argument to
        % optionally show a progress bar.
        function MeanSpec = GetMeanSpectrum(this, showProgress)
            if(showProgress)
                h = waitbar(0,'Calculating mean...');
            end
            MeanSpec = zeros(1,this.wdf.PointsPerSpectrum());
            TotalReadInSoFar = 0;
            this.wdf.StartChunkwiseReading();
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                Weight = (TotalReadInSoFar/(TotalReadInSoFar+NumberReadIn));
                MeanSpec = (1-Weight).*mean(z)+Weight.*MeanSpec;
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        % Calculates the mean spectrum of the data. Use the input to
        % optionally pass in the mean spectrum (if previously calculated)
        % or show a progress bar.
        function VarianceSpec = GetVarianceSpectrum(this, varargin)
            if(nargin == 1)
                MeanSpec = this.wdf.GetMeanSpectrum(false);
                showProgress = false;
            elseif(nargin == 2)
                MeanSpec = varargin{1};
                showProgress = false;
            elseif(nargin == 3)
                MeanSpec = varargin{1};
                showProgress = varargin{2};
            end
            if(showProgress)
                h = waitbar(0,'Calculating variance...');
            end
            VarianceSpec = zeros(1,this.wdf.PointsPerSpectrum());
            this.wdf.StartChunkwiseReading();
            TotalReadInSoFar = 0;
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                VarianceSpec = VarianceSpec + sum((z-MeanSpec).*(z-MeanSpec));
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        
        % Generates a map of the total intensity of each spectrum (sum of
        % the spectrum). Note this total intensity map is trivially related
        % to the mean value. 
        % Optionally show a progress bar by passing in true or false.
        function TotalIntensityMap = CalculateIntegratedIntensity(this, showProgress)
            TotalIntensityMap = zeros(this.wdf.Count(),1);
            TotalReadInSoFar = 1;
            if(showProgress)
                h = waitbar(0,'Calculating intensity map...');
            end
            this.wdf.StartChunkwiseReading();
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                TotalIntensityMap(TotalReadInSoFar:TotalReadInSoFar+NumberReadIn-1) = sum(z,2);
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        % Generates a map of the value of the spectrum at the nearest
        % wavenumber in the xlist to the chosen one.
        % Optionally show a progress bar by passing in true or false along
        % with the desired wavenumber.
        function TotalIntensityMap = PeakIntensityMap(this, Wavenumber, showProgress)
            xList = this.wdf.GetXList();
            Position = find(xList < Wavenumber, 1);
            TotalIntensityMap = zeros(this.wdf.Count(),1);
            TotalReadInSoFar = 1;
            if(showProgress)
                h = waitbar(0,'Calculating peak intensity map...');
            end
            this.wdf.StartChunkwiseReading();
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                TotalIntensityMap(TotalReadInSoFar:TotalReadInSoFar+NumberReadIn-1) = z(:,Position);
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        % MapIn is a 1D row vector consisting of the map values in
        % collected spectrum order as obtained from wdf.GetSpectra.  MapOut
        % is a 2D array of map values with x increasing l to r along the 
        % columns and y increasing down the rows.  Plotting image/sc(MapOut) 
        % produces an image in the same orientation as WiRE.
        % Currently supports 2D contiguous raster scans of all types from
        % inVia or RA802.  Does NOT support snake or linefocus maps.  
        
        % Steps to use ReorderMapForMatlabPlot:
        %         wdf = WdfReader(WdfName); % Read in wdf data.
        %         Utils = WdfUtils; % Set up a WdfUtils object
        %         Utils.wdf = wdf; % Set the wdf in the utils object
        %         [MapData, MapNames] = wdf.getMapData(); % Read in map data
        %         % MapNames contains the name of each map, use to find the index of the
        %         % required map (or loop through all of them)
        %         MapIndex = 1; % e.g find the first map
        %         MapOut= Utils.ReorderMapForMatlabPlot(MapData(MapIndex,:)); % Store the map in MapOut
        
        function MapOut = ReorderMapForMatlabPlot(this, MapIn)
            % Need the map to be passed in
            if(nargin < 2)
                throwAsCaller(WdfError('ReorderMapForMatlabPlot: The map data to reorder, MapIn, must be specified.'));
            end
            % Check that the size of the map passed in corresponds to the
            % number of spectra in the wdf  
            if (size(MapIn,2) ~= this.wdf.Count)
                throwAsCaller(WdfError('ReorderMapForMatlabPlot: mismatch between input map size and number of spectra in the WDF file.'));
            end
            % Check that the X and Y origin exists
            try
                wMap = this.wdf.GetWMapblock;
            catch
                throwAsCaller(WdfError('ReorderMapForMatlabPlot: Spatial information not found in WDF file.'));
            end
                        
            % Line type
            lineType = 0;
            fLength = length(wMap.flag);
            if strcmp(wMap.flag(fLength-3:fLength),'Line')
                if (wMap.StepSize(2) == 0)
                    lineType = 1;   % X Line
                elseif (wMap.StepSize(1) == 0)
                    lineType = 2;   % Y Line
                else 
                    lineType = 3;   % XY Line
                end
            end
            
            % Map dims from numPoints
            % For line maps numPoints(1) always contains the number of
            % points regardless of whether they are along X, Y or XY.  
            if (lineType == 2)
                numX = 1;
                numY = wMap.numPoints(1);
            else
                numX = wMap.numPoints(1);
                numY = wMap.numPoints(2);
            end
            
            % Check that the data collection is complete and spatial area 
            % not volume.  If Line then must be in X or Y direction.
            if ((numX*numY ~= this.wdf.Count) || (lineType == 3))
                throwAsCaller(WdfError('ReorderMapForMatlabPlot: only completed rectangular XY area, X and Y line mapping files are supported.'));
            end
                       
            % Setup collection order and direction flags
            xFirst = false; xDecreasing = false; yDecreasing = false; 
            if (((lineType == 0) && ~isempty(strfind(wMap.flag,'RowMajor'))) || lineType == 1)
                xFirst = true;  % RowMajor => x first, ColumnMajor => y first
            end
            if (wMap.StepSize(1) < 0)
                xDecreasing = true;
            end
            if (wMap.StepSize(2) < 0)
                yDecreasing = true;
            end
            
            % Reshape into a 2D array, depending on collection order   
            snakeFlag = ~isempty(strfind(wMap.flag,'Snake')); % Check for snake or linefocus snake
             if wMap.linefocus_size <= 1  % 0 if not LF, but also plots OK if size is 1.    
                if xFirst
                    MapOut = (reshape(MapIn,[numX, numY]))';
                else
                    MapOut = reshape(MapIn,[numY, numX]);
                end
                % If the map is collected in snake order, flip alternate rows
                if snakeFlag
                    MapOut(2:2:numY,:) = fliplr(MapOut(2:2:numY,:));
                end
            
            else % Linefocus, create map in multiple blocks   
                LF_size = wMap.linefocus_size;
                numBlocks = numY / LF_size; % Assume linefocus maps are always yFirst
                MapOut = [];
                for i = 1: numBlocks
                    MapBlock = MapIn((i-1)*numX*LF_size+1: i*numX*LF_size);
                    MapBlockOut = reshape(MapBlock,[LF_size, numX]); 
                    if snakeFlag && (rem(i,2)==0) % If linefocus snake, flip even blocks
                        MapBlockOut = fliplr(MapBlockOut);
                    end
                    MapOut = [MapOut; MapBlockOut];
                end
                % If the final block is incomplete, get the remaining data.
                % Note this is an unexpected scenario. 
                if numBlocks ~= floor(numBlocks) 
                    MapBlock = MapIn(numBlocks*numX*LF_size+1: end);
                    MapBlockOut = reshape(MapBlock,[LF_size*(numBlocks-floor(numBlocks)), numX]);
                    if snakeFlag && (rem(i+1,2)==0) % Flip if linefocus snake and even numbered block
                        MapBlockOut = fliplr(MapBlockOut);
                    end
                    MapOut = [MapOut; MapBlockOut];
                end
             end

            % Flip if x or y decreasing
            if xDecreasing
                MapOut = fliplr(MapOut);
            end
            if yDecreasing
                MapOut = flipud(MapOut);
            end
        end
        % End of public methods
    end
end