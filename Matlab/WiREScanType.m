% WiREScanType  Represents a WiRE scan type value and associated flags
%
% The scan-type is the method used to acquire a single spectrum / dataset.
% The basic scan type (see WiREScanBasicType) can be combined with the
% following additional flags (which are exposed as read-only properties of
% a WiREScanType object):
%   IsMultitrackStitched  If TRUE, multiple tracks of data were acquired,
%                         but are 'stitched' together to produce a single
%                         spectrum / dataset per scan.
%   IsMultitrackDiscrete  If TRUE, multiple tracks of data were acquired
%                         and are stored separately using the multi-track
%                         file format (not currently supported in Matlab).
%   IsLineFocusMapping    If TRUE, line-focus mapping mode was used to
%                         acquire the data; certain adjacent map points
%                         were collected simultaneously.
%
% See also: WiREScanBasicType

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

classdef (Sealed = true) WiREScanType 
    properties (Constant, GetAccess = private)
        MultitrackStitched = uint32(256);
        MultitrackDiscrete = uint32(512);
        LineFocusMapping   = uint32(1024);
    end;
    
    properties (SetAccess = private)
        BasicType = WiREScanBasicType.Unspecified;
        IsMultitrackStitched = false;
        IsMultitrackDiscrete = false;
        IsLineFocusMapping = false;
    end;
    
    methods
        function v = WiREScanType(bitFlags)
            bitFlags = uint32(bitFlags);
            v.BasicType = WiREScanBasicType(bitand(uint32(255), bitFlags));
            v.IsMultitrackStitched = bitand(WiREScanType.MultitrackStitched, bitFlags) ~= 0;
            v.IsMultitrackDiscrete = bitand(WiREScanType.MultitrackDiscrete, bitFlags) ~= 0;
            v.IsLineFocusMapping = bitand(WiREScanType.LineFocusMapping, bitFlags) ~= 0;
        end;
    end;
end