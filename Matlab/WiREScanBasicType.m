% WiREScanBasicType  Enumeration used to indicate the basic WiRE scan type
%
% The scan-type is the method used to acquire a single spectrum / dataset.
% The basic scan type can be combined with various additional flags; see
% WiREScanType for details.
%
% The following basic scan types are used by WiRE:
%   Static         A single read-out off the detector where the grating is
%                  not moved.  Can be a spectrum or a CCD image.
%   Continuous     The grating is moved whilst data is read from the CCD,
%                  producing a spectrum over an extended range.
%   StepRepeat     Multiple overlapping static scans are taken and then
%                  'stitched' together to provide an extended range.
%   FilterScan   \ FilterScan and FilterImage are both provided for
%   FilterImage  / historical reasons and are unlikely to be encountered.
%   StreamLine     Renishaw's fast-mapping mode.
%   StreamLineHR   Renishaw's high confocality fast-mapping mode.
%   PointDetector  The scan is performed using a point detector.
%
% See also: WiREScanType

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

classdef WiREScanBasicType < uint32
    enumeration
        Unspecified        (0),
        Static             (1),
        Continuous         (2),
        StepRepeat         (3),
        FilterScan         (4),
        FilterImage        (5),
        StreamLine         (6),
        StreamLineHR       (7),
        PointDetector      (8)
    end;
end