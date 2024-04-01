% WiREMeasurementType  Enumeration indicating WiRE measurement type
%
% There are three types of WiRE measurement:
%   Single    The file contains a single spectrum / dataset.
%   Series    The file contains a series of spectra / datasets, but these
%             spectra are related via something other than spatial position
%             (for example: time, temperature, pressure, etc).
%   Map       The file contains multiple spectra / datasets that relate to
%             different spatial positions in the sample.

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

classdef WiREMeasurementType < uint32
    enumeration
        Unspecified        (0),
        Single             (1),
        Series             (2),
        Map                (3)
    end;
end