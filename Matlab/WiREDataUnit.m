% WiREDataUnit  Enumeration used by WiRE to indicate units of measurement

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

classdef WiREDataUnit < uint32
    enumeration
        Arbitrary          (0),
        RamanShift         (1),
        Wavenumber         (2),
        Nanometre          (3),
        ElectronVolt       (4),
        Micron             (5),
        Counts             (6),
        Electrons          (7),
        Millimetres        (8),
        Metres             (9),
        Kelvin            (10),
        Pascal            (11),
        Seconds           (12),
        Milliseconds      (13),
        Hours             (14),
        Days              (15),
        Pixels            (16),
        Intensity         (17),
        RelativeIntensity (18),
        Degrees           (19),
        Radians           (20),
        Celsius           (21),
        Fahrenheit        (22),
        KelvinPerMinute   (23),
        FileTime          (24),
        Microseconds      (25),
        Volts             (26),
        Amps              (27),
        MilliAmps         (28),
        Strain            (29),
        Ohms              (30),
        DegreesR          (31),
        Coulombs          (32),
        PicoCoulombs      (33)
    end;
end