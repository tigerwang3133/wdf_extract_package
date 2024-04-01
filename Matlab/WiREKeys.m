% WiREKeys  Enumeration used by WiRE to indicate the p-set type.

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

classdef WiREKeys < uint32
    enumeration
        WiRE2Version (1),
        ErrorCode (2),
        CreationTime (3),
        StartTime (4),
        EndTime (5),
        ETA (6),
        wizardclsid (7),
        ExpectedDatasetCount (8),
        AcquisitionCount (9),
        NumberOfPoints (10),
        FileHandlerVersion (11),
        autoSaveInterval (12),
        MeasurementType (13),
        responseCalibration (14),
        restoreInstrumentState (15),
        closeLaserShutter (16),
        FocusTrackEnabled (17),
        FocusTrackInterval (18),
        LineFocusMode (19),
        yStepSize (20),
        DepthSeriesInterval (201),
        DepthSeriesStartPos (202),
        DepthSeriesFinalPos (203),
        ScanCompletionStatus (301),
        usingPixelIntensityVariationFunction (302),
        Results (401),
        Property (402),
        Properties( 403),
        Data(410),
        Label (411),
        MapType (412),
        DataList0 (420),
        DataList1 (421),
        DataList2 (422),
        DataList3 (423),
        DataList4 (424),
        Operator (430),
        Time (431),
        Version (432),
        Name (1001),
        Description (1002),
        Language (1003),
        Code (1004),
        Status (1005),
        System (1006),
        DataFileFormat (1007),
        DataSaveMode (1008),
        DataSaveFile (1009),
        Result (1010),
        NamedItems (1011),
        Image (1012),
        AreaKey (1013),
        InstrumentState (1014),
        LaserName (1015),
        BeamPath (1016),
        FocusMode (1017),
        GratingName (1018),
        Laser (1019),
        Camera (1020),
        Instrument (1021),
        Objectivemagnification (1022),
        Objectivename (1023),
    end;
end