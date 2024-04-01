function Exception = WdfError(Msg, varargin)

% WdfError  Initialise an MException with "WdfError" message-ID
%
% WE = WdfError(ERRMSG, V1, V2, ... VN) captures information about a
% specific error and stores it in WdfError object WE, which has the
% message-ID "Renishaw:SPD:WiRE:WdfError".
%
% See also: MException

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

Exception = MException('Renishaw:SPD:WiRE:WdfError', Msg, varargin{:});