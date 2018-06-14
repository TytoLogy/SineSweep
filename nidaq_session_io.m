function [resp, varargout] = nidaq_session_io(iodev, stim, sweepduration)
%--------------------------------------------------------------------------
% [resp, index, timestamps, triggertime] = 
% 										nidaq_session_io(iodev, stim, sweepduration)
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% input/output for calibration using SESSION interface
%------------------------------------------------------------------------
% Input Arguments:
%	iodev					input/output struct
%	stim					stimulus vector
%	sweepduration		time to to collect data (seconds)
% 
% Output Arguments:
%	resp				collected data {2X1} cell with vectors (1Xinpts) in size
%	index				# points collected
%	timestamps		sample timestamps
%	triggertime		DAQ trigger timing
%------------------------------------------------------------------------
% See also: NICal, nidaq_aiao_init, DAQ Toolbox (MATLAB)
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 1 February 2017 from nidaq_calibration_io (SJS)
% 
% Revisions:
%--------------------------------------------------------------------------

% load stimulus onto NI memory
queueOutputData(iodev.NI.S, stim');
% START ACQUIRING
[rawdata, timestamps, triggertime] = startForeground(iodev.NI.S);
% stop acquiring
stop(iodev.NI.S);
% reformat output into cell (legacy)
resp{1} = rawdata(:, 1)';
resp{2} = rawdata(:, 2)';
% assign outputs by determining how many outputs were requested
nout = max(nargout,1) - 1;
if nout >= 1
	% return "index"
	varargout{1} = length(rawdata(:, 1));
end
if nout >= 2
	varargout{2} = timestamps;
end
if nout >= 3
	varargout{3} = triggertime;
end

