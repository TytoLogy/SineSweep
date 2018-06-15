function [iodev, init_status] = nidaq_open(device_number, sample_rate)
%--------------------------------------------------------------------------
% nidaq_open.m
%--------------------------------------------------------------------------
%  ->  -> nidaq tools
%--------------------------------------------------------------------------
% NI data acquisition toolbox parameters
% sets up nidaq system 
%------------------------------------------------------------------------
% Input Arguments:
% 	device_number	 usually 1 (identifies NI-DAQmx hardware)
%	sample_rate		 sample rate
% Output Arguments:
% 	iodev		initialized NI struct:
%						NI.S		session interface
%						NI.chO	analog output channel objects
%						NI.chI	analog input channel objects
% init_status	0 if unsuccessful (handles.NI will be empty)
% 					1 if successful
%------------------------------------------------------------------------
% See also: nidaq_close, nidaq_aiao_init, nidaq_session_io
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 14 June 2018 (SJS)
% 				Created from NICal -> NICal_NIinit.m
% 
% Revisions:
%--------------------------------------------------------------------------

fprintf('%s: starting NI hardware...\n', mfilename);

%-----------------------------------------------------------------------
% Settings/Constants
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% Initialize the NI device
%-----------------------------------------------------------------------
try
	iodev.NI = nidaq_aiao_init('NI-SESSION', device_number);
catch errMsg
	errordlg('error initializing NI device')
	fprintf('%s: %s\n%s\n', mfilename, errMsg.identifier, errMsg.message);
	init_status = 0;
	return
end
%-----------------------------------------------------------------------
% set sample rate to value specified in cal settings
%-----------------------------------------------------------------------
iodev.NI.S.Rate = sample_rate;
% check actual rate for mismatch
ActualRate = iodev.NI.S.Rate;
if sample_rate ~= ActualRate
	warning('nidaq_init:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', ...
					sample_rate, ActualRate);
end
% store actual rate
iodev.Fs = sample_rate;

%-----------------------------------------------------------------------
% input, output channel properties
%-----------------------------------------------------------------------
% range needs to be in [RangeMin RangeMax] format
aiaoRange = 5 * [-1 1];
for n = 1:length(iodev.NI.chI)
	% set analog input range
	iodev.NI.chI(n).Range = aiaoRange;
	% set input TerminalConfig to 'SingleEnded' (default is
	% 'Differential')
	iodev.NI.chI(n).TerminalConfig = 'SingleEnded';
end
% set analog output range
for n = 1:length(iodev.NI.chO)
	iodev.NI.ao.chO(n).Range = aiaoRange;
end

%------------------------------------------------------------------------
% HARDWARE TRIGGERING
%------------------------------------------------------------------------
% only 1 "sweep" per trigger event 
iodev.NI.S.TriggersPerRun = 1;
%-------------------------------------------------------
% set init_status to 1
%-------------------------------------------------------
init_status = 1;
%-------------------------------------------------------
% Done!!!!
%-------------------------------------------------------

fprintf('done\n');
