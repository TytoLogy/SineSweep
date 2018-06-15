function [iodev, status] = nidaq_close(iodev)
%--------------------------------------------------------------------------
% nidaq_close.m
%--------------------------------------------------------------------------
%  ->  -> nidaq tools
%--------------------------------------------------------------------------
% NI data acquisition toolbox parameters
% closes NI devices nicely
%------------------------------------------------------------------------
% Input Arguments:
% 	iodev		initialized NI struct:
%						NI.S		session interface
%						NI.chO	analog output channel objects
%						NI.chI	analog input channel objects
% Output Arguments:
% 	iodev		closed NI struct:
%						NI.S		session interface
%						NI.chO	analog output channel objects
%						NI.chI	analog input channel objects
%	status	0 if closed successfully
% 				1 if error
%------------------------------------------------------------------------
% See also: nidaq_open, nidaq_aiao_init, nidaq_session_io
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 15 June 208 (SJS)
% 				Created from NICAL_NIexit.m
% 
% Revisions:
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Clean up the IO things
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
disp('...closing NI devices...');

try
	% stop session
	stop(iodev.NI.S);
	% release hardware
	release(iodev.NI.S);
	% reset
	daqreset
catch errMsg
	errordlg('error closing NI device')
	fprintf('%s: %s\n%s\n', mfilename, errMsg.identifier, errMsg.message);
	status = 1;
	return
end
disp('...closed successfully')
status = 0;


