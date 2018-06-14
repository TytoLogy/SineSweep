%--------------------------------------------------------------------------
% NICal_NIexit.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
%
% closes NI devices nicely
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 9 July 2012 (SJS)
% 				Created from SpeakerCal_tdtexit.m
% 
% Revisions:
%	9 July, 2012 (SJS) renamed for NICal project
%	18 Jan 2017 (SJS): updated comments
%	1 Feb 2017 (SJS): added session interface things
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Clean up the IO things
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
disp('...closing NI devices...');

if handles.DAQSESSION
	% stop session
	stop(handles.iodev.NI.S);
	% release hardware
	release(handles.iodev.NI.S);
	% reset
	daqreset

else
	% get event log
	EventLogAI = showdaqevents(handles.iodev.NI.ai);
	EventLogAO = showdaqevents(handles.iodev.NI.ao);
	% delete and clear ai and ch0 object
	delete(handles.iodev.NI.ai);
	delete(handles.iodev.NI.ao);
	delete(handles.iodev.NI.chI);
	delete(handles.iodev.NI.chO);
	clear handles.iodev.NI.ai handles.iodev.NI.ao 
	clear handles.iodev.NI.chI handles.iodev.NI.chO
	% save settings information to mat file
	save(fullfile(pwd, 'NICal_EventLogs.mat'), ...
			'EventLogAI'			, ...
			'EventLogAO'			, ...
			'-MAT' );
end

	
