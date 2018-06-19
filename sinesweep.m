function varargout = sinesweep
%---------------------------------------------------------------------
% sinesweep.m
%---------------------------------------------------------------------
%{

Function for frequency response measurement with swept sinusoid
requires Data acquisition toolbox, TytoLogy toolboxes, NI-DAQ drivers,
National Instruments hardware

Process:

(1) User confirms settings
(2) User specifies output file
(3) Initialize NIDAQ interface
(4) Generate sweep 
(5) Play sweep
(6) Play tones to determine SPL scaling (may want to do this before the
sweep so that output level can be examined to see if it is too loud)
(7) Write data to file (or do this after each signal I/O event)
(8) Close NIDAQ interface

%}
%---------------------------------------------------------------------


%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: June 2018 (SJS)
% 
% Revisions:
%	19 Jun, 2018 (SJS):
% 		- created as function (script is sinesweep_script)
% 		- tested, seems to be working
%							
%--------------------------------------------------------------------------

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% check paths for necessary toolboxen
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% directory when using installed version:
pdir = ['C:\TytoLogy\Toolboxes\TytoLogySettings\' getenv('USERNAME')];
if isempty(which('ms2samples'))
	run(fullfile(pdir, 'tytopaths'))
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Settings
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%----------------------------------------------------
% I/O settings
%----------------------------------------------------
% input, output sampling rate (samples/second)
desiredFs = 500000;
% inter-stimulus interval (milliseconds)
ISI = 500;
% NIDAQ device number (from daq.getDevices)
DevNum = 'Dev1';
%----------------------------------------------------
% Constants, precomputed
%----------------------------------------------------
% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
RMSsin = 1/sqrt(2);
%----------------------------------------------------
% file information
%----------------------------------------------------
% set default data directory
if ispc
	if exist('E:\', 'dir')
		defaultPath = 'E:\Calibrate\SweepData';
		if ~exist(defaultPath, 'dir')
			defaultPath = 'E:\';
		end
	elseif exist('F:\', 'dir')
		defaultPath = fullfile('F:\Users', username);
	else
		defaultPath = fullfile('C:\Users', username);
	end
else
	defaultPath = ['~' username];
end
% get data file and directory
[dataName, dataPath] = uiputfile('*.mat', 'Sweep data output file', defaultPath);
if dataName == 0
	fprintf('Cancelled\n');
	return
else
	fprintf('Sweep data will be written to %s\n', fullfile(dataPath, dataName));
end

%----------------------------------------------------
%----------------------------------------------------
% Signal/Stimulus/Microphone Settings
%----------------------------------------------------
%----------------------------------------------------

%----------------------------------------------------
% sweep parameters
%----------------------------------------------------
% properties of sweep (aka chirp) test signal
%	dur		duration (ms)
%	acq_dur	length of data to acquire (should be longer than sweep)	(ms)
%	start		start frequency (Hz)
%	end		end frequency (Hz)
%	mode		sweep mode: 'linear', 'log'
%	peak		peak level of output sweep (Volts)
%	ramp		ramp onset/offset duration (ms)
%	reps		# of times to present sweep
%	S			signal vector (empty for now)
%	R			{1, reps} cell array of responses to S
%	Fs			sample rate (samples/s) - determined after NIDAQ init
%----------------------------------------------------
% Some things to note for use of this in frequency response measurements:
%	- sweep frequencies should start slightly below and end slightly above
%	the desired characterization range (if device-under-test supports or
%	will not be damaged at these frequency limits!!!)
%	- Note that even though there is a reps setting, ONLY THE FIRST REP will
%	be processed/used in this script!!!!!
%----------------------------------------------------
sweep.dur = 1000;
sweep.acq_dur = sweep.dur + 10;
sweep.start = 3500;
sweep.end = 50000;
sweep.mode = 'log';
sweep.peak = 1;
sweep.ramp = 0.1;
sweep.reps = 1;
sweep.S = [];
sweep.R = cell(1, sweep.reps);
sweep.Fs = [];
%----------------------------------------------------
% test tone parameters
%----------------------------------------------------
% tones are used to calibrate the levels in dB SPL
%	dur		duration (ms)
%	acq_dur	length of data to acquire (should be longer than sweep)	(ms)
%	freq		frequency or frequencies to test (Hz)
%	peak		peak level of output sweep (Volts)
%	ramp		ramp onset/offset duration (ms)
%	reps		# of times to present sweep
%	S			{# frequencies, 1} cell array of tone signals (empty for now)
%	R			{# frequencies, reps} cell array of responses to S
%	Fs			sample rate (samples/s) - determined after NIDAQ init
tone.dur = 200;
tone.acq_dur = tone.dur + 10;
tone.freq = 1000*[10 16 32];
tone.peak = 1;
tone.ramp = 1;
tone.reps = 3;
tone.S = cell(length(tone.freq), 1);
tone.R = cell(length(tone.freq), tone.reps);
tone.Fs = [];
%----------------------------------------------------
% these are parameters for calibration microphone
%----------------------------------------------------
%	gain		mic gain (dB)
%	sense		mic sensitivity (Volts/Pascal)
%	VtoPa		conversion factor
mic.gain = 0;
mic.sense = 1;
mic.VtoPa = (1/invdb(mic.gain)) * (1 / mic.sense);
%----------------------------------------------------
% parameters for input data filtering
%----------------------------------------------------
% 	Fc_lo		lowpass filter cutoff frequency (high frequency limit)
% 	Fc_hi		highpass filter cutoff frequency (low frequency limit)
% 	order		filter order
% 	b, a		filter coefficients (empty for now)
dfilt.Fc_lo = 55000;
dfilt.Fc_hi = 1500;
dfilt.order = 3;
dfilt.b = [];
dfilt.a = [];

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Tell user whats happening
%---------------------------------------------------------------------
%---------------------------------------------------------------------
fprintf('Running sinesweep to measure frequency response\n')
fprintf('\tTo process data after collection, use process_sweep\n');

if query_user('Continue') == 0
	fprintf('sinesweep Exiting\n');
	return
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Initialize NIDAQ
%---------------------------------------------------------------------
%---------------------------------------------------------------------
[iodev, status] = nidaq_open(DevNum, desiredFs);
if ~status
	error('Bad return status for nidaq_open');
end
fprintf('Initialized nidaq with Fs = %.2f\n', iodev.Fs)
% Nyquist freq and time interval
Fnyq = iodev.Fs / 2;
dt = 1/iodev.Fs;
% save sample rate in signal structs
sweep.Fs = iodev.Fs;
tone.Fs = iodev.Fs;

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% build some things
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%----------------------------------------------------
% Define a bandpass filter for processing the data
%----------------------------------------------------
% passband definition
dfilt.band = [dfilt.Fc_hi dfilt.Fc_lo] ./ Fnyq;
% filter coefficients using a butterworth bandpass filter
[dfilt.b, dfilt.a] = butter(dfilt.order, dfilt.band, 'bandpass');

%----------------------------------------------------
% synthesize sweep
%----------------------------------------------------
% generate time vector
t = 0:dt:( (0.001*sweep.dur) - dt);
% generate raw signal
sweep.S = chirp(t, sweep.start, sweep.dur/1000, sweep.end, 'logarithmic');
% ramp on/off if requested
if sweep.ramp > 0
	sweep.S = sin2array(sweep.S, sweep.ramp, iodev.Fs);
end

%----------------------------------------------------
% synthesize test tones
%----------------------------------------------------
% loop through test frequencies
for f = 1:length(tone.freq)
	% generate raw signal for output on Channel 0 ("left" or 'L')
	tmp = syn_calibrationtone2(tone.dur, iodev.Fs, tone.freq(f), 0, 'L');
	% save null stimulus (nidaq i/o expects two channels of output data)
	if f == 1
		nullstim = tmp(2, :);
	end
	tone.S{f, 1} = tmp(1, :);
	% scale the sound
	tone.S{f, 1} = tone.peak * tone.S{f, 1};
	if tone.ramp > 0
		% apply the sin^2 amplitude envelope to the stimulus
		tone.S{f, 1} = sin2array(tone.S{f, 1}, tone.ramp, iodev.Fs);
	end
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% tone input/output
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%----------------------------------------------------
% setup plotting things
%----------------------------------------------------
figure(1)
% generate time vector for tone plots
tt = 0:dt:( (0.001*tone.dur) - dt);
%----------------------------------------------------
% loop through tones and play/record
%----------------------------------------------------
fprintf('Playing calibration tones to determine SPL level...\n');
for f = 1:length(tone.freq)
	% build output matrix
	%	first column goes to channel 1, second column to channel 2, so need to
	%	append two row vectors (nidaq_session_io will do the transposition)
	outstim = [tone.S{f, 1}; nullstim];
	% then play stimulus, record response
	fprintf('Playing %d kHz tone, %.2f V peak...\n', tone.freq(f), tone.peak);
	for n = 1:tone.reps
		fprintf('\t\trep # %d ...', n);
		[tmp, ~] = nidaq_session_io(iodev, outstim, 0.001*tone.acq_dur);
		% keep only the channel 0 (left) data in tmp{1};
		tone.R{f, n} = tmp{1};
		clear tmp;
		% window, filter input data
		tmp = sin2array(tone.R{f, n}, 1, iodev.Fs);
		tmp = filtfilt(dfilt.b, dfilt.a, tmp);
		% determine the magnitude and phase of the response
		mag = fitsinvec(tmp, 1, iodev.Fs, tone.freq(f));
		% compute dB SPL, using precomputed conversion factor
		magdB = dbspl(RMSsin * mic.VtoPa * mag);
		% report results
		plot(tt, tmp);
		tstr1 = sprintf('%d kHz, Rep %d', 0.001*tone.freq(f), n);
		tstr2 = sprintf('Level = %.2f dB SPL\t(%.4f mV)', magdB, 1000*mag);
		title({tstr1, tstr2});
		xlabel('time (ms)')
		ylabel('V');
		fprintf('%s\t%s\n', tstr1, tstr2);
		% pause for ISI
		pause(0.001*ISI);
	end
	fprintf('...done\n');
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% sweep input/output
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%----------------------------------------------------
% build output matrix
%----------------------------------------------------
% first column goes to channel 1, second column to channel 2, so need to
% append two row vectors (nidaq_session_io will do the transposition)
outstim = [sweep.S; 0*sweep.S];
%----------------------------------------------------
% Odd things were happening with contamination of 
% sweep with tone stimulus... kludgey fix is to 
% play 0s..... ???????? 
%----------------------------------------------------
fprintf('flushing I/O buffer\n');
[~, ~] = nidaq_session_io(iodev, [0*sweep.S; 0*sweep.S], 0.001*sweep.acq_dur);
% pause for ISI
pause(0.001*ISI);
%----------------------------------------------------
% then play sweep stimulus, record response
%----------------------------------------------------
for n = 1:sweep.reps
	fprintf('Playing sweep, rep # %d ...', n);
	[tmp, ~] = nidaq_session_io(iodev, outstim, 0.001*sweep.acq_dur);
	% keep only the channel 0 (left) data in tmp{1};
	sweep.R{n, 1} = tmp{1};
	fprintf('...done\n');
	% pause for a ISI
	pause(0.001*ISI);
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% plot sweep data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%----------------------------------------------------
% process sweeps using filtfilt to remove noise
%----------------------------------------------------
R = sweep.R;
for n = 1:sweep.reps
	R{1, n} = filtfilt(dfilt.b, dfilt.a, sweep.R{1, n});
end

%----------------------------------------------------
% plot signals
%----------------------------------------------------
% local copy of sweep signal
S = sweep.S;
% plot input (stim) ...
figure(2)
subplot(211)
plot(t, S);
xlim([min(t) max(t)])
ylim(sweep.peak*[-1.1 1.1])
title('Sweep Signal');
xlabel('time (ms)');
% plot response (resp)
figure(3)
subplot(211)
plot(t, R{1, 1});
xlim([min(t) max(t)]);
ylim(max(abs(R{1, 1}))*[-1.1 1.1]);
title('Sweep Response');
xlabel('time (ms)');

% process data, using length as NFFT
% first, make sure length is even so that NFFT/2 is an integer
if ~even(length(sweep.S))
	% if not, pad with a 0
	S = [sweep.S 0];
end
if ~even(length(R{1, 1}))
	% if not, pad with a 0
	R{1, 1} = [R{1, 1} 0];
end
Nfft_S = length(S);
Sfft = fft(S, Nfft_S);
Nfft_R = length(R{1, 1});
Rfft = fft(R{1, 1}, Nfft_R);

% compute magnitude spectrum
% non-redundant points are kept
Nunique_S = (Nfft_S/2) + 1;
Smag = Sfft(1:Nunique_S);
% and magnitude computed
Smag = 2*abs(Smag)/Nfft_S;
Sphi = angle(Sfft(1:Nunique_S)); %#ok<NASGU>

% non-redundant points are kept
Nunique_R = (Nfft_R/2) + 1;
Rmag = Rfft(1:Nunique_R);
% and magnitude computed
Rmag = 2*abs(Rmag)/Nfft_R;
Rphi = angle(Rfft(1:Nunique_R)); %#ok<NASGU>

figure(2)
subplot(212)
Freq = Fnyq*linspace(0, 1, length(Smag));
semilogx(Freq, db(Smag));
xlim([0 Fnyq]);
ylim([-100 0]);
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');
figure(3)
subplot(212)
Freq = Fnyq*linspace(0, 1, length(Rmag));
semilogx(Freq, db(Rmag));
xlim([0 Fnyq]);
ylim([-100 0]);
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% save data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
fprintf('Saving data to %s\n', fullfile(dataPath, dataName));
save(fullfile(dataPath, dataName), 'sweep', 'tone', 'dfilt', 'mic', '-MAT')

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Close NIDAQ device
%---------------------------------------------------------------------
%---------------------------------------------------------------------
fprintf('Closing NI-DAQ interface...\n');
[iodev, status] = nidaq_close(iodev);

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% assign outputs
%---------------------------------------------------------------------
%---------------------------------------------------------------------
if nargout
	if any(nargout == [1 2 3])
		varargout{1} = sweep;
	end
	if any(nargout == [2 3])
		varargout{2} = iodev;
	end
	if nargout == 3
		varargout{3} = status;
	end
end





