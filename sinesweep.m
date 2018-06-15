% function varargout = sinesweep
%---------------------------------------------------------------------
% sinesweep.m
%---------------------------------------------------------------------
%{

Testing "notebook" for frequency response measurement with swept sinusoid

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
% desired sample rate
%----------------------------------------------------
desiredFs = 100000;
DevNum = 'Dev1';
% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
RMSsin = 1/sqrt(2);

%----------------------------------------------------
% file information
%----------------------------------------------------
testPath = pwd;

[dataName, dataPath] = uiputfile('*.mat', 'Sweep data output file', testPath);
if dataName == 0
	fprintf('Cancelled\n');
	return
else
	fprintf('Sweep data will be written to %s\n', fullfile(dataPath, dataName));
end

%----------------------------------------------------
% sweep parameters
%----------------------------------------------------
% properties of sweep (aka chirp) test signal
%	dur		duration (ms)
%	acq_dur	length of data to acquire (should be longer than sweep)	(ms)
%	start		start frequency (Hz)
%	end		end frequency (Hz)
%	mode		sweep mode: 'linear', 'log'
%	mag		peak level of output sweep (Volts)
%	ramp		ramp onset/offset duration (ms)
%	reps		# of times to present sweep
%	S			signal vector (empty for now)
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
sweep.end = 45000;
sweep.mode = 'log';
sweep.mag = 1;
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
%	mag		peak level of output sweep (Volts)
%	ramp		ramp onset/offset duration (ms)
%	reps		# of times to present sweep
tone.dur = 200;
tone.acq_dur = tone.dur + 10;
tone.freq = 1000*[10 16];
tone.mag	= 1;
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
dfilt.Fc_lo = 47000;
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
fprintf('\tTo process data after collection, use the process_sweep.m script.\n');

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
for f = 1:length(tone.freq)
	% generate raw signal
	tmp = syn_calibrationtone2(tone.dur, iodev.Fs, tone.freq(f), 0, 'L');
	% save null stimulus
	if f == 1
		nullstim = tmp(2, :);
	end
	tone.S{f, 1} = tmp(1, :);
	% scale the sound
	tone.S{f, 1} = tone.mag * tone.S{f, 1};
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
figure(1)
% generate time vector for tone plots
tt = 0:dt:( (0.001*tone.dur) - dt);
% loop through tones
fprintf('Playing calibration tones to determine SPL level...\n');
for f = 1:length(tone.freq)
	% build output matrix
	%	first column goes to channel 1, second column to channel 2, so need to
	%	append two row vectors (nidaq_session_io will do the transposition)
	outstim = [tone.S{f, 1}; nullstim];
	% then play stimulus, record response
	fprintf('Playing %d kHz tone...\n', tone.freq(f));
	for n = 1:tone.reps
		fprintf('\t\trep # %d ...', n);
		[tmp, ~] = nidaq_session_io(iodev, outstim, 0.001*tone.acq_dur);
		% keep only the channel 0 (left) data in tmp{1};
		tone.R{f, n} = tmp{1};
		clear tmp;
		% calculate db SPL
		tmp = sin2array(tone.R{f, n}, 1, iodev.Fs);
		tmp = filtfilt(dfilt.b, dfilt.a, tmp);
		% determine the magnitude and phase of the response
		mag = fitsinvec(tmp, 1, iodev.Fs, tone.freq(f));
		% compute dB SPL, using precomputed conversion factor
		magdB = dbspl(RMSsin * mic.VtoPa * mag);
		% report results
		plot(tt, tmp);
		tstr = sprintf('Level = %.2f dB SPL\t(%.4f mV)', magdB, 1000*mag);
		title({sprintf('%d kHz, Rep %d', tone.freq(f), n), tstr});
		xlabel('time (ms)')
		ylabel('V');
		fprintf('%s\n', tstr);
		% pause for a bit
		pause(0.5);
	end
	fprintf('...done\n');
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% sweep input/output
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% build output matrix
%	first column goes to channel 1, second column to channel 2, so need to
%	append two row vectors (nidaq_session_io will do the transposition)
outstim = [sweep.S; 0*sweep.S];
fprintf('flushing buffer\n');
[tmp, ~] = nidaq_session_io(iodev, [0*sweep.S; 0*sweep.S], 0.001*sweep.acq_dur);

% then play stimulus, record response
for n = 1:sweep.reps
	fprintf('Playing sweep, rep # %d ...', n);
	[tmp, ~] = nidaq_session_io(iodev, outstim, 0.001*sweep.acq_dur);
	% keep only the channel 0 (left) data in tmp{1};
	sweep.R{n, 1} = tmp{1};
	fprintf('...done\n');
	% pause for a bit
	pause(0.5);
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
% plot input (stim) ...
figure(2)
subplot(211)
plot(t, sweep.S);
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Sweep Signal');
xlabel('time (ms)')
% plot response (resp)
figure(3)
subplot(211)
plot(t, R{1, 1});
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Sweep Response')
xlabel('time (ms)')

% process data, using length as NFFT
% first, make sure length is even
if ~even(length(sweep.S))
	% if not, pad with a 0
	sweep.S = [sweep.S 0];
end
if ~even(length(R{1, 1}))
	% if not, pad with a 0
	R{1, 1} = [R{1, 1} 0];
end
sweep.Nfft = length(sweep.S);
sweep.Sfft = fft(sweep.S, sweep.Nfft);
resp.Nfft = length(R{1, 1});
resp.Sfft = fft(R{1, 1}, resp.Nfft);

% compute magnitude spectrum
% non-redundant points are kept
Nunique = (sweep.Nfft/2) + 1;
sweep.mag = sweep.Sfft(1:Nunique);
% and magnitude computed
sweep.mag = 2*abs(sweep.mag)/sweep.Nfft;
sweep.phi = angle(sweep.Sfft(1:Nunique));

% non-redundant points are kept
Nunique = (resp.Nfft/2) + 1;
resp.mag = resp.Sfft(1:Nunique);
% and magnitude computed
resp.mag = 2*abs(resp.mag)/resp.Nfft;
resp.phi = angle(resp.Sfft(1:Nunique));

figure(2)
subplot(212)
Freq = Fnyq*linspace(0, 1, length(sweep.mag));
semilogx(Freq, db(sweep.mag));
xlim([0 Fnyq]);
ylim([-100 0]);
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');
figure(3)
subplot(212)
Freq = Fnyq*linspace(0, 1, length(resp.mag));
semilogx(Freq, db(resp.mag));
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
save(fullfile(dataPath, dataName), 'sweep', 'tone', 'dfilt', 'mic', '-MAT')

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Close NIDAQ device
%---------------------------------------------------------------------
%---------------------------------------------------------------------
[iodev, status] = nidaq_close(iodev);

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% assign outputs
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% if nargout
% 	if any(nargout == [1 2 3])
% 		varargout{1} = sweep;
% 	end
% 	if any(nargout == [2 3])
% 		varargout{2} = iodev;
% 	end
% 	if nargout == 3
% 		varargout{3} = status;
% 	end
% end
% 
% 



