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
%% Settings
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%----------------------------------------------------
% desired sample rate
%----------------------------------------------------
desiredFs = 100000;
DevNum = 1;

%----------------------------------------------------

%----------------------------------------------------
dataPath = pwd;
dataName = 'testdata';

%----------------------------------------------------
% sweep parameters
%----------------------------------------------------
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
%----------------------------------------------------
sweep.dur = 1000;
sweep.acq_dur = sweep.dur + 10;
acq_dur = sweep.dur + 20;
sweep.start = 4000;
sweep.end = 45000;
sweep.mode = 'log';
sweep.mag = 1;
sweep.ramp = 0.1;
sweep.reps = 1;
sweep.S = [];
sweep.R = cell(1, sweep.reps);

%----------------------------------------------------
% test tone parameters
%----------------------------------------------------
% tones are used to calibrate the levels in dB SPL
%	dur		duration (ms)
%	acq_dur	length of data to acquire (should be longer than sweep)	(ms)
%	freq		frequency or frequencies to test (kHz)
%	mag		peak level of output sweep (Volts)
%	ramp		ramp onset/offset duration (ms)
%	reps		# of times to present sweep
tone.dur = 200;
tone.acq_dur = tone.dur + 10;
tone.freq = [10 16];
tone.mag	= 1;
tone.ramp = 1;
tone.reps = 3;
tone.S = cell(length(tone.freq), 1);
tone.R = cell(length(tone.freq), tone.reps);

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
dfilt.Fc_hi = 3500;
dfilt.order = 3;
dfilt.b = [];
dfilt.a = [];

%----------------------------------------------------
% other things
%----------------------------------------------------


%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Initialize NIDAQ
%---------------------------------------------------------------------
%---------------------------------------------------------------------
[iodev, status] = nidaq_init(DevNum, desiredFs);
if ~status
	error('Bad return status for nidaq_init');
else
	fprintf('Initialized nidaq with Fs = %.2f\n', iodev.Fs)
end
% Nyquist freq and time interval
Fnyq = iodev.Fs / 2;
dt = 1/iodev.Fs;

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
[dfilt.b, dfilt.a] = butter(dfilt.forder, dfilt.fband, 'bandpass');

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


%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% sweep input/output
%---------------------------------------------------------------------
%---------------------------------------------------------------------
for n = 1:sweep.reps
	fprintf('Playing sweep, rep # %d ...', n);
	[sweep.R{1, n}, indx] = nidaq_session_io(iodev, sweep.S, sweep.acq_dur);
	fprintf('...done\n');
	% pause for a bit
	pause(0.5);
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% process sweep data
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%----------------------------------------------------
% process sweeps using filtfilt
%----------------------------------------------------
for n = 1:sweep.reps
	sweep.R{1, n} = filtfilt(dfilt.b, dfilt.a, sweep.R{1, n});
end



*******************start work here!~~~~~~

%----------------------------------------------------
% plot signals
%----------------------------------------------------
% plot input (stim) ...
figure(1)
subplot(211)
plot(t, sweep.S);
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Input');
xlabel('time (ms)')
% plot response (resp)
figure(2)
subplot(211)
plot(t, sweep.R);
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Output')
xlabel('time (ms)')

% process data, using length as NFFT
% first, make sure length is even
if ~even(length(sweep.S))
	% if not, pad with a 0
	sweep.S = [sweep.S 0];
end
if ~even(length(sweep.R))
	% if not, pad with a 0
	sweep.R = [sweep.R 0];
end
sweep.Nfft = length(sweep.S);
sweep.Sfft = fft(sweep.S, sweep.Nfft);
resp.Nfft = length(sweep.R);
resp.Sfft = fft(sweep.R, resp.Nfft);

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

figure(1)
subplot(212)
Freq = Fnyq*linspace(0, 1, length(sweep.mag));
semilogx(Freq, db(sweep.mag));
xlim([0 Fnyq]);
ylim([-100 0]);
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');
figure(2)
subplot(212)
Freq = Fnyq*linspace(0, 1, length(resp.mag));
semilogx(Freq, db(resp.mag));
xlim([0 Fnyq]);
ylim([-100 0]);
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');


%% build compensation for -3dB/octave gain characteristic of sweep
figure(3)

% correction for -3dB/octave (or -10dB/decade)
correct3dB = 0.5*db(Freq);
% set any ±Inf values to 0
correct3dB(isinf(correct3dB)) = 0;
% apply correction only to appropriate region of sweep
fbins = find((Freq >= sweep.start) & (Freq <= sweep.end));
sweep.mag_corr3dB = sweep.mag;
sweep.mag_corr3dB(fbins) = sweep.mag(fbins) .* invdb(correct3dB(fbins));

% plot raw, correction factor, corrected sweep spectrum
semilogx(Freq, db(sweep.mag), '.-');
hold on
	semilogx(Freq, correct3dB, '.-');
	semilogx(Freq, db(sweep.mag) + correct3dB , '.-');
	semilogx(Freq, db(sweep.mag_corr3dB), '.-');
hold off
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');

%%



