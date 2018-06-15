%---------------------------------------------------------------------
% process_sweep.m
%---------------------------------------------------------------------
%{

process data for frequency response measurement with swept sinusoid

Process:

(1) load data from .mat file

%}
%---------------------------------------------------------------------

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Settings
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%----------------------------------------------------
% file information
%----------------------------------------------------
testPath = pwd;

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% load data from file
%---------------------------------------------------------------------
%---------------------------------------------------------------------
[dataName, dataPath] = uigetfile('*.mat', 'Select sweep data file', testPath);
if dataName == 0
	fprintf('Cancelled\n');
	return
else
	fprintf('Sweep data will be read from %s\n', fullfile(dataPath, dataName));
end
load(fullfile(dataPath, dataName))

%----------------------------------------------------
% Information in mat file from sinesweep program
%----------------------------------------------------
% sweep struct:
% 				properties of sweep (aka chirp) test signal
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
% tone struct:
% 				test tone parameters
%----------------------------------------------------
% tones are used to calibrate the levels in dB SPL
%	dur		duration (ms)
%	acq_dur	length of data to acquire (should be longer than sweep)	(ms)
%	freq		frequency or frequencies to test (kHz)
%	mag		peak level of output sweep (Volts)
%	ramp		ramp onset/offset duration (ms)
%	reps		# of times to present sweep
%	S			cell array of tone stimulis {length(freq), 1}
%----------------------------------------------------
% mic struct:
% 				these are parameters for calibration microphone
%----------------------------------------------------
%	gain		mic gain (dB)
%	sense		mic sensitivity (Volts/Pascal)
%	VtoPa		conversion factor
%----------------------------------------------------
% dfilt struct:
%				parameters for input data filtering
%----------------------------------------------------
% 	Fc_lo		lowpass filter cutoff frequency (high frequency limit)
% 	Fc_hi		highpass filter cutoff frequency (low frequency limit)
% 	order		filter order
% 	b, a		filter coefficients (empty for now)

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% process sweep data
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%----------------------------------------------------
% process sweeps using filtfilt
%----------------------------------------------------
R = sweep.R;
for n = 1:sweep.reps
	R{1, n} = filtfilt(dfilt.b, dfilt.a, sweep.R{1, n});
end

%----------------------------------------------------
% plot signals
%----------------------------------------------------
% plot input (stim) ...
figure(1)
subplot(211)
plot(t, sweep.S);
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Sweep Signal');
xlabel('time (ms)')
% plot response (resp)
figure(2)
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
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');


%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% build compensation for -3dB/octave gain characteristic of sweep
%---------------------------------------------------------------------
%---------------------------------------------------------------------
figure(3)

% correction for -3dB/octave (or -10dB/decade)
correct3dB = 0.5*db(Freq);
% set any ±Inf values to 0
correct3dB(isinf(correct3dB)) = 0;
% apply correction only to appropriate region of sweep
fbins = find((Freq >= sweep.start) & (Freq <= sweep.end));
sweep.mag_corr3dB = sweep.mag;
sweep.mag_corr3dB(fbins) = sweep.mag(fbins) .* power(10, correct3dB(fbins)./20);

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
