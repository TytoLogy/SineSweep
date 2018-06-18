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
RMSsin = sqrt(2)/2;

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
fprintf('Loading data...\n');
load(fullfile(dataPath, dataName))
fprintf('...done\n');

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
%	R			{1, reps} cell array of responses to S
%	Fs			sample rate (samples/s)
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
%	S			{# frequencies, 1} cell array of tone signals (empty for now)
%	R			{# frequencies, reps} cell array of responses to S
%	Fs			sample rate (samples/s) - determined after NIDAQ init
%----------------------------------------------------
% mic struct:
% 				these are parameters for calibration microphone
%----------------------------------------------------
%	gain		mic gain (dB)
%	sense		mic sensitivity (Volts/Pascal)
%	VtoPa		conversion factor (accounts for gain and mic sense)
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
%% process tone data
%---------------------------------------------------------------------
%---------------------------------------------------------------------

% storage for amplitude (mV, Pascal, dB SPL), phases (us)
tone.amp = zeros(length(tone.freq), tone.reps);
tone.ampPa = zeros(length(tone.freq), tone.reps);
tone.ampdB = zeros(length(tone.freq), tone.reps);
tone.phi = zeros(length(tone.freq), tone.reps);

% loop through frequencies
for f = 1:length(tone.freq)
	for r = 1:tone.reps
		% filter data after applying a sin2 window to it (prevents
		% onset/offset transients)
		fData = filtfilt(dfilt.b, dfilt.a, sin2array(tone.R{f, r}, 1, tone.Fs));
		% determine the magnitude and phase of the response
		[tone.amp(f, r), tone.phi(f, r)] = ...
												fitsinvec(fData, 1, tone.Fs, tone.freq(f));
		% Pascal
		tone.ampPa(f, r) = mic.VtoPa*tone.amp(f, r);
		% dB SPL
		tone.ampdB(f, r) = 20*log10( (RMSsin*tone.ampPa(f, r)) / 20e-6);
		% phase in us
		tone.phi(f, r) = 1.0e6 * unwrap(tone.phi(f, r)) ./ (2*pi*tone.freq(f));
	end
end

if tone.reps > 1
	% compute mean, std. dev
	tone.amp_mean = mean(tone.amp, 2);
	tone.amp_std = std(tone.amp, 0, 2);	
	tone.ampPa_mean = mean(tone.ampPa, 2);
	tone.ampPa_std = std(tone.ampPa, 0, 2);
	tone.ampdB_mean = mean(tone.ampdB, 2);
	tone.ampdB_std = std(tone.ampdB, 0, 2);
	tone.phi_mean = mean(tone.phi, 2);
	tone.phi_std = std(tone.phi, 0, 2);
else
	% just use single values
	tone.amp_mean = tone.amp;
	tone.amp_std = 0;	
	tone.ampPa_mean = tone.ampPa;
	tone.ampPa_std = 0;
	tone.ampdB_mean = tone.ampdB;
	tone.ampdB_std = 0;
	tone.phi_mean = tone.phi;
	tone.phi_std = 0;
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% process sweep data
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%----------------------------------------------------
% process sweeps using filtfilt
%----------------------------------------------------
% filtering sweep data
R = sweep.R;
for n = 1:sweep.reps
	R{1, n} = filtfilt(dfilt.b, dfilt.a, sweep.R{1, n});
end

%----------------------------------------------------
% spectra
%----------------------------------------------------
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
stim.Nfft = length(sweep.S);
stim.Sfft = fft(sweep.S, stim.Nfft);
resp.Nfft = length(R{1, 1});
resp.Sfft = fft(R{1, 1}, resp.Nfft);

% compute magnitude spectrum
% non-redundant points are kept
Nunique = (stim.Nfft/2) + 1;
stim.mag = stim.Sfft(1:Nunique);
% and magnitude computed
stim.mag = 2*abs(stim.mag)/stim.Nfft;
stim.phi = angle(stim.Sfft(1:Nunique));

% non-redundant points are kept
Nunique = (resp.Nfft/2) + 1;
resp.mag = resp.Sfft(1:Nunique);
% and magnitude computed
resp.mag = 2*abs(resp.mag)/resp.Nfft;
resp.phi = angle(resp.Sfft(1:Nunique));
Fnyq = sweep.Fs/2;

% compensation for -3dB/octave (or -10dB/decade) gain characteristic of sweep
% need vector of Frequencies
Freq = Fnyq*linspace(0, 1, length(stim.mag));
% correction for -3dB/octave
correct3dB = 0.5*db(Freq);
% set any �Inf values to 0
correct3dB(isinf(correct3dB)) = 0;
% apply correction only to appropriate region of sweep (apply correction in
% linear units (not dB)
fbins = find((Freq >= sweep.start) & (Freq <= sweep.end));
stim.mag_corr = stim.mag;
stim.mag_corr(fbins) = stim.mag(fbins) .* power(10, correct3dB(fbins)./20);
resp.mag_corr = resp.mag;
resp.mag_corr(fbins) = resp.mag(fbins) .* power(10, correct3dB(fbins)./20);


%% Scale to dB SPL using info from tone tests
% first, need to find closest FFT frequency bin to test tone frequencies
minbin = zeros(size(tone.freq));
minfreq = zeros(size(tone.freq));
minval = zeros(size(tone.freq));
for f = 1:length(tone.freq)
	% compute difference between FFT freq and test tone freq
	deltaF = Freq - tone.freq(f);
	% find min absolute difference and return bin in Freq
	[~, minbin(f)] = min(abs(deltaF));
	% get frequency for this bin
	minfreq(f) = Freq(minbin(f));
	% and FFT value at this bin
	minval(f) = resp.mag_corr(minbin(f));
end
% compute ratio of tone.magPa_mean and minval
pafactors = tone.ampPa_mean ./ minval';

% check if anything is odd (values are above tolerance)
if min(diff(pafactors)) > 2
	error('pafactors exceeds limit of 2')
else
	resp.Pa_scale = mean(pafactors);
end
% apply correction to response spectrum
resp.magPa_corr = resp.Pa_scale * resp.mag_corr;
resp.magdBSPL_corr = dbspl(RMSsin * resp.magPa_corr);

%----------------------------------------------------
%% plot signals
%----------------------------------------------------
% generate time vector for sweep plotting
t = 0:(1/sweep.Fs):( (0.001*sweep.dur) - (1/sweep.Fs));
% plot input (stim) ...
figure(1)
subplot(211)
plot(t, sweep.S);
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Sweep Stimulus');
xlabel('time (ms)')
ylabel('V')
% plot response (resp)
figure(2)
subplot(211)
plot(t, R{1, 1});
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Sweep Response')
xlabel('time (ms)')
ylabel('V');
%----------------------------------------------------
% plot uncorrected spectra
%----------------------------------------------------
figure(1)
subplot(212)
semilogx(Freq, db(stim.mag));
xlim([0 Fnyq]);
ylim([-100 0]);
title('Stimulus Spectrum')
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');
figure(2)
subplot(212)
semilogx(Freq, db(resp.mag));
xlim([0 Fnyq]);
ylim([-100 0]);
title('Response Spectrum')
xlabel('Frequency (Hz)')
ylabel('dB')
grid('on');
%----------------------------------------------------
% plot -3dB corrected spectra
%----------------------------------------------------
figure(1)
subplot(212)
hold on
	semilogx(Freq, db(stim.mag_corr));
hold off
legend({'uncorrected', 'corrected'})
figure(2)
subplot(212)
hold on
	semilogx(Freq, db(resp.mag_corr));
hold off
legend({'uncorrected', 'corrected'})

%----------------------------------------------------
%% plot response in dB SPL
%----------------------------------------------------
figure(3)
% plot response (resp)
subplot(211)
plot(t, R{1, 1});
xlim([min(t) max(t)])
ylim([-1.1 1.1])
title('Sweep Response')
xlabel('time (ms)')
ylabel('Pa');
% plot spectrum (dB SPL, corrected)
subplot(212)
semilogx(Freq, resp.magdBSPL_corr);
xlim([0 Fnyq]);
title('Response dB SPL')
xlabel('Frequency (Hz)')
ylabel('dB SPL')
grid('on');
