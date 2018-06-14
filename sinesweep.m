%---------------------------------------------------------------------
%{

Testing "notebook" for frequency response measurement with swept sinusoid

Process:

(1) User confirms settings


%}
%---------------------------------------------------------------------

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Settings
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% sample rate and Nyquist frequency
Fs = 100000;
Fnyq = Fs/2;

% sweep parameters
%	dur		duration (ms)
%	start		start frequency (Hz)
%	end		end frequency (Hz)
%	mode		sweep mode: 'linear', 'log'
%	ramp		ramp onset/offset duration (ms)
%	S			signal vector (empty for now)
%
% Some things to note for use of this in frequency response measurements:
%	- sweep frequencies should start slightly below and end slightly above
%	the desired characterization range (if device-under-test supports or
%	will not be damaged at these frequency limits!!!)
sweep.dur = 1000;
sweep.start = 4000;
sweep.end = 45000;
sweep.mode = 'log';
sweep.ramp = 0.1;
sweep.Fs = Fs;
sweep.S = [];

% these are dummy parameters for microphone
%	gain		mic gain (dB)
%	sense		mic sensitivity (Volts/Pascal)
%	VtoPa		conversion factor
mic.gain = 0;
mic.sense = 1;
mic.VtoPa = (1/invdb(mic.gain)) * (1 / mic.sense);

% time interval
dt = 1/Fs;

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% build some things
%---------------------------------------------------------------------
%---------------------------------------------------------------------

% test filter parameters
[b, a] = butter(1, [1000 1500]./Fnyq, 'stop');


%----------------------------------------------------
% synthesize sweep
%----------------------------------------------------
% generate time vector
t = 0:dt:( (0.001*sweep.dur) - dt);
% generate raw signal
sweep.S = chirp(t, sweep.start, sweep.dur/1000, sweep.end, 'logarithmic');
% ramp on/off if requested
if sweep.ramp > 0
	sweep.S = sin2array(sweep.S, sweep.ramp, Fs);
end

%----------------------------------------------------
% process sweep using filter to simulate response
%----------------------------------------------------
resp.S = filter(b, a, sweep.S);

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
plot(t, resp.S);
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
if ~even(length(resp.S))
	% if not, pad with a 0
	resp.S = [resp.S 0];
end
sweep.Nfft = length(sweep.S);
sweep.Sfft = fft(sweep.S, sweep.Nfft);
resp.Nfft = length(resp.S);
resp.Sfft = fft(resp.S, resp.Nfft);

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



