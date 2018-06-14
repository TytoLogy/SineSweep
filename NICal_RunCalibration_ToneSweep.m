%--------------------------------------------------------------------------
% NICal_RunCalibration_ToneSweep.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% Runs the speaker calibration
% if FR Correction is selected, apply mic correction using data from
% MicrophoneCal program (earphone fr data)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	16 Oct 2014 from NICal_RunCalibration_ToneStack,	SJS
%
% Revisions:
%	1 Feb 2017 (SJS): updated for session interface
%	9 Feb 2017 (SJS): modifying stimulus to include pre/post stim time
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Global Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_Constants;
% local settings
% set the COMPLETE flag to 0
COMPLETE = 0;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialization Scripts
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%---------------------------------------------
% Load the settings and constants 
%---------------------------------------------
NICal_settings;
% save the GUI handle information
guidata(hObject, handles);

%-----------------------------------------------------------------------
% check output  file - if it exists, check with user
%-----------------------------------------------------------------------
calfile = check_output_file(handles);
if isequal(calfile, 0)
	return
else
	handles.cal.calfile = calfile;
	update_ui_str(handles.CalFileCtrl, handles.cal.calfile);
	guidata(hObject, handles);
end


%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Start DAQ things
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[handles, initFlag] = NICal_NIinit(handles);
guidata(hObject, handles);
if initFlag == 0
	warning('NICAL:HW', '%s: NIinit failure', mfilename)
	return
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = handles.iodev.Fs / 2;
% passband definition
handles.cal.fband = [handles.cal.InputHPFc handles.cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.cal.fcoeffb, handles.cal.fcoeffa] = ...
					butter(handles.cal.forder, handles.cal.fband, 'bandpass');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% if raw data are to be saved, initialize the file
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if handles.cal.SaveRawData
	[pathstr, fname, fext] = fileparts(handles.cal.calfile);
	rawfile = fullfile(pathstr, [fname '.dat']);
	fp = fopen(rawfile, 'w');
	writeStruct(fp, handles.cal, 'cal');
	fclose(fp);
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup caldata struct for storing the calibration data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_caldata_init;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Preallocate some arrays that are used locally
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set the start and end bins for the calibration
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[start_bin, end_bin] = startendbins(handles);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% create null stimulus and time vector for plots, set up plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% calculate difference between SweepDuration and StimDelay + StimDuration
% this will be used to pad end of stimulus
% this is due to the session DAQ interface using # of cued output samples to
% determine the number of samples to read in. For the legacy interface,
% this shouldn't make a difference
PostDuration = handles.cal.SweepDuration - ...
						(handles.cal.StimDelay + handles.cal.StimDuration);
% make sure PostDuration is ok (greater than or equal to 0)
if PostDuration < 0
	errordlg('SweepDuration must be greater than StimDelay + StimDuration');
	NICal_NIexit;
	COMPLETE = 0;
	return
else
	% if ok, create poststim
	poststim = syn_null(PostDuration, handles.iodev.Fs, 0);
end
% create null stimulus
zerostim = syn_null(handles.cal.StimDuration, handles.iodev.Fs, 0);
% insert stim delay
zerostim = insert_delay(zerostim, handles.cal.StimDelay, handles.iodev.Fs);
% append post-stim
zerostim = [zerostim poststim];
% downsample (no need to plot all points)
zerostim = downsample(zerostim, handles.cal.deciFactor);
% downsample-factor adjusted sample interval
dt = handles.cal.deciFactor/handles.iodev.Fs;
% time vector for stimulus plots
tvec_stim = 1000*dt*(0:(length(zerostim)-1));
% fake acquired data
zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
zeroacq = downsample(zeroacq, handles.cal.deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
tvec_acq = 1000*dt*(0:(acqpts-1));
% compute # of points per sweep
SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
% stimulus start and end points
stim_start = ms2bin(handles.cal.StimDelay, handles.iodev.Fs);
stim_end = stim_start + outpts - 1;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Build null output array
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% synthesize the L sine wave;
Nullstim = syn_null(handles.cal.StimDuration, handles.iodev.Fs, 1);
% scale the sound
Nullstim = 0 * Nullstim;
% insert delay
Nullstim = insert_delay(Nullstim, handles.cal.StimDelay, handles.iodev.Fs);
Nullstim_downsample =  downsample(Nullstim(1, :), handles.cal.deciFactor);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup Plots
%-----------------------------------------------------------------------
% to speed up plotting, the vectors Lacq, Racq, tvec_acq, L/Rfft, fvec
% are pre-allocated and then those arrys are used as XDataSource and
% YDataSource for the respective plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-------------------------------------------------------
% create arrays for plotting and plot them
%-------------------------------------------------------
% stim
Lstim = zerostim;
Rstim = zerostim;
% acq
Lacq = zeroacq;
Racq = zeroacq;
% FFT
nfft = length(start_bin:end_bin);
tmp = zeros(1, nfft);
[fvec, Lfft] = daqdbfft(tmp, handles.iodev.Fs, nfft); %#ok<ASGLU>
[fvec, Rfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
% convert fvec to kHz
fvec = 0.001 * fvec;
clear tmp

%-------------------------------------------------------
% plot null data, save handles for time-domain plots
%-------------------------------------------------------
% stimulus
H.Lstim = plot(handles.Lstimplot, tvec_stim, Lstim, 'g');
set(H.Lstim, 'XDataSource', 'tvec_stim', 'YDataSource', 'Lstim');
ylabel(handles.Lstimplot, 'V');
H.Rstim = plot(handles.Rstimplot, tvec_stim, Rstim, 'r');
set(H.Rstim, 'XDataSource', 'tvec_stim', 'YDataSource', 'Rstim');
% response
H.Lacq = plot(handles.Lmicplot, tvec_acq, Lacq, 'g');
set(H.Lacq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Lacq');
ylabel(handles.Lmicplot, 'V')
H.Racq = plot(handles.Rmicplot, tvec_acq, Racq, 'r');
set(H.Racq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Racq');

%-------------------------------------------------------
% plot null data, save handles for frequency-domain plots
%-------------------------------------------------------
H.Lfft = plot(handles.Lfftplot, fvec, Lfft);
set(H.Lfft, 'XDataSource', 'fvec', 'YDataSource', 'Lfft');
xlabel(handles.Lfftplot, 'Frequency (kHz)')
ylabel(handles.Lfftplot, 'dBV')
H.Rfft = plot(handles.Rfftplot, fvec, Rfft);
set(H.Rfft, 'XDataSource', 'fvec', 'YDataSource', 'Rfft');
xlabel(handles.Rfftplot, 'Frequency (kHz)');

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% setup attenuation
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% check if AttenFix is set and that it is in the range of [0, MAX_ATTEN]
if handles.cal.AttenFix && between(handles.cal.AttenFixValue, 0, MAX_ATTEN)
	Latten = handles.cal.AttenFixValue;
	Ratten = handles.cal.AttenFixValue;
else
	% use StartAtten
	if ~between(handles.cal.AttenFixValue, 0, MAX_ATTEN)
		warning('NICal:Atten', [mfilename ...
					': AttenFixValue out of range, using default StartAtten value'])
	end
	% set the adjustable starting attenuator values	
	Latten = handles.cal.StartAtten;
	Ratten = handles.cal.StartAtten;
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Now initiate sweeps
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
stopFlag = 0;
%---------------------------------------------------
% make a local copy of the cal settings structure
%---------------------------------------------------
cal = handles.cal;
%---------------------------------------------------
% make local copy of iodev  control struct
%---------------------------------------------------
iodev = handles.iodev;
% update the frequency display value
update_ui_str(handles.FreqValText, 'Tone Sweep');

% check for abort button press
if read_ui_val(handles.AbortCtrl) == 1
	% if so, stop
	disp('abortion detected')
	return
end

%------------------------------------------------------------------
% if cal.Side is 1 or 3 (LEFT or BOTH), calibrate L channel
% Side is set by the SideCtrl pulldown under CalibrationSettings
% 		cal.Side == 1 is LEFT
% 		cal.Side == 2 is RIGHT
% 		cal.Side == 3 is BOTH channels, 
%------------------------------------------------------------------
if cal.Side == 1 || cal.Side == 3
	%LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
	% set input channel: if InputChannel is set to 3 (both), 
	% use left input when calibrating Left side
	%LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
	if handles.cal.InputChannel == 3
		inChan = 1;
	else
		inChan = handles.cal.InputChannel;
	end
	
	%-------------------------------------------------------
	% Build stimulus output array for this frequency
	%-------------------------------------------------------
	% synthesize the sweep
	% need time vector
	t = (1/iodev.Fs) * (0:(ms2samples(cal.StimDuration, iodev.Fs)-1));
	csig = chirp(t, Freqs(1), 0.001*cal.StimDuration, Freqs(end));
	S = [csig; 0*csig];
	% scale the sound
	S = cal.DAscale * S;
	% apply the sin^2 amplitude envelope to the stimulus before adding 
	% pre and post zeros
	S = sin2array(S, cal.StimRamp, iodev.Fs);
	% insert delay, add zeros to pad end
	S = [insert_delay(S, cal.StimDelay, iodev.Fs) ...
								syn_null(PostDuration, iodev.Fs, 1)];
	% save in Satt
	Satt = S;
	% plot the stimuli - set R stim to zero
	Lstim = downsample(S(1, :), handles.cal.deciFactor);
	Rstim = zerostim;
	refreshdata(H.Lstim, 'caller');
	refreshdata(H.Rstim, 'caller');
	%-------------------------------------------------------
	% set the L attenuator value.
	%-------------------------------------------------------
	% no need to test attenuation but, 
	% do need to set the attenuators
	Satt(1, :) = handles.attfunction(S(1, :), Latten);
	Satt(2, :) = handles.attfunction(S(2, :), MAX_ATTEN);
	update_ui_str(handles.LAttenText, Latten);
	update_ui_str(handles.RAttenText, MAX_ATTEN);
	pause(0.001*cal.ISI);
	%-------------------------------------------------------
	% now, collect the data for frequency FREQ, LEFT channel
	%-------------------------------------------------------
	for rep = 1:cal.Nreps
		% update the reps display value
		update_ui_str(handles.RepNumText, sprintf('%d L', rep));
		% play the sound;
		if handles.DAQSESSION
			[resp, indx] = handles.iofunction(iodev, Satt, ...
													handles.cal.SweepDuration);
		else
			[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
		end
		% filter the data if asked
		if handles.cal.InputFilter
			tmp = sin2array(resp{L}, 1, iodev.Fs);
			resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
			if handles.cal.MeasureLeak
				tmp = sin2array(resp{R}, 1, iodev.Fs);
				resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
			end
			clear tmp
		end
		% plot the response and FFT
		Lacq = downsample(resp{L}, handles.cal.deciFactor);
		refreshdata(H.Lacq, 'caller');
		[tmpf, Lfft] = daqdbfft(resp{L}(start_bin:end_bin), ...
										iodev.Fs, length(resp{L}(start_bin:end_bin)));
		refreshdata(H.Lfft, 'caller');
		if handles.cal.MeasureLeak
			Racq = downsample(resp{R}, handles.cal.deciFactor);
			refreshdata(H.Racq, 'caller');
			[tmpf, Rfft] = daqdbfft(resp{R}(start_bin:end_bin), ...
											iodev.Fs, length(resp{R}(start_bin:end_bin)));
			refreshdata(H.Rfft, 'caller');
		end
		drawnow
		% draw spectrogram
		axes(handles.Lspecgram); %#ok<*LAXES>
		myspectrogram(resp{L}, iodev.Fs, ...
								[10 5], @hamming, handles.SpectrumWindow, ...
								[-100 -1], false, 'default', false, 'per');
		if handles.cal.MeasureLeak
			axes(handles.Rspecgram);
			myspectrogram(resp{R}, iodev.Fs, ...
									[10 5], @hamming, handles.SpectrumWindow, ...
									[-100 -1], false, 'default', false, 'per');			
		end
		% save raw data
		if handles.cal.SaveRawData
			fp = fopen(rawfile, 'a');
			if handles.cal.MeasureLeak
				writeCell(fp, resp);
			else
				writeCell(fp, resp(L));
			end
			fclose(fp);
		end
		% Pause for ISI
		pause(0.001*cal.ISI);

		% check for abort button press
		if read_ui_val(handles.AbortCtrl) == 1
			% if so, stop
			disp('abortion detected')
			break
		end
	end
	%---------------------------------------------------------------------
	% now, collect the background data for frequency FREQ, LEFT channel
	%---------------------------------------------------------------------
	if handles.cal.CollectBackground
		Lstim = Nullstim_downsample;
		Rstim = Nullstim_downsample;
		refreshdata(H.Lstim, 'caller');
		refreshdata(H.Rstim, 'caller');
		for rep = 1:cal.Nreps
			% update the reps display value
			update_ui_str(handles.RepNumText, sprintf('%d L (bg)', rep));
			% play the sound;
			if handles.DAQSESSION
				[resp, indx] = handles.iofunction(iodev, Nullstim, ...
														handles.cal.SweepDuration);
			else
				[resp, indx] = handles.iofunction(iodev, Nullstim, SweepPoints);
			end
			% filter the data if asked
			if handles.cal.InputFilter
				tmp = sin2array(resp{L}, 1, iodev.Fs);
				resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				tmp = sin2array(resp{R}, 1, iodev.Fs);
				resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				clear tmp
			end
			% plot the response
			Lacq = downsample(resp{L}, handles.cal.deciFactor);
			Racq = downsample(resp{R}, handles.cal.deciFactor);
			refreshdata(H.Lacq, 'caller');
			refreshdata(H.Racq, 'caller');
			% plot fft
			figure(10)
			[tmpf, tmpm] = daqdbfft(resp{L}(start_bin:end_bin), ...
									iodev.Fs, length(resp{L}(start_bin:end_bin)));
			plot(tmpf, tmpm);
			title('Left Background')
			% Pause for ISI
			if handles.cal.SaveRawData
				fp = fopen(rawfile, 'a');
				writeCell(fp, resp); 				
				fclose(fp);
			end
			pause(0.001*cal.ISI);
			% check for abort button press
			if read_ui_val(handles.AbortCtrl) == 1
				% if so, stop
				disp('abortion detected')
				break
			end
		end		
	end
end		% END OF L CHANNEL	

%------------------------------------
% check for abort button press
%------------------------------------
if read_ui_val(handles.AbortCtrl) == 1
	% if so, stop
	disp('abortion detected')
	return
end
%------------------------------------
% pause for ISI
%------------------------------------
pause(0.001*cal.ISI);

%------------------------------------------------------------------
% if cal.Side is 2 or 3 (RIGHT or BOTH), calibrate *R* channel
%------------------------------------------------------------------	
if cal.Side == 2 || cal.Side == 3
	%RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
	% set input channel: if INputChannel is set to 3 (both), use right
	if handles.cal.InputChannel == 3
		inChan = 2;
	else
		inChan = handles.cal.InputChannel;
	end
	%-------------------------------------------------------------
	% synthesize the R sine wave;
	%-------------------------------------------------------------
	% synthesize 
	% need time vector
	t = 0:(1/iodev.Fs):0.001*cal.StimDuration;
	csig = chirp(t, Freqs(1), 0.001*cal.StimDuration, Freqs(end));
	S = [0*csig; csig];
	% scale the sound
	S = cal.DAscale * S;
	% apply the sin^2 amplitude envelope to the stimulus
	S = sin2array(S, cal.StimRamp, iodev.Fs);
	% insert delay, add zeros to pad end
	S = [insert_delay(S, cal.StimDelay, iodev.Fs) ...
								syn_null(PostDuration, iodev.Fs, 1)];
	% save in Satt
	Satt = S;
	% plot the stimulus arrays
	Lstim = zerostim;
	Rstim = downsample(S(2, :), handles.cal.deciFactor);
	refreshdata(H.Lstim, 'caller');
	refreshdata(H.Rstim, 'caller');
	% set R attenuator value.
	% no need to test attenuation but, 
	% do need to set the attenuators
	Satt(1, :) = handles.attfunction(S(1, :), MAX_ATTEN);
	Satt(2, :) = handles.attfunction(S(2, :), Ratten);
	update_ui_str(handles.LAttenText, MAX_ATTEN);
	update_ui_str(handles.RAttenText, Ratten);
	pause(0.001*cal.ISI);
	% now, collect the data for frequency FREQ, RIGHT headphone
	for rep = 1:cal.Nreps
		% update the reps display value
		update_ui_str(handles.RepNumText, sprintf('%d R', rep));
		% play the sound;
		if handles.DAQSESSION
			[resp, indx] = handles.iofunction(iodev, Satt, ...
													handles.cal.SweepDuration);
		else
			[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
		end
		% filter the data if asked
		if handles.cal.InputFilter
			if handles.cal.MeasureLeak
				tmp = sin2array(resp{L}, 1, iodev.Fs);
				resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
			end
			tmp = sin2array(resp{R}, 1, iodev.Fs);
			resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
			clear tmp
		end
		% plot the response
		if handles.cal.MeasureLeak
			Lacq = downsample(resp{L}, handles.cal.deciFactor);
			refreshdata(H.Lacq, 'caller');
			[~, Lfft] = daqdbfft(resp{L}(start_bin:end_bin), iodev.Fs, ...
												length(resp{L}(start_bin:end_bin)));
			refreshdata(H.Lfft, 'caller');
		end
		Racq = downsample(resp{R}, handles.cal.deciFactor);
		refreshdata(H.Racq, 'caller');

		[tmpf, Rfft] = daqdbfft(resp{R}(start_bin:end_bin), iodev.Fs, ...
											length(resp{R}(start_bin:end_bin)));
		refreshdata(H.Rfft, 'caller');
		drawnow
		% draw spectrogram
		axes(handles.Rspecgram);
		myspectrogram(resp{R}, iodev.Fs, ...
								[10 5], @hamming, handles.SpectrumWindow, ...
								[-100 -1], false, 'default', false, 'per');
		% save raw data
		if handles.cal.SaveRawData
			fp = fopen(rawfile, 'a');
			if handles.cal.MeasureLeak
				writeCell(fp, resp);
			else
				writeCell(fp, resp(R));
			end
			fclose(fp);
		end
		% pause for ISI (convert to seconds)
		pause(0.001*cal.ISI);
		% check for abort button press
		if read_ui_val(handles.AbortCtrl) == 1
			% if so, stop
			disp('abortion detected')
			break
		end
	end

	%---------------------------------------------------------------------
	% now, collect the background data for frequency FREQ, RIGHT channel
	%---------------------------------------------------------------------
	if handles.cal.CollectBackground
		Lstim = Nullstim_downsample;
		Rstim = Nullstim_downsample;
		refreshdata(H.Lstim, 'caller');
		refreshdata(H.Rstim, 'caller');
		for rep = 1:cal.Nreps
			% update the reps display value
			update_ui_str(handles.RepNumText, sprintf('%d R (bg)', rep));
			% play the sound;
			if handles.DAQSESSION
				[resp, indx] = handles.iofunction(iodev, Nullstim, ...
														handles.cal.SweepDuration);
			else
				[resp, indx] = handles.iofunction(iodev, Nullstim, SweepPoints);
			end
			% filter the data if asked
			if handles.cal.InputFilter
				tmp = sin2array(resp{L}, 1, iodev.Fs);
				resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				tmp = sin2array(resp{R}, 1, iodev.Fs);
				resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				clear tmp
			end
			% plot the response
			Lacq = downsample(resp{L}, handles.cal.deciFactor);
			Racq = downsample(resp{R}, handles.cal.deciFactor);
			refreshdata(H.Lacq, 'caller');
			refreshdata(H.Racq, 'caller');
			% plot fft
			figure(10)
			[tmpf, tmpm] = daqdbfft(resp{R}(start_bin:end_bin), ...
										iodev.Fs, length(resp{R}(start_bin:end_bin)));
			plot(tmpf, tmpm);
			title('Right Background')
			% save raw data
			if handles.cal.SaveRawData
				fp = fopen(rawfile, 'a');
				writeCell(fp, resp); 				
				fclose(fp);
			end
			% Pause for ISI
			pause(0.001*cal.ISI);
			% check for abort button press
			if read_ui_val(handles.AbortCtrl) == 1
				% if so, stop
				disp('abortion detected')
				break
			end
		end
	end
% END OF R CAL
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close TDT objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_NIexit;
disp('Finished.')


