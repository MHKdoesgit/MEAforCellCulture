


st = stimdata{2}.Recording{1};


if not(isempty(st.EventStream))
    
   es = st.EventStream{1}.Events; 
    
end

ops.root = 'R:\MEA\20200605\exp1_chamb_ea1\h5';
h5filenames = dir([ops.root,filesep,'*.h5']);
[~, reindex]=sort(str2double(regexp(({h5filenames(:).name}),'\d+','match','once')));
h5filenames={h5filenames(reindex).name}'; Nfiles=numel(h5filenames);

h5file = h5filenames{4};
cfg = []; cfg.dataType = 'single';
stimdata = McsHDF5.McsData([ops.root,filesep,h5file],cfg);
h5dat = stimdata.Recording{1}.AnalogStream{1};
toVoltsFac = 10^double(h5dat.Info.Exponent(1));
mcddat = h5dat.ChannelData*toVoltsFac;
fs = round(h5dat.getSamplingRate);
stimsamples = h5dat.ChannelDataTimeStamps(end);


% 
% rawdatstreams = st.AnalogStream{1};
% 
% rawdat = rawdatstreams.ChannelData;
% rawdattime = rawdatstreams.ChannelDataTimeStamps;


% plot(h5dat.ChannelDataTimeStamps, mcddat(20,:))
signal = mcddat(20,:);


plotint = .4e-3;
madTh = 20;

%signal = rawdat(20,:);
%fs = 25e4;

%nsmooth = round(plotint*fs);
%smooths = smoothdata(signal, 'gaussian', nsmooth);
Thres = madTh * mad(diff(signal));

fonsets  = find( diff(diff(signal) >  Thres) <0 ) + 1;
foffsets = find( diff(diff(signal) < -Thres) >0 ) + 1;

pulses =  sort(unique([fonsets,foffsets]),'ascend');


% hold on;
% plot(h5dat.ChannelDataTimeStamps(pulses(1:2:end)), ypulses,'ro')
% plot(h5dat.ChannelDataTimeStamps(pulses(2:2:end)), ypulses,'kx')

newpulse = [1,find(diff(pulses) > 1e3)+1];
pulseonset = pulses(newpulse);
pulseoffset = pulses([newpulse(2:end)-1,length(pulses)]);

%pulsedur = cell(length(pulseonset),1);
for ii = 1:length(pulseoffset)
    f1 = signal(pulseoffset(ii):pulseoffset(ii)+25);
    %[~,idx] = min(abs(f1)-abs(signal(pulseonset(ii))));
    idx = interp1(unique(f1),1:length(unique(f1)),signal(pulseonset(ii)),'nearest');
    if isnan(idx)
        f1 = signal(pulseoffset(ii):pulseoffset(ii)+100);
        idx = interp1(unique(f1),1:length(unique(f1)),signal(pulseonset(ii)),'nearest');
    end
    pulseoffset(ii) = pulseoffset(ii)+idx-1;
    
    %pulsedur{ii} = pulseonset(ii):pulseoffset(ii);
end


ind = interp1(unique(f1),1:length(unique(f1)),signal(pulseonset(ii)),'nearest');


hold on;
plot(h5dat.ChannelDataTimeStamps(pulseonset), ypulses,'ro')
plot(h5dat.ChannelDataTimeStamps(pulseoffset), ypulses,'kx')

plot(h5dat.ChannelDataTimeStamps, mcddat(20,:));
hold on;

mcsignalnostim =  mcddat(20,:);
md = median(signal);
for ii = 1:length(pulseonset)
   mcsignalnostim(pulseonset(ii):pulseoffset(ii)) = md;
    plot(h5dat.ChannelDataTimeStamps(pulseonset(ii):pulseoffset(ii)),  mcddat(20,pulseonset(ii):pulseoffset(ii)),'r')
end
plot(h5dat.ChannelDataTimeStamps, mcsignalnostim,'k');


ypulsesnew = min(mcsignalnostim)+range(mcsignalnostim)/2;
plot(h5dat.ChannelDataTimeStamps, mcsignalnostim);
hold on;
plot(h5dat.ChannelDataTimeStamps(pulses(1:2:end)), ypulsesnew,'ro')
plot(h5dat.ChannelDataTimeStamps(pulses(2:2:end)), ypulsesnew,'kx')


%medfiltLoopVoltage = medfilt1(signal,8);


hold on;
plot(h5dat.ChannelDataTimeStamps(fonsets), ypulses,'ro')
plot(h5dat.ChannelDataTimeStamps(foffsets), ypulses,'kx')




%=========================================


%D = load('dataee.txt');
D = signal;
T = 1:length(D);                                        % Time Vector
Ddt = detrend(D);
Dmean = mean(Ddt);
Ddt = Ddt-Dmean;
ia =  Ddt > 10;
Ddt(ia) = Dmean;                                        % Set Spikes = Mean
Dnew = Ddt;
Tnew = T;
figure(1)
plot(Tnew, Dnew)                                        % Time Domain Plot Of Detrended & Despiked Signal
Ts = mean(diff(T));                                     % Sampling Interval
Fs = 1/Ts;                                              % Sampling Frequency
L = length(Dnew);
FTD = fft(Dnew)/L;                                      % Fourier Transform
Fv = linspace(0, 1, fix(L/2)+1)*Fs;                     % Frequency Vector
Iv = 1:length(Fv);                                      % Index Vector
figure(2)
plot(Fv, abs(FTD(Iv))*2)
axis([0  0.1    ylim])
signal = Dnew;                                          % Design & Implement Filter
Ts = mean(diff(T));
Fs = 1/Ts;                                              % Sampling Frequency (Hz)
Fn = Fs/2;                                              % Nyquist Frequency (Hz)
Wp = 0.030/Fn;                                          % Passband Frequency (Normalised)
Ws = 0.032/Fn;                                          % Stopband Frequency (Normalised)
Rp =   1;                                               % Passband Ripple (dB)
Rs =  50;                                               % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
[z,p,k] = cheby2(n,Rs,Ws, 'high');                      % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability
figure(3)
freqz(sosbp, 2^16, Fs)                                  % Filter Bode Plot
filtered_signal = filtfilt(sosbp, gbp, double(signal));         % Filter Signal
figure(4)
plot(Tnew, filtered_signal)                             % Plot Result




signal = mcddat(20,:);
mfSignal = medfilt1(signal);
badIndexes = mfSignal > 0.1;
signal(badIndexes) = mfSignal(badIndexes);







