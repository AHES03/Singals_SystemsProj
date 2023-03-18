% CW2 Tutorial
load 1-AM.mat
% for example, in this case, I choose 
%1. channel 5 ECG
myECG = b(:,5);
%2. channel 2 PCG
myPCG = b(:,2);
%3. channel 4 or 9 respiratory
respSig1 = b(:,4);
respSig2 = b(:,9);
myTime = b(:,1);
Fs = 20000; 
% Visualise
figure(1)
subplot(221), plot(myTime,myECG), title('ECG')
subplot(222), plot(myTime,myPCG), title('PCG')
subplot(223), plot(myTime,respSig1), title('raw respiration')
subplot(224), plot(myTime,respSig2), title('Filtered respiration')

% high pass filter
[b,a] = butter(3,0.1/(Fs/2),'high');
myECG = filtfilt(b,a,myECG);
myPCG = filtfilt(b,a,myPCG);

figure(2)
subplot(211), plot(myTime,myECG), title('filtered ECG')
subplot(212), plot(myTime,myPCG), title('filtered PCG')

% Find peaks (one of the many methods)
[pks,locs] = findpeaks(myECG);
subplot(211), plot(myTime,myECG), title('filtered ECG')
hold on
subplot(211), plot(myTime(locs),pks,'r*')
close all
% add parameter
[pks1,locs1] = findpeaks(myECG, Fs, 'MinPeakDistance',0.5);
figure(3)
plot(myTime,myECG), hold on
plot(locs1,pks1,'r*')

close all
plot(myTime,respSig2)
hold on
plot(locs1,pks1,'r--')

% Using correlation
% N = round(length(respSig2)/length(locs1))
% ECGpeaks = interp(pks1,N);
% corrcoef(respSig2,ECGpeaks(1:length(respSig2)))
%using peak count
close all
newfs = length(pks1)/max(myTime)
[pks2,locs2] = findpeaks(pks1,newfs, 'MinPeakDistance',3);
plot(locs1,pks1)
hold on
plot(locs2,pks2,'*k')
ecgRate = 60*length(pks2)/max(myTime)

close all
% peaks of respiration 
[pks3,locs3] = findpeaks(respSig2,Fs, 'MinPeakDistance',3);
plot(myTime,respSig2)
hold on
plot(locs3,pks3,'*k')
trueRate = 60*length(pks3)/max(myTime)

Error = abs(ecgRate-trueRate)
close all
% What about the frequency
% What is the spectrum of respiration
L = length(respSig2);
f = Fs*(0:(L/2))/L;
respFFT = fft(respSig2);
close all
if mod(L,2) == 1
    L= L+1;
end 
Px = abs(respFFT(1:L/2));
close all
plot(log10(f(1:L/2)),Px)

% mean frequency
mn = meanfreq(respSig2,Fs);
freqrate = mn*60 % does this look correct to you?
md = medfreq(respSig2,Fs);
freqrate2 = md*60 % does this look correct to you?