%%
% 
%  PREFORMATTED
%  TEXT
% 
% CW2 Tutorial
load ('C:\Users\dabst\OneDrive\Desktop\signals\Data\Data\1\1-AP.mat')
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
 plot(myTime,myECG), title('ECG')
 plot(myTime,myPCG), title('PCG')
 plot(myTime,respSig2), title('Filtered respiration')

%finding cut off frequency of ECG
L = length(myECG);
f = Fs*(0:(L/2))/L;
ecgFFT = fft(myECG);
if mod(L,2) == 1
    L= L+1;
end 
ecgPower = abs(ecgFFT(1:L/2));
totalpower_ecg=sum(ecgPower);
plot(log10(f(1:L/2)),ecgPower),title("ECG in f")
energy=0
for i= 1:length(f)
    energy=energy+ecgPower(i);
    if energy>=0.5*totalpower_ecg
        cutoff=f(i)
        break
    end
end

%finding bandwidth
energy=0
for i= 1:length(f)
    energy=energy+ecgPower(i);
    if energy>=0.95*totalpower_ecg
        bandwidth=f(i)
        break
    end
end


% high pass filter
[b,a] = butter(3,cutoff/(Fs/2),'high');
myECG = filtfilt(b,a,myECG);
myPCG = filtfilt(b,a,myPCG);

plot(myTime,myECG), title('filtered ECG')
%plot(myTime,myPCG), title('filtered PCG')


%[pks,locs] = findpeaks(myECG);
% plot(myTime,myECG), title('filtered ECG')
%hold on
%plot(myTime(locs),pks,'r*')
%hold off
% add parameter
% Find peaks (one of the many methods)
[pks1,locs1] = findpeaks(myECG, Fs, 'MinPeakHeight',0.15);
plot(myTime,myECG) 
hold on
plot(locs1,pks1,'r*'),title('filtiered ECG')

hold off
plot(myTime,respSig2)
hold on
plot(locs1,pks1,'r--')


% Using correlation
% N = round(length(respSig2)/length(locs1))
% ECGpeaks = interp(pks1,N);
% corrcoef(respSig2,ECGpeaks(1:length(respSig2)))
%using peak count
close all
newfs = length(pks1)/max(myTime) %new sampling frequency is length of peaks from ecg over time
[pks2,locs2] = findpeaks(pks1,newfs, 'MinPeakDistance',3);
plot(locs1,pks1)
hold on
plot(locs2,pks2,'*k')
ecgRate = 60*length(pks2)/max(myTime)

hold off


% peaks of respiration 
[pks3,locs3] = findpeaks(respSig2,Fs, 'MinPeakHeight',0.14);
plot(myTime,respSig2)
hold on
plot(locs3,pks3,'*k')
trueRate = 60*length(pks3)/max(myTime)

Error = abs(ecgRate-trueRate)
hold off
% What about the frequency
% What is the spectrum of respiration
L = length(respSig2);
f = Fs*(0:(L/2))/L;
respFFT = fft(respSig2);
if mod(L,2) == 1
    L= L+1;
end 
Px = abs(respFFT(1:L/2));
plot(log10(f(1:L/2)),Px)

% mean frequency
mn = meanfreq(respSig2,Fs);
freqrate = mn*60 % does this look correct to you?
md = medfreq(respSig2,Fs);
freqrate2 = md*60 % does this look correct to you?
