%% fMRI data plotting
clear all; close all; clc;

% Import data
load('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\Data\Bayesian_topology_results\fmri_newBf_results_20.mat');

Ts = 1.986622; % Sample time of time series

% Plot time series
figure(1)
for i = 1:20
    subplot(5,4,i)
    plot(Xb(:,i))
    grid on
    xlabel('Time [s]')
    ylabel(['Signal ',num2str(i)])
end
suptitle('Timeseries signals')

% Take one signal to do fft
if false    
    Signal = Xb(:,19);
    fftSignal = fft(Signal);
    figure()
    plot(1:length(fftSignal),db(fftSignal));
    grid on
    xlim([0 length(fftSignal)/2])
    ylim([0 50])
end

% Create iddata() objects from all timeseries 
for i = 1:20
   signal{i} = iddata(Xb(:,i),[],Ts); 
end



