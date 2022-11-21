%% Rank-reduced noise identification on fMRI data using the SLS algorithm
% Created by Job Meijer - j.b.t.meijer@student.tue.nl
clear all; close all; clc;

%% Load fMRI data (insert folder path in dataFolderLocation)
importSettings.dataFolderLocation = 'C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\5. Mozart_ICA_timeseries_data\Mozart_ICA_timeseries_data\Timeseries';
importSettings.dataFilesLocation = [importSettings.dataFolderLocation, '\*.txt'];       % Add \*.txt to select all txt files in dir
importSettings.matFiles = dir(fullfile(importSettings.dataFilesLocation));
importSettings.initDir = cd;    % Save initial directory
cd(importSettings.dataFolderLocation);

% Import all data and sort on person and scan
importSettings.nPersons = 16;
importSettings.nScansPerPerson = 4;
for nPerson = 1:importSettings.nPersons
    for nScan = 1:importSettings.nScansPerPerson
        fMriPerson{nPerson}.Scan{nScan} = importdata(importSettings.matFiles((nPerson-1)*4+(nScan)).name);
    end
end
clear nPerson nScan;

cd(importSettings.initDir);     % return to initial directory

%% Plot time series of all 16 persons (4 scans with timeseries of 20 signals with length of 300 samples)
if false
    desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    myGroup = desktop.addGroup('fMRI ICA signals');
    desktop.setGroupDocked('fMRI ICA signals', 0);
    myDim   = java.awt.Dimension(4, 2);   % 4 columns, 2 rows
    % 1: Maximized, 2: Tiled, 3: Floating
    desktop.setDocumentArrangement('fMRI ICA signals', 2, myDim)
    figH    = gobjects(1, 8);
    bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    plotPerson = 1;
    plotScan = 1;
    for iFig = 1:64
        figH(iFig) = figure('WindowStyle', 'docked', ...
            'Name', sprintf('Figure %d', iFig), 'NumberTitle', 'off');
        drawnow;
        %     pause(0.02);  % Magic, reduces rendering errors
        set(get(handle(figH(iFig)), 'javaframe'), 'GroupName', 'fMRI ICA signals');
        
        % Plotting data
        for i = 1:20
            subplot(5,4,i)
            plot(fMriPerson{plotPerson}.Scan{plotScan}(:,i))
            grid on
            xlabel('Time [s]')
            ylabel(['Signal ',num2str(i)])
        end
        suptitle(['Timeseries signals - Person ', num2str(plotPerson) ,' scan ', num2str(plotScan)])
        
        if plotScan ~= 4
            plotScan = plotScan + 1;
        else
            plotScan = 1;
            plotPerson = plotPerson + 1;
        end
    end
    warning(bakWarn);
    clear plotScan plotPerson;
end

%% Run SLS algorithm
slsSettings.toolboxLocation = 'C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\2. Harm rank reduced noise algorithm\toolbox';
addpath(slsSettings.toolboxLocation);
% slsSettings.sampleTime = 1.986622;
slsSettings.sampleTime = 2; % debug
addpath('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\2. Harm rank reduced noise algorithm\simulations'); % debug

% Interpolate data to get rid of ARX error (didnt work!)
interpFactor = 1;
for i = 1:20
    interpData(:,i) = interp(fMriPerson{1}.Scan{1}(:,i),interpFactor);
end

% Convert timeseries data to iddata object
% data = iddata(fMriPerson{1}.Scan{1},[],slsSettings.sampleTime);
% data = iddata(fMriPerson{1}.Scan{1},fMriPerson{1}.Scan{1},slsSettings.sampleTime);
data = iddata(interpData,[],slsSettings.sampleTime/interpFactor);

% Simulation patameters
% N = 300; % number of datapoints 
% var_e = 0.3; % power of noise
% Ts = slsSettings.sampleTime; % sample time

L = 20; % number of nodes
K = 20; % number of external excitations
original_order = 2; % order of original network modules

% addpath(genpath('../toolbox')) % access the path with the functions


% Generate network
% while(1)
%     [G0,R0,H0,T0,unstable] = generate_physical_network(L,K,Ts,1);
%     if unstable + sum(abs(eig(G0)) > 0.95) == 0 && norm(G0) > 0
%         break;
%     end
% end

% figure()
% bodemag(G0)
% grid on


% Original network topology
% G0_imp = impulse(G0);
% original_structure.NG = squeeze(G0_imp(2,:,:)) ~= 0; 

% R0_imp = impulse(R0);
% original_structure.NR = squeeze(R0_imp(1,:,:)) ~= 0; 

% H0_imp = impulse(H0);
% original_structure.NH = squeeze(H0_imp(2,:,:)) ~= 0; 

% Fixed dynamics
% fixed.R = R0;


% orders.D = original_order;
% orders.NG = original_order;
% orders.NR = original_order;
% orders.NH = original_order;

% Generate data
% [w,r,e] = generate_data(T0,N,L,K,var_e);


% make the true noise model monic, and adjust the noise source accordingly
% for ii = 1:1:L
%     scale = H0.num{ii,ii}(2);
%     H0.num{ii,ii} = [H0.num{ii,ii}(2:end) 0] ./ scale;
%     e(ii,:) = [0 e(ii,1:end-1)] .* scale;
% end

% Create idnetw structure
% slsSettings.topo = 0;   % TODO
% slsSettings.n = 0;      % TODO
% slsSettings.fixed = 0;  % TODO
% slsSettings.netw = idnetw();
% slsSettings.netw = idnetw(original_structure,orders,fixed);
slsSettings.netw = idnetw('diagR',L,K,2);

% Set slsOptions
slsSettings.opt = slsOptions();

% Run algorithm
slsOutput = SLS(data, slsSettings.netw ,slsSettings.opt)
% slsOutput = SLS(data);




