%%% Script to launch a 2D tracking of particles on micrographies clear all close all %%% Architecture and experiment parameters 
session.input_path='/home/ybs8124/Documents/chips3.2/Clog_I0.1/'; 
session.output_path='/home/ybs8124/Documents/chips3.2/ProcessedDATA/'; 

expdate='251125';
DP=30;   
I=0.1; 

manipName=['Clogging_DP' num2str(DP) '_I' num2str(I) '_' expdate]; 
acquisName=['unclogging_DP' num2str(DP) '_I' num2str(I) '_' expdate]; % do not forget final /
pictureName=acquisName;

%%% Calibraton parameters 

period=50; %µm, between two pores (center-to-center) 
size_pore=5; %µm, pore width 
extension='png'; %optional 

%%% Background parameters 
Step=20; 
format='%04d.png'; %optional, default '%06d.png' 

%%% Centenr findings parameters 
th=30; %threshold for particle detection 
sz=10; %typical size of the particles 
nframes=[]; %optional, total number of frames 
Test_CF=false; %optional, true-> test mode, false (default)-> classic mode 
PartialSave=0; %if PartialSave>0 it removes background %%% for the first PartialSave frames and save them in %%% folderout/TestThreshold. Can be usefull to check if %%% a particle moves
BackgroundType="BackgroundMean"; %optional, determine which background is substracted to pictures. By defaut is equal to BackgroundMean %%% 
FileName=['centers_' pictureName]; %Name of matlab file containing centers without extension 
maxdist=100; % : maximum travelled distance between two successive frames (px) 
lmin=500; % : minimum length of a trajectory (number of frames) 
flag_pred=1; % : 1 for predictive tracking, 0 otherwise 
npriormax=5; % : maximum number of prior frames used for predictive tracking 
flag_conf=1; % : 1 for conflict solving, 0 otherwise 
minFrame=1; % : (optional) number of the first frame. Default = 1. 
Test_TR=false; % : (optional) allows you to not save data when you are doing %tests to find best parameters 

%%% Stitiching parameters 
dfmax=500; % maximum number of tolerated missing frames to reconnect to trajectories 
dxmax=2; % maximum tolerated distance (in norm) between projected % point after the first trajectory and the real beginning position of the % stitched one, 
dvmax=2000; % maximum tolerated relative velocity difference between 

lmin_stitch=100; % minimum length for a trajectory to be stitched 
h5filename=['tracks_centers_' pictureName]; % Name of the tracks h5 filename without extension 

%%% Plot sequence=false; %tracé ou dont de la séquence d'images pour illustration 
%onlytracks=false; 

%% Calibration - if not already done 
calib = MakeCalibration(session, manipName, acquisName, period, size_pore, extension);

%% Background - if not already done 
close all 
BackgroundComputation(session,manipName,acquisName,pictureName,Step,format, extension)
%% Center findings 
CC = CenterFinding2D(session,manipName,acquisName,pictureName,th,sz,nframes,Test_CF,PartialSave,BackgroundType,format,extension)

%% Tracking 
[tracks,traj]=track2d(session,manipName,acquisName,FileName,maxdist,lmin,flag_pred,npriormax,flag_conf,Test_TR); 

%% Stitching 
StitchedTraj = Stitching(session,manipName,acquisName,pictureName,h5filename,dfmax,dxmax,dvmax,lmin_stitch);
%% Stitched file reading 
filepath=fullfile(session.output_path, manipName, acquisName, ['stitched_dfmax' num2str(dfmax) '_' pictureName]);

% path of the stitching h5 file without extension 
[stitchtraj, stitchingparameters] = h52stitch(filepath); 
%% Plot trajectories without stitching 
close all
onlytracks=true; %plot_trajectoires_fromh5 %(à lancer en autonomie après avoir chargé les % paramètres plus haut) 
addpath(genpath('/home/ybs8124/Documents/brewermap'));
plot_trajectoires_fromh5 % ou à lancer dans la foulée (attention génère une séquence d'images) 

% %% Plot trajectories with stitching onlytracks=false; %plot_trajectoires_fromh5 %(à lancer en autonomie après avoir chargé les % paramètres plus haut) 
% % plot_trajectoires_fromh5 % ou à lancer
%%
addpath(genpath('/home/ybs8124/Documents/Tracking_PTV2D_1D_IMFT/vitesse'))
OUT = vitesse(tracks_unstuck_filtered,calib)

%%
options.smoothWindow = 5;
options.nbins = 200;
analysis=analyzeTrajectory(tracks_unstuck_filtered,calib,options);
%% Lagrangian
[OUT, tracks_unstuck_filtered] = vitesse(tracks_unstuck_filtered, calib);

run('run_lagrangian_statistics.m');
%% Moment scaling spectrum

addpath(genpath('/home/ybs8124/Documents/Tracking_PTV2D_1D_IMFT/DC-MSS/DC-MSS-master/Code'));  % dossier où tu as mis tous les fichiers du dépôt

traj_sub = traj(:, :, 1:10:end);  % exemple
trackMatrix_sub = makeTrackMatrixForDCMSS_raw(traj_sub);
[trackClass, mssSlope, genDiffCoef, scalingPower, normDiffCoef, estimError] = trackMSSAnalysis(trackMatrix_sub, 2);
