function [tracks, traj] = track2d(session,manipName,acquisName,FileName,maxdist,lmin,flag_pred,npriormax,flag_conf,test)
% Load matches data and give it to track3d_manualfit to track particles.
% ____________________________________________________________________________
% INPUTS
% session   : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
% ManipName : Name of the experiment
% FileName  : Name of matlab file containing centers without extension
% maxdist   : maximum travelled distance between two successive frames
% lmin      : minimum length of a trajectory (number of frames)
% flag_pred : 1 for predictive tracking, 0 otherwise
% npriormax : maximum number of prior frames used for predictive tracking
% flag_conf : 1 for conflict solving, 0 otherwise
% minFrame : (optional) number of the first frame. Default = 1.
% test      : (optional) allows you to not save data when you are doing
% tests to find best parameters
%
% OUTPUTS
% traj(kt).ntraj  : trajectory index
% traj(kt).L      : trajectory length
% traj(kt).frames : trajectory frames
% traj(kt).x      : x-position
% traj(kt).y      : y-position
% traj(kt).nmatch : element indices in tracks
% tracks          : trajectory raw data
%
% 2020-2021 D. Dumont (adapted from M. Bourgoin)
% adapted for IMFT/MPB by O. Liot march 2023
% ____________________________________________________________________________

if ~exist('test','var') || isempty(test)
    test = false;
end

folder = fullfile(session.output_path, manipName, acquisName);
if ~isfolder(folder)
    mkdir(folder);
end

filename = fullfile(folder, [FileName '.mat']);
if ~exist(filename, 'file')
    error('Centers file not found: %s', filename);
end

% disp('Loading centers...');
load(filename, 'CC');
%data = h52matches(filename,NbFrame,minFrame);
for i = 1:numel(CC)
    if ~isfield(CC(i), 'X')
        CC(i).X = [];
    end
    if ~isfield(CC(i), 'Y')
        CC(i).Y = [];
    end
    if ~isfield(CC(i), 'frame')
        CC(i).frame = [];
    end
end

%% Call for track2d_manualfit function
[tracks, traj] = track2d_manualfit(folder, FileName, CC, maxdist, lmin, flag_pred, npriormax, flag_conf, test);

% %% Call for track2d_polyfit function
% [tracks,traj] = track2d_polyfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf);
end
