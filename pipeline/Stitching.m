function StitchedTraj = Stitching(session,ManipName,acquisName,pictureName,h5filename,dfmax,dxmax,dvmax,lmin)
% Create automatic filepath, load data and call for stitchTracks function.
% To use after track3d.m function.
% 04/2020 - David Dumohhhyhhnt
%----------------------------------------------------------------------------------------
% Parameters:
%   session        : session.path contains MyPath, (2 fields: session.input_path
% and session.output_path),
%   ManipName      : Name of the experiment,
%   acquisName       : Name of the acquisition,
%   h5fileame       : Name of the tracks h5 filename without extension,
%   dfmax          : maximum number of tolerated missing frames to
%   reconnect to trajectories,
%   dxmax          : maximum tolerated distance (in norm) between projected
%   point after the first trajectory and the real beginning position of the
%   stitched one,
%   dvmax          : maximum tolerated relative velocity difference between
%   the end of the first trajectory and the beginning of the stitched one,
%   lmin           : minimum length for a trajectory to be stitched,
% window_length : size of the window used for savgol filtering (to compute speed),
% poly_order    : order of the polynom used for savgol filtering (to
% compute velocity),
% framerate     : framerate of the experiment.
% ------------------------------------------------------------------------------------------

datapath = fullfile(session.output_path,ManipName,acquisName,h5filename);
filepath = fullfile(session.output_path, ManipName,acquisName, ['stitched_dfmax' num2str(dfmax) '_' pictureName]);

%% Data loading
traj = h52tracks(datapath);

fprintf("Before stitching, we have %d trajectories\n",numel(traj))

%% Let's stitch trajectories !
StitchedTraj = stitchTracks(traj,dfmax,filepath,dxmax,dvmax,lmin);

end



