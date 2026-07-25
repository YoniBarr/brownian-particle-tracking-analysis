function trackMatrix = makeTrackMatrixForDCMSS_raw(input)
% Build DC-MSS trackMatrix WITHOUT any preprocessing
% INPUT:
%   input : traj struct (from track2d_manualfit) OR tracks struct (from h52tracks)
%           OR string path to .h5 file (with or without .h5)
% OUTPUT:
%   trackMatrix : Ntracks x (Nframes*8)  (empty [] if no valid trajectories)

% Load / convert input into traj struct array with fields: frames,x,y,L 
traj = [];

if ischar(input) || isstring(input)
    fname = char(input);
    if ~endsWith(fname, '.h5')
        fname = [fname '.h5'];
    end
    if ~isfile(fname)
        error('File not found: %s', fname);
    end
    % try to use provided loader h52tracks if available (returns tracks struct)
    try
        [tracks_struct,~] = h52tracks(fname);
    catch
        % fallback: try to read using h5read and rebuild tracks_struct
        info = h5info(fname);
        % attempt to reconstruct simple tracks_struct with fields frames,x,y,L,ntraj
        % This fallback assumes that h5 file stores per-field concatenated arrays and /L and /ntraj
        try
            L = h5read(fname,'/L');
            ntraj = h5read(fname,'/ntraj');
        catch
            error('Cannot read .h5 with h52tracks and fallback failed. Ensure file format matches.');
        end
        % try to read possible fields
        possibleFields = {'frames','x','y'};
        dataMap = struct();
        for f = 1:numel(possibleFields)
            nm = possibleFields{f};
            try
                dataMap.(nm) = h5read(fname, ['/' nm]);
            catch
                dataMap.(nm) = [];
            end
        end
        % Reconstruct tracks_struct
        tracks_struct = struct([]);
        c = 1;
        for kt = 1:numel(L)
            temp = L(kt);
            if temp <= 0
                continue;
            end
            tracks_struct(end+1).frames = dataMap.frames(c:c+temp-1);
            tracks_struct(end).x = dataMap.x(c:c+temp-1);
            tracks_struct(end).y = dataMap.y(c:c+temp-1);
            tracks_struct(end).L = temp;
            c = c + temp;
        end
    end

    % convert tracks_struct - traj
    if isempty(tracks_struct)
        trackMatrix = [];
        return;
    end
    for k = 1:numel(tracks_struct)
        % assume tracks_struct(k) has fields frames,x,y or similar
        if isfield(tracks_struct(k),'frames')
            fr = tracks_struct(k).frames;
        else
            fr = [];
        end
        if isfield(tracks_struct(k),'x')
            xx = tracks_struct(k).x;
        else
            xx = [];
        end
        if isfield(tracks_struct(k),'y')
            yy = tracks_struct(k).y;
        else
            yy = [];
        end
        traj(k).frames = fr;
        traj(k).x = xx;
        traj(k).y = yy;
        if isfield(tracks_struct(k),'L')
            traj(k).L = tracks_struct(k).L;
        else
            traj(k).L = numel(fr);
        end
    end

elseif isstruct(input)
    % input may already be traj or tracks
    s = input;
    % detect if s is an array where each element already has frames,x,y
    if isfield(s, 'frames') && isfield(s,'x') && isfield(s,'y')
        traj = s;
    else
        % maybe tracks style with fields like L,x,y,frames
        % try to map each element
        for k = 1:numel(s)
            fr = [];
            xx = [];
            yy = [];
            if isfield(s(k),'frames'), fr = s(k).frames; end
            if isfield(s(k),'x'), xx = s(k).x; end
            if isfield(s(k),'y'), yy = s(k).y; end
            traj(k).frames = fr;
            traj(k).x = xx;
            traj(k).y = yy;
            if isfield(s(k),'L'), traj(k).L = s(k).L; else traj(k).L = numel(fr); end
        end
    end
else
    error('Unsupported input type. Provide traj struct, tracks struct, or .h5 filename.');
end

% Normalize fields: force row vectors, handle empties
for k = 1:numel(traj)
    if ~isfield(traj(k),'frames') || isempty(traj(k).frames)
        traj(k).frames = [];
        traj(k).x = [];
        traj(k).y = [];
        traj(k).L = 0;
    else
        traj(k).frames = double(traj(k).frames(:)');  % row vector
        traj(k).x = double(traj(k).x(:)');            % row vector
        traj(k).y = double(traj(k).y(:)');            % row vector
        traj(k).L = numel(traj(k).frames);
    end
end

% Remove empty trajectories only (technical safety)
valid = arrayfun(@(t)t.L > 0, traj);
traj = traj(valid);

if isempty(traj)
    trackMatrix = [];
    return;
end

% Robust computation of max frame (avoid cell2mat pitfalls)
maxFrame = 0;
for k = 1:numel(traj)
    if ~isempty(traj(k).frames)
        maxFrame = max(maxFrame, max(traj(k).frames));
    end
end
if maxFrame <= 0
    trackMatrix = [];
    return;
end
Nframes = maxFrame;
Ntracks = numel(traj);

% Allocate and fill trackMatrix
trackMatrix = zeros(Ntracks, Nframes * 8);

for k = 1:Ntracks
    fr = traj(k).frames;
    x  = traj(k).x;
    y  = traj(k).y;
    if isempty(fr)
        continue;
    end
    % ensure lengths match
    n = min([numel(fr), numel(x), numel(y)]);
    if n < numel(fr)
        % if mismatch, trim to common length and warn
        warning('Trajectory %d: frames/x/y length mismatch, trimming to %d entries.', k, n);
        fr = fr(1:n); x = x(1:n); y = y(1:n);
    end
    for j = 1:numel(fr)
        f = fr(j);
        if f < 1 || f > Nframes
            continue;
        end
        col = (f-1)*8 + 1;
        trackMatrix(k, col)   = x(j);
        trackMatrix(k, col+1) = y(j);
    end
end

end
