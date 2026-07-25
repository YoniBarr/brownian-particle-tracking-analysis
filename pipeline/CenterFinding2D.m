function CC = CenterFinding2D(session,ManipName,acquisName,pictureName,th,sz,nframes,Test,PartialSave,BackgroundType,format,extension)
%%% Detect particles position in picture and provides their positions in 
%%% px.
%--------------------------------------------------------------------------------
%%% Parameters :
%%%     session                    : Path to the achitecture root (2 fields: session.input_path 
% and session.output_path)
%%%     ManipName                  : Name of the folder experiment
%%%     nframes  (optional)                  : total number of pictures
%%%     acquisName                 : Acquisition folder name
%%%     pictureName                : Prefix of picture files
%%%     th                         : Detection threshold
%%%     sz                         : typical size of the particles
%%%     Test                       : true-> test mode, false-> classic mode
%%%     PartialSave (optional)     : Number of frames to save with background removed
%%%     BackgroundType (optional)  : Background type ('BackgroundMean', 'BackgroundMax', 'BackgroundMin')
%%%     format (optional)          : Picture name format (e.g. '%06d.png')
%%%     extension (optional)       : File extension ('png', 'tif', etc.)
%%%     for the first PartialSave frames and save them in
%%%     folderout/TestThreshold. Can be usefull to check if a particle
%%%     moves.
%--------------------------------------------------------------------------------
% 2020-2021 : D. Dumont (adapted from M. Bourgoin)
% adapted for IMFT/MPB by O. Liot (march 2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%% Test if Test exists or not
if ~exist('Test','var'), 
    Test=false; 
end
%% Test if BackgroundType exists or not
if ~exist('BackgroundType','var'), 
    BackgroundType="BackgroundMean"; 
end

% By defaut format='%06d.png'
if ~exist('format','var'), 
    format = '%04d.png'; 
end
nb_0 = str2double(format(3)); 

%% Test if PartialSave exists or not
if ~exist('PartialSave','var'), PartialSave = 0; end
if ~exist('extension','var'), extension = 'png'; end


BaseName=[pictureName '_' format];


%% Definition of folders

folderin = fullfile(session.input_path, ManipName, acquisName);
folderout = fullfile(session.output_path, ManipName, acquisName);
if ~isfolder(folderout)
    mkdir(folderout);
end

fprintf('Processing acquisition: %s\n', acquisName);

BackgroundFile = fullfile(folderout, ['Background_' acquisName '.mat']);
CalibFile = fullfile(session.output_path, ManipName, 'calib.mat');

if ~isfile(BackgroundFile)
    error('Background file not found: %s', BackgroundFile);
end
if ~isfile(CalibFile)
    error('Calibration file not found: %s', CalibFile);
end
%% Load background and calib
load(BackgroundFile, 'BackgroundMin', 'BackgroundMax', 'BackgroundMean');
load(CalibFile, 'calib');

%% Determine available frames

fileList = dir(fullfile(folderin, [pictureName '_*.' extension]));

if isempty(fileList)
    error('No images found in %s matching %s', folderin, fullfile(folderin, [pictureName '_*.' extension]));
end

[~, idx] = sort({fileList.name});
fileList = fileList(idx);

% Extract numeric frame indices
frames = [];
for i = 1:numel(fileList)
    tokens = regexp(fileList(i).name, '\d+', 'match');
    if ~isempty(tokens)
        frames(end+1) = str2double(tokens{end});
    end
end
firstFrame = min(frames);
lastFrame = max(frames);

if isempty(nframes) || nframes > lastFrame
    nframes = lastFrame;
end


%% Choice of background type
if BackgroundType=="BackgroundMean"
    Background=BackgroundMean;
elseif BackgroundType=="BackgroundMax"
    Background=BackgroundMax;
elseif BackgroundType=="BackgroundMin"
    Background=BackgroundMin;
end

%% Find centers

if ~Test
    CC = struct();
    if firstFrame < 1, firstFrame = 1; end

    for kframe = firstFrame:nframes
        fprintf('Processing frame %d/%d\n', kframe, nframes);
        ImgName = fullfile(folderin, sprintf(BaseName, kframe));

        if exist(ImgName, 'file')
            Iraw= imread(ImgName);
            Im=imsubtract(Iraw, Background);
            Im(Im<0)=0;


        if calib.orientation=='r'
                Im=imrotate(Im,-90);
            elseif calib.orientation=='d'
                Im=Im;
            elseif calib.orientation=='u'
                Im=imrotate(Im,-180);
            elseif calib.orientation=='l'
                Im=imrotate(Im,90);
        end


            Nx = size(Im,2);
            Ny = size(Im,1);
            out = pkfnd(Im, th, sz); % Provides intensity maxima positions
            fprintf("Frame %d → detected %d particles\n", i, size(out,1));
            npar = size(out,1);

            %% We keep only spots with a gaussian shape
            x = []; y = [];
            for j = 1:npar
                Nwidth = 1;
                if (out(j,2)-Nwidth > 0) && (out(j,1)-Nwidth > 0) &&(out(j,2)+Nwidth < Ny) && (out(j,1)+Nwidth < Nx)
                    Ip = double(Im(out(j,2)-Nwidth:out(j,2)+Nwidth,out(j,1)-Nwidth:out(j,1)+Nwidth));
                    x(end+1)=out(j,1)+0.5*log(Ip(2,3)/Ip(2,1))/log((Ip(2,2)^2)/(Ip(2,1)*Ip(2,3)));
                    y(end+1)=out(j,2)+0.5*log(Ip(3,2)/Ip(1,2))/log((Ip(2,2)^2)/(Ip(1,2)*Ip(3,2)));

                end
            end

            if isempty(x)
    CC(kframe).X = [];
    CC(kframe).Y = [];
    CC(kframe).frame = [];
else
    CC(kframe).X = x;
    CC(kframe).Y = y;
    CC(kframe).frame = kframe * ones(1, numel(x));
end


        else
            disp('No Picture')
           CC(kframe).frame = kframe;
        end
    end

    %% Centers saving into a .mat file
    save(fullfile(folderout, ['centers_' pictureName '.mat']), 'CC', 'nframes', '-v7.3');
    % fprintf(CC(1200).X);
    fprintf('Centers saved successfully in %s\n', folderout);

end
end

