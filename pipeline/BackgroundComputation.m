%% 
function BackgroundComputation(session, manipName, acquisName, pictureName, Step, format, extension)
%%% Compute the background for the pictures of the camera NumCam in the
%%% experiment manipName
%%% Between picture StartFrame and EndFrame ( every Step frame) it takes
%%% the maximal/minimal/mean intensity for each pixel.
%----------------------------------------------------------------------------
%%% Parameters :
%%%     session      : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%%%     manipName    : Name of the folder experiment
%%%     acquisName    : Name of the folder acquisition
%%%     Step (optional)         : step between 2 frames taken for computation
%%%     format (optional)     : picture names. By defaut it is '%06d.png'.

%------------------------------------------------------------------------------
% 2020-2021 : D. Dumont (adapted from M. Bourgoin)
% adapted for IMFT/MPB by O. Liot march 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% By defaut format='%05d.tif'-------------------------------

if nargin < 7
    extension = 'png';
end

%% Definitions of folder
folderIn  = fullfile(session.input_path, manipName, acquisName);
folderOut = fullfile(session.output_path, manipName, acquisName);

%% Creation of output Folder if it does not exist
if ~isfolder(folderOut)
    mkdir(folderOut);
end

%% Images used for background
% If Step is not defined, take a step such that the average will be done over ~1000 pictures.

% List all images in the acquisition folder 
fileList = dir(fullfile(folderIn, [pictureName '_*.' extension]));
if isempty(fileList)
    error(' No images found in %s matching pattern %s', folderIn, [pictureName '_*.' extension]);
end
fprintf(' Found %d images in %s\n', numel(fileList), folderIn);

% Sort by filename (important to keep correct chronological order)
[~, idx] = sort({fileList.name});
fileList = fileList(idx);

% Determine Step automatically if not provided
if ~exist('Step', 'var')
    if numel(fileList) > 1000
        Step = floor(numel(fileList) / 1000);
    else
        Step = 1;
    end
end
fprintf(' Computing background using every %d frames...\n', Step);

%Read the first image to initialize arrays
im0 = imread(fullfile(folderIn, fileList(1).name));
imgClass = class(im0);
imSize = size(im0);
fprintf('Image type: %s | Size: [%d x %d]\n', imgClass, imSize(1), imSize(2));

% Initialize background matrices
switch imgClass
    case 'uint8'
        BackgroundMean = zeros(imSize, 'double');
        BackgroundMax  = zeros(imSize, 'double');
        BackgroundMin  = ones(imSize, 'double') * 255;
    case 'uint16'
        BackgroundMean = zeros(imSize, 'double');
        BackgroundMax  = zeros(imSize, 'double');
        BackgroundMin  = ones(imSize, 'double') * 65535;
    otherwise
        error(' Unsupported image class: %s (must be uint8 or uint16)', imgClass);
end

%Compute mean, max, and min backgrounds
count = 0;
for i = 1:Step:numel(fileList)
    img = imread(fullfile(folderIn, fileList(i).name));
    img = double(img);
    BackgroundMean = BackgroundMean + img;
    BackgroundMax  = max(BackgroundMax, img);
    BackgroundMin  = min(BackgroundMin, img);
    count = count + 1;
end

BackgroundMean = BackgroundMean / count;

% Convert back to the same data type as input images
switch imgClass
    case 'uint8'
        BackgroundMean = uint8(BackgroundMean);
        BackgroundMax  = uint8(BackgroundMax);
        BackgroundMin  = uint8(BackgroundMin);
    case 'uint16'
        BackgroundMean = uint16(BackgroundMean);
        BackgroundMax  = uint16(BackgroundMax);
        BackgroundMin  = uint16(BackgroundMin);
end

% Save mean background image 
backgroundName = sprintf('Background_%s.%s', pictureName, extension);
imwrite(BackgroundMean, fullfile(folderIn, backgroundName));
fprintf(' Background image saved to %s\n', fullfile(folderIn, backgroundName));

% Save all background matrices as .mat file
BackgroundFile = fullfile(folderOut, ['Background_' acquisName '.mat']);
save(BackgroundFile, 'BackgroundMean', 'BackgroundMax', 'BackgroundMin', '-v7.3');
fprintf(' Background MAT saved to %s\n', BackgroundFile);
