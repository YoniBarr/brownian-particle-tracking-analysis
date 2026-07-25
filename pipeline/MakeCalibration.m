function [calib] = MakeCalibration(session,manipName,acquisName,period,size_pore,extension)
% Make calibration file.
% With human help, detect the region of interest where tracking will be
% performed
% Store everything in a structure calib.mat saved in the same directory
%----------------------------------------------------------------------------------------
% INPUT Parameters:
%   dirIn                 : directory containning reference picture. calib.mat is saved in this directory,
%   period                : distance between two consecutive pores in µm
%   size_pore               : pore width in µm
%   extension (optional)  : pictures extension. By defaut extention = 'png', 
% 
% OUTPUT
%     a calib.mat file saved in dirIn which contains the structure 'calib' :
%     
%     calib.pix2um                   : pixel to µm calibration
%     calib.fps                   : framerate in fps
%     calib.orientation   : orientation of the flow (d for downwards, u for upwards, r for rightwards, l for leftwards)
%     calib.globalROI            : rectangle [xmin ymin width height] for the global ROI
%     calib.localROI(k) : local ROI above each pore numbered k

% ------------------------------------------------------------------------------------------

%% Definition of optional variables it they do not exist.
if ~exist('extension','var')
    extension = 'png';
end

dirIn=session.input_path;
dirOut=session.output_path;

%% Ask framerate
calib.fps=input('Framerate of the acquisition (fps)? ');

%% Load and display bright field picture
refFolder=fullfile(dirIn,manipName,acquisName);
refFileName=dir(fullfile(refFolder,['*.' extension]));

if isempty(refFileName)
    error('No image found in %s matching *.%s',refFolder,extension);
end

% Sort by name
[~, idx]=sort({refFileName.name});
refFileName=refFileName(idx);

% Always take the 2nd image
if numel(refFileName)<2
    warning('Only one image found; using the first one as reference.');
    refFileName=refFileName(1);
else
    refFileName=refFileName(2);
end

fprintf('Using reference image: %s\n',refFileName.name);

img_ref = imread(fullfile(refFolder,refFileName.name));
figure(1); clf;
imshow(imadjust(img_ref));
title('Bright field reference picture (2nd image)');

%% Ask number of pores and orientation
nb_nc=input('How many pores ? ');

rotation=input('Flow orientation ?\n (d for downwards, u for upwards, r for rightwards, l for leftwards)','s');
calib.orientation=rotation;

figure(2); clf; 
axis on;
if rotation=='r'
    img_rotate=imrotate(img_ref,-90);
    imshow(imadjust(img_rotate))
elseif rotation=='d'
    img_rotate=img_ref;
    imshow(imadjust(img_rotate))
elseif rotation=='u'
    img_rotate=imrotate(img_ref,-180);
    imshow(imadjust(img_rotate))
elseif rotation=='l'
    img_rotate=imrotate(img_ref,90);
    imshow(imadjust(img_rotate))
end
title('Picture rotated')

%% Select ROI
figure(2)
disp('Click on the top-left corner of the left pore');
[x1, y1]=ginput(1);
disp('Click on the top-left corner of the right pore');
[x2, y2]=ginput(1);

%% Compute pixel to µm calibration

calib.pix2um =(nb_nc-1)*period/abs(x1-x2);

%% Global ROI

xmin=1;width=size(img_rotate, 2);
ymin=1;height=max(y1, y2);
rect=[xmin ymin width height];
calib.globalROI=rect;

Icropped = imcrop(img_rotate, rect);
figure(3); clf;
imshow(imadjust(Icropped));
title('Global ROI');

%% ROI for each pore

for i=1:nb_nc
    x=x1+size_pore/calib.pix2um/2-period/calib.pix2um/2+(i-1)*period/calib.pix2um;
    height_local=y1+(i-1)*period/calib.pix2um*(y2-y1)/(x2-x1);
    calib.localROI{i}=[x,1,period/calib.pix2um,height_local];
end

figure(4); clf;
imshow(imadjust(Icropped));

hold on;
for i=1:nb_nc
rectangle('Position', calib.localROI{i}, 'EdgeColor', 'r');
end

%% save calib
outDir=fullfile(dirOut, manipName);
if ~isfolder(outDir)
    mkdir(outDir);
end

save(fullfile(outDir, 'calib.mat'), 'calib');


end
