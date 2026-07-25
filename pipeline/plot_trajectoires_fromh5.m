%% Film détection pour particules trackées proprement sur trajectoire
% A lancer en autonomie après avoir chargé les paramètres dans
% submission_fulltracking.m

%clear all
close all
set(0,'DefaultAxesColorOrder',brewermap(NaN,'RdYlBu'))
sequence=true
%% données manip
ProcessedDataDirectory=session.output_path; %Dossier postprocessing
DataDirectory=session.input_path; %Dossier initiales (images)

acquisName=acquisName; %Dossier travail

picturesPrefix=acquisName(1:end-1);

refFilepath=fullfile(session.output_path, manipName, 'calib.mat');

tracksDirectory=fullfile(ProcessedDataDirectory, manipName, acquisName);
pictureDirectory=fullfile(DataDirectory, manipName, acquisName);

extensionFrame_movie='.png';

%% Chargement trajectoires et calibration

load(refFilepath)

%onlytracks=false;
if ~onlytracks
    dfmax=num2str(100);   
    tracksPath=fullfile(tracksDirectory, ['stitched_dfmax' dfmax '_' pictureName]); % sans le .h5
    tracks=h52stitch(tracksPath);

else
    tracksPath = fullfile(tracksDirectory, ['tracks_centers_' pictureName]); %fichier avec trajectoires, sans le .h5
    tracks=h52tracks(tracksPath);
end

%% Calcul nombre de frames

nb_0=4;%nombre de 0 dans le numéro de l'image
files=dir(fullfile(pictureDirectory, [picturesPrefix '*' extensionFrame_movie]));
PicturesList={files.name}';
frames=cellfun(@(f) str2double(f(end-(3+nb_0):end-4)), PicturesList);
firstFrame=frames(1);
if ~exist('nFrames','var')
    nframes=frames(end);
end

%% suppression des morceaux de trajectoires immobiles

% remove stopped particles
disp('Remove stuck particles...')
tic
limiteur_dist=2;
if exist('tracks_unstuck','var')
    clear tracks_unstuck
end
for ii=1:numel(tracks)
    tracks_unstuck(1,ii) = struct('L',[],'x',[],'y',[],'frames',[]);
end
for i=1:numel(tracks)
    
    if (calib.pix2um)*sqrt((tracks(i).y(1)-tracks(i).y(end)).^2+(tracks(i).x(1)-tracks(i).x(end)).^2)>limiteur_dist
        tracks_unstuck(1,i).y=tracks(i).y;
        tracks_unstuck(1,i).x=tracks(i).x;
        tracks_unstuck(1,i).frames=tracks(i).frames;
        tracks_unstuck(1,i).L=numel(tracks_unstuck(1,i).x);
    else
        tracks_unstuck(1,i).y=[];
        tracks_unstuck(1,i).x=[];
        tracks_unstuck(1,i).T=[];
        tracks_unstuck(1,i).L=0;
    end

end

tracks_unstuck([tracks_unstuck.L]==0)=[];
disp('...finished')
toc

%% conservation des trajectoires longues et qui démarrent au début

% remove stopped particles
disp('Keep long trajectories starting at the unclogging start...')
tic
long_min=1500;
if exist('tracks_unstuck_long_begin','var')
    clear tracks_unstuck_long_begin
end
for ii = 1:numel(tracks_unstuck)
    tracks_unstuck_long_begin(1,ii) = struct('L',[],'x',[],'y',[],'frames',[]);
end
for i=1:numel(tracks_unstuck)
    %if abs(nanmean(calibration/dt*diff([tracks(i).X])))>limiteur
    % On choisit les trajectoires qui ont réellement bougé
    if tracks_unstuck(i).L>long_min %&& tracks_unstuck(i).frames(1)==0
        tracks_unstuck_long_begin(1,i).y=tracks_unstuck(i).y;
        tracks_unstuck_long_begin(1,i).x=tracks_unstuck(i).x;
        tracks_unstuck_long_begin(1,i).frames=tracks_unstuck(i).frames;
        %tracks_nostop(1,i).L=tracks(i).L;
        tracks_unstuck_long_begin(1,i).L=numel(tracks_unstuck_long_begin(1,i).x);
    else
        tracks_unstuck_long_begin(1,i).y=[];
        tracks_unstuck_long_begin(1,i).x=[];
        tracks_unstuck_long_begin(1,i).T=[];
        tracks_unstuck_long_begin(1,i).L=0;
    end

end

tracks_unstuck_long_begin([tracks_unstuck_long_begin.L]==0)=[];
disp('...finished')
toc

%% filtrage savgol des trajectoires uncut (robuste)
disp('Filtering unstuck trajectories starts...')
tic

tracks_unstuck_filtered = repmat(struct('L',[],'x',[],'y',[],'r',[],'frames',[]), 1, numel(tracks_unstuck));

for i = 1:numel(tracks_unstuck)
    x = tracks_unstuck(i).x;
    y = tracks_unstuck(i).y;
    n = numel(x);

    if n < 3
        % trop court pour filtrer
        tracks_unstuck_filtered(i).x = x;
        tracks_unstuck_filtered(i).y = y;
    else
        % choisir framelen : impair, <=101 et <= n
        framelen = min(101, n);
        if mod(framelen,2) == 0
            framelen = framelen - 1;
        end
        if framelen < 3
            % fallback minimal
            tracks_unstuck_filtered(i).x = x;
            tracks_unstuck_filtered(i).y = y;
        else
            % appliquer sgolayfilt
            try
                tracks_unstuck_filtered(i).y = sgolayfilt(y, 1, framelen);
                tracks_unstuck_filtered(i).x = sgolayfilt(x, 1, framelen);
            catch ME
                warning('sgolayfilt failed for track %d (n=%d, framelen=%d). Copy raw. Error: %s', i, n, framelen, ME.message);
                tracks_unstuck_filtered(i).x = x;
                tracks_unstuck_filtered(i).y = y;
            end
        end
    end

    tracks_unstuck_filtered(i).frames = tracks_unstuck(i).frames;
    tracks_unstuck_filtered(i).L = numel(tracks_unstuck_filtered(i).x);
end

% nettoyer
tracks_unstuck_filtered([tracks_unstuck_filtered.L] == 0) = [];
disp('...finished')
toc


%% Tracé trajectoires complètes


%close all
lastPicture_temp = imread(fullfile(pictureDirectory, PicturesList{end}));
if calib.orientation=='r'
    lastPicture=imrotate(lastPicture_temp,-90);
elseif calib.orientation=='d'
    lastPicture=lastPicture_temp;
elseif calib.orientation=='u'
    lastPicture=imrotate(lastPicture_temp,-180);
elseif calib.orientation=='l'
    lastPicture=imrotate(lastPicture_temp,90);
end
%figure
%clf

subplot(2,2,1);
imshow(lastPicture*4)
axis on
hold on
title('Full tracks')

for j=1:numel(tracks)
    plot(tracks(j).x(tracks(j).frames<=frames(end)),tracks(j).y(tracks(j).frames<=frames(end)),'-')
end

%% Tracé trajectoires unstuck

%close all
lastPicture_temp = imread(fullfile(pictureDirectory, PicturesList{end}));
if calib.orientation=='r'
    lastPicture=imrotate(lastPicture_temp,-90);
elseif calib.orientation=='d'
    lastPicture=lastPicture_temp;
elseif calib.orientation=='u'
    lastPicture=imrotate(lastPicture_temp,-180);
elseif calib.orientation=='l'
    lastPicture=imrotate(lastPicture_temp,90);
end
% figure
% clf
subplot(2,2,2);
imshow(lastPicture*4)
axis on
hold on
title('Tracks unstuck')

for j=1:numel(tracks_unstuck)
    plot(tracks_unstuck(j).x(tracks_unstuck(j).frames<=frames(end)),tracks_unstuck(j).y(tracks_unstuck(j).frames<=frames(end)),'-')
end

%% Tracé trajectoires unstuck filtrees

%close all
lastPicture_temp = imread(fullfile(pictureDirectory, PicturesList{end}));
if calib.orientation=='r'
    lastPicture=imrotate(lastPicture_temp,-90);
elseif calib.orientation=='d'
    lastPicture=lastPicture_temp;
elseif calib.orientation=='u'
    lastPicture=imrotate(lastPicture_temp,-180);
elseif calib.orientation=='l'
    lastPicture=imrotate(lastPicture_temp,90);
end
% figure
% clf
subplot(2,2,3);
imshow(lastPicture*4)
axis on
hold on
title('Tracks unstuck filtrees')

for j=1:numel(tracks_unstuck)
    plot(tracks_unstuck_filtered(j).x(tracks_unstuck_filtered(j).frames<=frames(end)),tracks_unstuck_filtered(j).y(tracks_unstuck_filtered(j).frames<=frames(end)),'-')
    %(tracks(j).X(tracks(j).nFrame<),tracks(j).Y(tracks(j).nFrame==i),'x')
end

%% Tracé trajectoires unstuck long begin

%close all
lastPicture_temp = imread(fullfile(pictureDirectory, PicturesList{end}));
if calib.orientation=='r'
    lastPicture=imrotate(lastPicture_temp,-90);
elseif calib.orientation=='d'
    lastPicture=lastPicture_temp;
elseif calib.orientation=='u'
    lastPicture=imrotate(lastPicture_temp,-180);
elseif calib.orientation=='l'
    lastPicture=imrotate(lastPicture_temp,90);
end
% figure
% clf
subplot(2,2,4);
imshow(lastPicture*4)
axis on
hold on
title('Tracks unstuck long begin')

for j=1:numel(tracks_unstuck_long_begin)
    plot(tracks_unstuck_long_begin(j).x(tracks_unstuck_long_begin(j).frames<=frames(end)),tracks_unstuck_long_begin(j).y(tracks_unstuck_long_begin(j).frames<=frames(end)),'-')
end
%% Séquence images avec trajectoires
if sequence
    %close all
    tracks_chosen=tracks_unstuck_filtered;
    for k=10:500:round(length(PicturesList))
        Picture_temp = imread(fullfile(pictureDirectory, PicturesList{k}));
        if calib.orientation=='r'
            Picture=imrotate(Picture_temp,-90);
        elseif calib.orientation=='d'
            Picture=Picture_temp;
        elseif calib.orientation=='u'
            Picture=imrotate(Picture_temp,-180);
        elseif calib.orientation=='l'
            Picture=imrotate(Picture_temp,90);
        end

        figure(k)

        rect=calib.localROI{2};

        imshow(Picture*4);
        hold on
        title(['Tracks unstuck filtered- t=' num2str(k/calib.fps) 's' ])

        for j=1:numel(tracks_chosen)
            plot(tracks_chosen(j).x(tracks_chosen(j).frames<=frames(k)),tracks_chosen(j).y(tracks_chosen(j).frames<=frames(k)),'-m','LineWidth',2)
        end
        SavePictureName=fullfile(tracksDirectory, ['picture_traj_t=' num2str(k/calib.fps) 's.png']);
        saveas(gcf,SavePictureName)


    end
end
