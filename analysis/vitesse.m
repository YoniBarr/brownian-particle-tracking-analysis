function [OUT, tracks] = vitesse(tracks, calib)

% Paramètres
px2um = 0.1625; %au 40x
dt = 1 / calib.fps;

num_tracks = numel(tracks);
Vtot_time = cell(num_tracks,1);
Vmean_perTrack = nan(num_tracks,1);
frame_start = nan(num_tracks,1);
frame_end   = nan(num_tracks,1);
Vtot_all = [];

% Boucle principale
for i = 1:num_tracks

    % longueur
    if isfield(tracks(i),'L') && ~isempty(tracks(i).L)
        len = tracks(i).L;
    else
        len = numel(tracks(i).x);
    end

    if len <= 1
        frame_start(i) = 1;
        frame_end(i)   = 0;
        Vtot_time{i} = [];
        % ensure fields exist
        tracks(i).vx = NaN(len,1);
        tracks(i).vy = NaN(len,1);
        tracks(i).vtot= Nan(len,1);
        continue
    end

    % frame de début si existant
    if isfield(tracks(i),'frames') && ~isempty(tracks(i).frames)
        startF = tracks(i).frames(1);
        frames_i = tracks(i).frames(:)';
    else
        startF = 1;
        frames_i = (1:len);
        tracks(i).frames = frames_i(:);
    end

    len = numel(frames_i);
    frame_start(i) = startF;
    frame_end(i) = startF + len - 1;

    % Positions (assure ligne)
    X = tracks(i).x(:)';
    Y = tracks(i).y(:)';

    % Vitesses (Vx,Vy length = len-1)
    Vx = diff(X) / dt * px2um;
    Vy = diff(Y) / dt * px2um;
    Vtot = sqrt(Vx.^2 + Vy.^2);

    % aligner à la longueur des positions : prepend NaN
    vx_full = [NaN, Vx];  % 1 x len
    vy_full = [NaN, Vy];
    v_total_full = [NaN, Vtot];

    % injecter dans tracks (col vectors)
    tracks(i).vx = vx_full(:);
    tracks(i).vy = vy_full(:);
    tracks(i).vtot= v_total_full(:);
    % stockages pour OUT
    OUT.Vx_time{i} = vx_full(:);
    OUT.Vy_time{i} = vy_full(:);
    OUT.Vtot_time{i} = Vtot(:);
    Vtot_time{i} = Vtot(:);
    frame_end(i) = startF + numel(Vtot) - 1;
    Vmean_perTrack(i) = mean(Vtot, 'omitnan');

    % accumulation pour histogramme
    Vtot_all = [Vtot_all; Vtot(:)];
end

% Construction du temps absolu global
global_minF = min(frame_start);
global_maxF = max(frame_end);
frame_vec = global_minF : global_maxF;
nFrames = numel(frame_vec);

% Matrice alignée [nTracks x nFrames]
V_matrix = nan(num_tracks, nFrames);
for i = 1:num_tracks
    V = Vtot_time{i};
    if isempty(V), continue; end
    f0 = frame_start(i);
    frames_idx = f0 : (f0 + numel(V) - 1);
    cols = frames_idx - global_minF + 1;
    valid = cols >= 1 & cols <= nFrames;
    V_matrix(i, cols(valid)) = V(valid);
end

% Moyenne temporelle
Vmean_time = mean(V_matrix, 1, 'omitnan');
count_time = sum(~isnan(V_matrix), 1);
time_vec = (frame_vec - global_minF) / calib.fps;
 
%% PDF simplifiée de la vitesse totale (remplace la section précédente)
window_size = 601;   % frames par fenêtre (info visuelle)
deltaT = 10;          % intervalle en frames pour estimer la vitesse
nbins = 200;
pdf_start = 400;
pdf_end = pdf_start + window_size - 1;

all_vtot = [];

nTracks = numel(tracks);
for k = 1:nTracks
    X = tracks(k).x(:);
    Y = tracks(k).y(:);
    if numel(X) <= deltaT
        continue;
    end

    % Extraire la trajectoire correspondant à [pdf_start, pdf_end]
    if isfield(tracks(k),'frames') && ~isempty(tracks(k).frames)
        frames_k = tracks(k).frames(:)';
        % indices des frames de la trajectoire qui tombent dans la fenêtre globale
        idx = find(frames_k >= pdf_start & frames_k <= pdf_end);
        if numel(idx) <= deltaT
            continue;
        end
        Xw = X(idx);
        Yw = Y(idx);
    else
        % pas de décalage de frames : on prend l'intervalle d'indice
        s = max(1, pdf_start);
        e = min(pdf_end, numel(X));
        if (e - s + 1) <= deltaT
            continue;
        end
        Xw = X(s:e);
        Yw = Y(s:e);
    end

    % calcul vectorisé des déplacements sur deltaT frames
    dX = Xw((1+deltaT):end) - Xw(1:(end-deltaT));
    dY = Yw((1+deltaT):end) - Yw(1:(end-deltaT));

    % vitesse totale (µm/s)
    vtot = sqrt(dX.^2 + dY.^2) * px2um / (deltaT * dt);

    all_vtot = [all_vtot; vtot(:)];  % concaténation 
end


% histogramme + normalisation manuelle (PDF)
    [counts, edges] = histcounts(all_vtot, nbins);
    binwidth = edges(2) - edges(1);            % suppose bins uniformes
    centers = edges(1:end-1) + binwidth/2;
    PDF_vtot = counts / (numel(all_vtot) * binwidth);

    % plot
    figure;
    plot(centers, PDF_vtot, 'LineWidth', 1.5);
    xlabel('Vitesse totale (µm/s)','FontSize', 24);
    ylabel('PDF','FontSize', 24);
    title(sprintf('PDF de la vitesse totale (\\Delta t = %d frames)', deltaT),'FontSize', 30, 'FontWeight', 'bold');
    grid on;




% Courbe vitesse moyenne vs temps absolu 
figure; hold on;
idx = count_time > 0;
t = time_vec(idx);
plot(t, Vmean_time(idx), 'b-', 'LineWidth', 2);
xlabel('Temps absolu (s)','FontSize', 24);
ylabel('Vitesse moyenne (µm/s)','FontSize', 24);
set(gca, 'FontSize', 22)
title('Vitesse moyenne au cours du temps','FontSize',30, 'FontWeight', 'bold')
grid on;
hold off;


% Sortie OUT
OUT.V_matrix  = V_matrix;
OUT.Vmean_time = Vmean_time;
OUT.count_time = count_time;
OUT.Vtot_time = Vtot_time;
OUT.Vmean_perTrack = Vmean_perTrack;
OUT.Vtot_all = Vtot_all;
OUT.frame_vec = frame_vec;
OUT.tracks_with_v = tracks;

end

