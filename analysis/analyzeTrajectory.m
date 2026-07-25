function analysis = analyzeTrajectory(tracks, calib, options)
pixel_size  = 0.1625;     % µm/pixel pour le 40x
nTracks = length(tracks); %nombre de traj
%% Options par défaut
if nargin < 3
    options.smoothWindow = 5;
    options.nbins = 100;
end
%% PDF des déplacements selon Y
window_size = 101;    % frames par fenêtre
deltaT = 10;          % intervalle de temps pour le déplacement (en frames)
nbins = 200;          % nombre de bins pour l'histogramme
frame_start = 1;      % début de la fenêtre absolue
frame_end=frame_start+window_size-1;

all_depY = [];

for k=1:length(tracks)
    Y=tracks(k).y;
    % on prend seulement les frames qui existent dans cette traj
    s=max(1, frame_start); % pour être sûr de ne pas dépasser le début
    e=min(frame_end, length(Y)); % ne pas dépasser la longueur de la trajectoire
    Yw=Y(s:e);
    depY_dt= -(Yw(1+deltaT:end) - Yw(1:end-deltaT))*pixel_size;
    all_depY = [all_depY; depY_dt];
end

% histogramme normalisé
[nY, yedges] = hist(all_depY, nbins);
% PDF_Y =nY/numel(all_depY);
PDF_Y = nY / ((numel(all_depY) * (yedges(2)-yedges(1))));


figure;
plot(yedges, PDF_Y, 'LineWidth',1.5);
xlabel('Déplacement Y (\mum)','FontSize', 24); ylabel('PDF','FontSize', 24);
set(gca, 'FontSize', 22)
title(sprintf('Distribution des déplacements en Y - frames %d à %d', frame_start, frame_end),'FontSize', 30, 'FontWeight', 'bold');
grid on;

%% Tortuosité classique

intervals = [2 400; 401 1000; 2000 5000];  % intervalles absolus
num_tracks = numel(tracks);
nIntervals = size(intervals,1);

% Start/end frame de chaque trajectoire 
frame_start = nan(num_tracks,1);
frame_end   = nan(num_tracks,1);

for i = 1:num_tracks
    frame_start(i) = tracks(i).frames(1);
    frame_end(i) = tracks(i).frames(1) + tracks(i).L- 1;
end
% Matrice tortuosité
Tortu_mat = nan(num_tracks, nIntervals);
count_tort_per_interval = zeros(nIntervals,1);
% Boucle principale
for it = 1:nIntervals
    Fstart = intervals(it,1);
    Fend   = intervals(it,2);
    for i = 1:num_tracks
        sF = frame_start(i);
        eF = frame_end(i);
        % Overlap
        ov_start = max(sF, Fstart);
        ov_end   = min(eF, Fend);
        % Meilleur précision
        if ov_end - ov_start + 1 < 5
            continue;
        end
        % Indices locaux dans la trajectoire
        idx1 = ov_start- sF + 1;
        idx2 = ov_end- sF + 1;
        
        Xw = tracks(i).x(idx1:idx2);
        Yw = tracks(i).y(idx1:idx2);
        
        if numel(Xw) < 2
            continue;
        end
        %Longueur totale
        Ltot = 0;
        for j= 1:(numel(Xw)-1)
            dx = Xw(j+1) - Xw(j);
            dy = Yw(j+1) - Yw(j);
            Ltot = Ltot + sqrt(dx^2 + dy^2);
        end        
        % Distance directe
        dmax = sqrt( (Xw(end)-Xw(1))^2 + (Yw(end)-Yw(1))^2 );
        if dmax == 0
            continue;
        end
        % Tortuosité
        Tortu_mat(i,it) = Ltot / dmax;
    end
    
    % Nombre de trajectoires présentes dans l'intervalle
    count_tort_per_interval(it) = sum(~isnan(Tortu_mat(:,it)));
end

% Tortuosité moyenne par intervalle
Tortu_mean = nanmean(Tortu_mat,1);

% Affichage
figure;
bar(1:nIntervals, Tortu_mean);
xticks(1:nIntervals);
xticklabels(arrayfun(@(r) sprintf('%d-%d', intervals(r,1), intervals(r,2)),1:nIntervals, 'UniformOutput', false))
xlabel('Intervalle (frames)','FontSize', 24);
ylabel('Tortuosité moyenne L / dmax','FontSize', 24);
set(gca, 'FontSize', 22)
title('Tortuosité moyenne par intervalle (temps absolu)','FontSize', 34, 'FontWeight', 'bold');
grid on;

for it = 1:nIntervals
    text(it, Tortu_mean(it), sprintf('n=%d', count_tort_per_interval(it)),'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontSize', 20, 'FontWeight','bold');
end
%% MSD (Mean Squared Displacement)
min_points_for_alpha = 3;    % nombre minimal de points pour la regression
min_frames_in_interval = 50; % seuil pour garder un segment 

tau1_frac = 0.325;  % (65% to 90% des grands tau)
tau2_frac = 0.45;  

%filtrage par nombre minimal de paires par tau
min_pairs_per_tau = 5; 

intervals = [2 400; 401 1000; 2000 4000; 4000 7000];
num_tracks = numel(tracks);

% Start/end frame de chaque trajectoire dans le référentiel absolu
frame_start= nan(num_tracks,1);
frame_end  = nan(num_tracks,1);
len_tracks = zeros(num_tracks,1);

% Définir nIntervals avant d'initialiser les cell
nIntervals = size(intervals,1);
MSD_store = cell(nIntervals,1);  
TAU_store = cell(nIntervals,1);  
for i = 1:num_tracks
    len_tracks(i) = tracks(i).L;
    frame_start(i) = tracks(i).frames(1);
    frame_end(i) = tracks(i).frames(1) + tracks(i).L-1;
end

% Calcul alpha par trajectoire pour chaque intervalle
alpha_mat = nan(num_tracks, nIntervals);
count_per_interval = zeros(nIntervals,1);

for it = 1:nIntervals
    Fstart = intervals(it,1); %min et max de l'interval
    Fend   = intervals(it,2);
    for i = 1:num_tracks
        sF = frame_start(i);
        eF = frame_end(i);
        % overlap entre la trajectoire et l'intervalle
        ov_start = max(sF, Fstart);
        ov_end   = min(eF, Fend);
        % minimum frames dans l'intervale
        if ov_end - ov_start + 1 < min_frames_in_interval
            continue
        end
        % indices locaux
        local_idx_start = ov_start - sF + 1;
        local_idx_end   = ov_end   - sF + 1;

        Xi = tracks(i).x(local_idx_start:local_idx_end) * pixel_size;
        Yi = tracks(i).y(local_idx_start:local_idx_end) * pixel_size;

        Lseg = numel(Xi);
        tau_max_seg = floor(Lseg/2);
        if tau_max_seg < 2, continue; end

        % MSD du segment (et nombre de paires par tau)
        msd_seg = nan(1, tau_max_seg);
        Npairs = zeros(1, tau_max_seg);
        for k = 1:tau_max_seg
            id1 = 1:(Lseg - k);
            id2 = 1+k : Lseg;
            d2 = (Xi(id1)-Xi(id2)).^2 + (Yi(id1)-Yi(id2)).^2;
            msd_seg(k) = mean(d2);
            Npairs(k) = numel(d2);
        end

        % Stockage MSD (alignement progressif)
        if isempty(MSD_store{it})
            MSD_store{it} = msd_seg(:)';   % première trajectoire
        else
            Lprev = size(MSD_store{it},2);
            Lnew  = numel(msd_seg);

            % aligner les tailles
            if Lnew > Lprev
                MSD_store{it}(:,end+1:Lnew) = NaN;
            elseif Lnew < Lprev
                msd_seg(end+1:Lprev) = NaN;
            end
            MSD_store{it} = [MSD_store{it}; msd_seg];
        end

        TAU_store{it} = (1:numel(msd_seg)) / calib.fps;

        % calculer bornes en frames (tau)
        tau1 = ceil(tau1_frac * Lseg);
        tau2 = floor(tau2_frac * Lseg);

        % sécurité
        tau1 = max(tau1, 3);               % éviter tau = 1,2
        tau2 = min(tau2, tau_max_seg);    % ne pas dépasser tau_max_seg
        if tau2 <= tau1, continue; end

        idxs = tau1:tau2;

        % filtrer par nombre de paires si demandé
        if min_pairs_per_tau > 0
            idxs = idxs(Npairs(idxs) >= min_pairs_per_tau);
        end

        if numel(idxs) < min_points_for_alpha, continue; end

        tps = idxs / calib.fps;
        msd_forfit = msd_seg(idxs);

        valid_fit = msd_forfit > 0 & isfinite(msd_forfit);
        if sum(valid_fit) < min_points_for_alpha, continue; end

        p = polyfit(log(tps(valid_fit)), log(msd_forfit(valid_fit)), 1);
        alpha_mat(i,it) = p(1);
    end
    count_per_interval(it) = sum(~isnan(alpha_mat(:,it)));
end

% alpha moyen pour chaque intervalle
alpha_mean = nanmean(alpha_mat,1);
figure;
for it = 1:nIntervals
    subplot(2, ceil(nIntervals/2), it);

    M = MSD_store{it};
    if isempty(M)
        title(sprintf('%d-%d : no data', intervals(it,1), intervals(it,2)));
        continue
    end
    tau = TAU_store{it};
    % MSD moyen
    msd_mean = nanmean(M,1);

% sécurité
valid = msd_mean > 0 & isfinite(msd_mean);

loglog(tau(valid), msd_mean(valid), 'o-', 'LineWidth', 1.8);
grid on;

xlabel('\tau (s)','FontSize', 24);
ylabel('MSD (\mum^2)','FontSize', 24);
title(sprintf('%d–%d frames | \\alpha = %.2f', intervals(it,1), intervals(it,2),alpha_mean(it)), 'FontSize', 20, 'FontWeight', 'bold');

end
sgtitle('MSD moyen (log–log) par intervalle de temps (grands τ, 25% des τ, extrêmes exclus)','FontSize', 34, 'FontWeight', 'bold');

figure;
b = bar(1:nIntervals, alpha_mean);
xticks(1:nIntervals);
xticklabels(arrayfun(@(r) sprintf('%d-%d', intervals(r,1), intervals(r,2)), 1:nIntervals, 'UniformOutput', false));
xlabel('Intervalle (frames)','FontSize', 24);
ylabel('\alpha moyen','FontSize', 24);
set(gca, 'FontSize', 22)
title('Alpha moyen par intervalle (grands τ, 25% des τ, valeurs extrêmes exclues)','FontSize', 34, 'FontWeight', 'bold');

grid on;

% Ajouter le nombre de trajectoires centré au-dessus de chaque barre
for it = 1:nIntervals
    text(it, alpha_mean(it), sprintf('n=%d', count_per_interval(it)), ...
        'HorizontalAlignment','center', ...   % centre le texte sur la barre
        'VerticalAlignment','bottom', ...
        'FontSize', 20, 'FontWeight','bold');
end


figure;
boxplot(alpha_mat, 'Labels', arrayfun(@(r) sprintf('%d-%d', intervals(r,1), intervals(r,2)), 1:nIntervals, 'UniformOutput', false));
ylabel('\alpha','FontSize', 24);
xlabel('Intervalle (frames)');
set(gca, 'FontSize', 22)
title('Distribution de α par intervalle de frames','FontSize', 34, 'FontWeight', 'bold');

grid on;

figure; hold on;
for it = 1:nIntervals
    xj = it + 0.1*(rand(1,count_per_interval(it))-0.5); % jitter
    plot(xj, alpha_mat(~isnan(alpha_mat(:,it)), it), 'o', 'MarkerFaceColor','b','MarkerEdgeColor','k');
end
plot(1:nIntervals, alpha_mean, 'r-', 'LineWidth',2); % moyenne
xticks(1:nIntervals);
xticklabels(arrayfun(@(r) sprintf('%d-%d', intervals(r,1), intervals(r,2)), 1:nIntervals, 'UniformOutput', false));
ylabel('\alpha','FontSize', 24);
xlabel('Intervalle (frames)','FontSize', 24);
set(gca, 'FontSize', 22)
title('Alpha par trajectoire avec moyenne par intervalle','FontSize', 34, 'FontWeight', 'bold');

grid on;

% nombre de trajectoires décalé légèrement à droite
dx = 0.2;  % décalage sur x
for it = 1:nIntervals
    text(it + dx, alpha_mean(it), sprintf('n=%d', count_per_interval(it)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 20, 'FontWeight','bold');
end

%% Résultat
analysis = struct();
analysis.PDF_Y = PDF_Y;
analysis.yedges = yedges;
OUT = struct();
OUT.alpha_mat = alpha_mat;
OUT.alpha_mean = alpha_mean;
OUT.count_per_interval = count_per_interval;
OUT.intervals = intervals;
OUT.min_points_for_alpha = min_points_for_alpha;

end