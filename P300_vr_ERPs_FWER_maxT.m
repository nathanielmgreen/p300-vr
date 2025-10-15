clc; clear all; close all;
%% =============================
%            SETUP
% ==============================
mat_file    = 'epochs_data_08302025_new.mat';
chan_labels = {'Cz','Pz','Oz','O1','O2'};
baseline_ms = [-300 0];
xlim_ms     = [-300 800];

stats_win_inner = [0 400];    % Eye shifting
stats_win_outer = [200 600];  % Head turning

alpha_raw   = 0.05;           % gray dots/bands
alpha_fwer  = 0.05;           % familywise rate for permutation max-T (black)
stats_mode  = 'trial_level';  % epochs has no file IDs
n_perm      = 5000;
tail        = 'both';         % 'both'|'right'|'left'

% visuals
pillar_width_ms   = [];    % [] -> auto (≈ 4*dt)
alpha_gray_pillar = 0.18;
alpha_blk_pillar  = 0.22;
color_gray        = [0.7 0.7 0.7];
color_black       = [0 0 0];

%% =============================
%        LOAD & BUILD ERPs
% ==============================
S = load(mat_file); assert(isfield(S,'epochs'),'Missing epochs');
epochs = S.epochs; times = epochs.times(:)';

chIdxMap  = map_channel_indices(epochs.chanlocs, chan_labels);
ERP_INNER = build_erps_from_epochs(epochs, 'inner', times, chIdxMap, chan_labels, baseline_ms);
ERP_OUTER = build_erps_from_epochs(epochs, 'outer', times, chIdxMap, chan_labels, baseline_ms);

%% =============================
%         PLOTTING
% ==============================
plot_perm_maxT_panel(ERP_INNER, 'Eye Shifting (Inner)', chan_labels, xlim_ms, stats_win_inner, ...
                     alpha_raw, alpha_fwer, stats_mode, n_perm, tail, ...
                     pillar_width_ms, alpha_gray_pillar, alpha_blk_pillar, color_gray, color_black);

plot_perm_maxT_panel(ERP_OUTER, 'Head Turning (Outer)', chan_labels, xlim_ms, stats_win_outer, ...
                     alpha_raw, alpha_fwer, stats_mode, n_perm, tail, ...
                     pillar_width_ms, alpha_gray_pillar, alpha_blk_pillar, color_gray, color_black);


%% ============ Core panels ============
function plot_perm_maxT_panel(ERP, panel_title, chan_labels, xlim_ms, stats_win, ...
                              alpha_raw, alpha_fwer, mode, n_perm, tail, ...
                              pillar_width_ms, aGray, aBlk, cGray, cBlk)

locks = {'gip','fix','stim'}; row_titles = {'GIP-locked','Fixation-locked','Stimulus-locked'};
ylims = []; % set if you want fixed y

figure('Name',['ERPs + Permutation FWER (max-T) — ' panel_title], ...
       'Units','normalized','Position',[0.05 0.1 0.9 0.8],'Color','w');

for r=1:numel(locks)
    lk = locks{r}; t = ERP.(lk).times; T=numel(t);
    inwin = (t>=stats_win(1) & t<=stats_win(2));

    for c=1:numel(chan_labels)
        lab = chan_labels{c};
        Xs = ERP.trials.(lk).std.(lab);  % trials x time
        Xd = ERP.trials.(lk).dev.(lab);
        ax = subplot(3,numel(chan_labels),(r-1)*numel(chan_labels)+c); hold(ax,'on'); grid(ax,'on');

        if isempty(Xs) || isempty(Xd)
            if r==1, title(ax, lab); end
            if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
            xline(ax,0,'k:'); xlim(ax,xlim_ms); continue;
        end

        % Means/SEMs & lines
        [mS,sS,mD,sD] = means_sems(Xs,Xd);
        % plot(ax,t,mS,'b','LineWidth',1.3); fill(ax,[t fliplr(t)],[mS-sS fliplr(mS+sS)],'b','FaceAlpha',0.15,'EdgeColor','none');
        % plot(ax,t,mD,'r','LineWidth',1.3); fill(ax,[t fliplr(t)],[mD-sD fliplr(mD+sD)],'r','FaceAlpha',0.15,'EdgeColor','none');
        % ----- Shaded SEM first (kept out of legend), then lines -----
        if ~isempty(sS)
            fill(ax, [t fliplr(t)], [mS-sS fliplr(mS+sS)], ...
                'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        if ~isempty(sD)
            fill(ax, [t fliplr(t)], [mD-sD fliplr(mD+sD)], ...
                'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        
        hS = plot(ax, t, mS, 'b-', 'LineWidth', 1.3);  % Standard (blue line)
        hD = plot(ax, t, mD, 'r-', 'LineWidth', 1.3);  % Deviant  (red  line)
        
        % Only put the top row’s legend; reference the *line* handles explicitly
        if r==1
            legend(ax, [hS hD], {'Standard','Deviant'}, 'Location','northeast', 'Box','off');
        end

        % Raw parametric p (gray)
        p_raw = pointwise_pvals(Xs,Xd,mode);
        raw_mask = false(1,T); raw_mask(inwin) = p_raw(inwin) < alpha_raw;

        % Permutation FWER via max-T
        [fwer_mask, ~] = perm_maxT_mask(Xs,Xd,t,stats_win,n_perm,tail,alpha_fwer);

        % Lock axes, compute baseline
        xline(ax,0,'k:','HandleVisibility','off'); xlim(ax,xlim_ms);
        if ~isempty(ylims), ylim(ax,ylims); end
        ax.XLimMode='manual'; ax.YLimMode='manual'; yb=ylim(ax); ybase=yb(1);

        % Pillars
        width_ms = pillar_width(ax,t,pillar_width_ms);
        draw_vertical_bands(ax, t(raw_mask), yb, cGray, aGray, width_ms);
        draw_vertical_bands(ax, t(fwer_mask), yb, cBlk,  aBlk,  width_ms);

        % Dots at baseline
        scatter(ax, t(raw_mask),  repmat(ybase,1,sum(raw_mask)),  8,  cGray,'filled','HandleVisibility','off');
        scatter(ax, t(fwer_mask), repmat(ybase,1,sum(fwer_mask)), 10, 'k',  'filled','HandleVisibility','off');

        if r==1, title(ax, lab); end
        if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
        if r==numel(locks), xlabel(ax,'Time (ms)'); end
        if r==1, legend(ax,{'Standard','Deviant'},'Location','northeast','Box','off'); end
        set(ax,'Layer','top');
    end
end
sgtitle(sprintf('%s  (raw p<%.2f, Perm FWER α=%.2f in %d–%d ms)', panel_title, alpha_raw, alpha_fwer, stats_win(1), stats_win(2)), 'FontWeight','bold');
end

%% ============ Stats helpers ============
function [mask,thr] = perm_maxT_mask(Xa,Xb,tvec,win_ms,n_perm,tail,alpha)
% Label-shuffle to get null distribution of max-|T| within window; threshold at (1-alpha)
T = numel(tvec); idx = find(tvec>=win_ms(1) & tvec<=win_ms(2));
na = size(Xa,1); pool=[Xa;Xb]; N=size(pool,1);

% Observed t within window
t_obs = zeros(1,numel(idx));
for i=1:numel(idx), col=idx(i); t_obs(i)=welch_t(Xa(:,col),Xb(:,col)); end

% Permute → max statistic
maxT = zeros(n_perm,1);
for k=1:n_perm
    rp = randperm(N);
    Ap = pool(rp(1:na), idx); Bp = pool(rp(na+1:end), idx);
    tp = zeros(1,numel(idx));
    for j=1:numel(idx), tp(j)=welch_t(Ap(:,j),Bp(:,j)); end
    switch lower(tail)
        case 'both',  maxT(k) = max(abs(tp));
        case 'right', maxT(k) = max(tp);
        case 'left',  maxT(k) = max(-tp);
    end
end
thr = quantile(maxT, 1-alpha);

% Significant mask within window
switch lower(tail)
    case 'both',  win_mask = abs(t_obs) >= thr;
    case 'right', win_mask = t_obs >= thr;
    case 'left',  win_mask = -t_obs >= thr;
end

mask = false(1,T); mask(idx(win_mask)) = true;
end

function p = pointwise_pvals(Xa,Xb,mode)
T=size(Xa,2); p=nan(1,T);
switch lower(mode)
    case 'trial_level'
        for tt=1:T
            xa=Xa(:,tt); xb=Xb(:,tt);
            xa=xa(~isnan(xa)); xb=xb(~isnan(xb));
            if numel(xa)>=2 && numel(xb)>=2
                [~,p(tt)] = ttest2(xa,xb,'Vartype','unequal');
            end
        end
    otherwise, error('Only trial_level supported.');
end
end

%% ============ ERP building & utilities ============
function ERP = build_erps_from_epochs(epochs, side, times, chIdxMap, chan_labels, baseline_ms)
locks_src={'gip_onset','fix_onset','stimulus_onset'}; locks_out={'gip','fix','stim'}; ERP=struct();
for L=1:numel(locks_src)
    src=locks_src{L}; out=locks_out{L};
    data = epochs.(src).([side '_data']); labels = epochs.(src).([side '_stimuli']);
    if ischar(labels)||isstring(labels), labels=cellstr(labels); end
    if iscategorical(labels), labels=cellstr(labels); end
    labels=labels(:).'; is_std=strcmpi(labels,'standard'); is_dev=strcmpi(labels,'deviant');

    % baseline
    if ~isempty(baseline_ms), data = baseline_epochs(data,times,baseline_ms); end

    ERP.(out).times = times;
    for c=1:numel(chan_labels)
        lab=chan_labels{c}; idx = chIdxMap.(lab);
        if isempty(idx), ERP.trials.(out).std.(lab)=[]; ERP.trials.(out).dev.(lab)=[]; continue; end
        X = squeeze(data(idx,:,:))'; % trials x time
        ERP.trials.(out).std.(lab)=X(is_std,:); ERP.trials.(out).dev.(lab)=X(is_dev,:);
    end
end
end

function chIdxMap = map_channel_indices(chanlocs, wanted_labels)
chIdxMap=struct(); have=arrayfun(@(x) string(x.labels),chanlocs,'UniformOutput',true);
for k=1:numel(wanted_labels)
    lab=wanted_labels{k}; idx=find(strcmpi(have,lab),1);
    if isempty(idx), warning('Channel %s not found.',lab); chIdxMap.(lab)=[]; else, chIdxMap.(lab)=idx; end
end
end

function data_bc = baseline_epochs(data,times,baseline_ms)
data_bc=data; t=times(:)'; bmin=max(baseline_ms(1),t(1)); bmax=min(baseline_ms(2),t(end));
if ~(bmin<bmax), warning('baseline: invalid window'); return; end
bmask=(t>=bmin & t<=bmax); if ~any(bmask), warning('baseline: no samples'); return; end
for tr=1:size(data,3)
    bl=mean(data(:,bmask,tr),2,'omitnan');
    data_bc(:,:,tr)=data(:,:,tr)-bl;
end
end

function [mS,sS,mD,sD]=means_sems(Xs,Xd)
mS=mean(Xs,1,'omitnan'); sS=std(Xs,0,1,'omitnan')/sqrt(max(1,size(Xs,1)));
mD=mean(Xd,1,'omitnan'); sD=std(Xd,0,1,'omitnan')/sqrt(max(1,size(Xd,1)));
end

function t=welch_t(x,y)
x=x(~isnan(x)); y=y(~isnan(y)); nx=numel(x); ny=numel(y);
if nx<2||ny<2, t=0; return; end
vx=var(x,0); vy=var(y,0); mx=mean(x); my=mean(y);
se=sqrt(vx/nx+vy/ny); t = (se==0) * 0 + (se~=0) * ((mx-my)/se);
end

function draw_vertical_bands(ax, x_positions, ylims, rgb, alpha, width_ms)
if isempty(x_positions), return; end
hb=gobjects(0);
for k=1:numel(x_positions)
    x0=x_positions(k)-width_ms/2; x1=x_positions(k)+width_ms/2;
    hb(end+1)=patch(ax,[x0 x1 x1 x0],[ylims(1) ylims(1) ylims(2) ylims(2)],rgb, ...
                    'FaceAlpha',alpha,'EdgeColor','none','HandleVisibility','off');
end
try uistack(hb,'bottom'); catch, end
end

function w = pillar_width(ax,t,pw)
if ~isempty(pw), w=pw; else
    if numel(t)>1, dt=median(diff(t)); else, dt=1; end
    w=4*dt;
end
end