clc; close all;

%% =============================
%            SETUP
% ==============================
mat_file    = 'epochs_data_08302025_new.mat';
chan_labels = {'Cz','Pz','Oz','O1','O2'};
baseline_ms = [-300 0];
xlim_ms     = [-300 800];

stats_win_inner = [0 400];
stats_win_outer = [200 600];

alpha_raw = 0.05;   % gray
q_fdr     = 0.05;   % BH-FDR (black)
tail      = 'both'; % used only for sign of t if you want; p is from ttest2 (two-tailed)
stats_mode= 'trial_level';

pillar_width_ms   = [];    alpha_gray_pillar=0.18; alpha_blk_pillar=0.22;
color_gray=[0.7 0.7 0.7]; color_black=[0 0 0];

%% =============================
%        LOAD & BUILD ERPs
% ==============================
S = load(mat_file); assert(isfield(S,'epochs'),'Missing epochs');
epochs = S.epochs; times = epochs.times(:)';

chIdxMap  = map_channel_indices(epochs.chanlocs, chan_labels);
ERP_INNER = build_erps_from_epochs(epochs, 'inner', times, chIdxMap, chan_labels, baseline_ms);
ERP_OUTER = build_erps_from_epochs(epochs, 'outer', times, chIdxMap, chan_labels, baseline_ms);

%% =============================
%          PLOTTING
% ==============================
plot_FDR_panel(ERP_INNER,'Eye Shifting (Inner)',chan_labels,xlim_ms,stats_win_inner,alpha_raw,q_fdr,stats_mode, ...
               pillar_width_ms,alpha_gray_pillar,alpha_blk_pillar,color_gray,color_black);

plot_FDR_panel(ERP_OUTER,'Head Turning (Outer)',chan_labels,xlim_ms,stats_win_outer,alpha_raw,q_fdr,stats_mode, ...
               pillar_width_ms,alpha_gray_pillar,alpha_blk_pillar,color_gray,color_black);

function plot_FDR_panel(ERP,panel_title,chan_labels,xlim_ms,stats_win,alpha_raw,q_fdr,mode, ...
                        pillar_width_ms,aGray,aBlk,cGray,cBlk)
locks={'gip','fix','stim'}; row_titles={'GIP-locked','Fixation-locked','Stimulus-locked'};
ylims=[];

figure('Name',['ERPs + FDR — ' panel_title],'Units','normalized','Position',[0.05 0.1 0.9 0.8],'Color','w');

for r=1:numel(locks)
    lk=locks{r}; t=ERP.(lk).times; T=numel(t);
    inwin=(t>=stats_win(1)&t<=stats_win(2));
    for c=1:numel(chan_labels)
        lab=chan_labels{c}; ax=subplot(3,numel(chan_labels),(r-1)*numel(chan_labels)+c); hold(ax,'on'); grid(ax,'on');

        Xs=ERP.trials.(lk).std.(lab); Xd=ERP.trials.(lk).dev.(lab);
        if isempty(Xs)||isempty(Xd)
            if r==1, title(ax,lab); end
            if c==1, ylabel(ax,sprintf('%s\n\\muV',row_titles{r})); end
            xline(ax,0,'k:'); xlim(ax,xlim_ms); continue;
        end

        [mS,sS,mD,sD]=means_sems(Xs,Xd);
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

        % pointwise parametric p
        p = pointwise_pvals(Xs,Xd,mode);
        raw_mask=false(1,T); raw_mask(inwin)=p(inwin) < alpha_raw;

        % FDR (BH) only within the window
        fdr_mask=false(1,T);
        p_sub=p(inwin); idx=find(inwin); good=~isnan(p_sub);
        if any(good)
            [~,~,adj_p]=fdr_bh_local(p_sub(good),q_fdr);
            sig = adj_p < q_fdr;
            fdr_mask(idx(good(sig))) = true;
        end

        xline(ax,0,'k:'); xlim(ax,xlim_ms); if ~isempty(ylims), ylim(ax,ylims); end
        ax.XLimMode='manual'; ax.YLimMode='manual'; yb=ylim(ax); ybase=yb(1);

        width_ms=pillar_width(ax,t,pillar_width_ms);
        draw_vertical_bands(ax,t(raw_mask),yb,cGray,aGray,width_ms);
        draw_vertical_bands(ax,t(fdr_mask), yb,cBlk, aBlk, width_ms);

        scatter(ax,t(raw_mask),repmat(ybase,1,sum(raw_mask)),8,cGray,'filled','HandleVisibility','off');
        scatter(ax,t(fdr_mask), repmat(ybase,1,sum(fdr_mask)),10,'k','filled','HandleVisibility','off');

        if r==1, title(ax,lab); end
        if c==1, ylabel(ax,sprintf('%s\n\\muV',row_titles{r})); end
        if r==numel(locks), xlabel(ax,'Time (ms)'); end
        if r==1, legend(ax,{'Standard','Deviant'},'Location','northeast','Box','off'); end
        set(ax,'Layer','top');
    end
end
sgtitle(sprintf('%s  (raw p<%.2f, BH-FDR q<%.2f in %d–%d ms)', ...
    panel_title, alpha_raw, q_fdr, stats_win(1), stats_win(2)), 'FontWeight','bold');
end

%% ==== Minimal helpers (same as Script #1) ====
function chIdxMap = map_channel_indices(chanlocs, wanted_labels)
chIdxMap=struct(); have=arrayfun(@(x) string(x.labels),chanlocs,'UniformOutput',true);
for k=1:numel(wanted_labels)
    lab=wanted_labels{k}; idx=find(strcmpi(have,lab),1);
    if isempty(idx), warning('Channel %s not found.',lab); chIdxMap.(lab)=[]; else, chIdxMap.(lab)=idx; end
end
end

function ERP = build_erps_from_epochs(epochs, side, times, chIdxMap, chan_labels, baseline_ms)
locks_src={'gip_onset','fix_onset','stimulus_onset'}; locks_out={'gip','fix','stim'}; ERP=struct();
for L=1:numel(locks_src)
    src=locks_src{L}; out=locks_out{L};
    data=epochs.(src).([side '_data']); labels=epochs.(src).([side '_stimuli']);
    if ischar(labels)||isstring(labels), labels=cellstr(labels); end
    if iscategorical(labels), labels=cellstr(labels); end
    labels=labels(:).'; is_std=strcmpi(labels,'standard'); is_dev=strcmpi(labels,'deviant');
    if ~isempty(baseline_ms), data=baseline_epochs(data,times,baseline_ms); end
    ERP.(out).times=times;
    for c=1:numel(chan_labels)
        lab=chan_labels{c}; idx=chIdxMap.(lab);
        if isempty(idx), ERP.trials.(out).std.(lab)=[]; ERP.trials.(out).dev.(lab)=[]; continue; end
        X=squeeze(data(idx,:,:))';
        ERP.trials.(out).std.(lab)=X(is_std,:); ERP.trials.(out).dev.(lab)=X(is_dev,:);
    end
end
end

function data_bc = baseline_epochs(data,times,baseline_ms)
data_bc=data; t=times(:)'; bmin=max(baseline_ms(1),t(1)); bmax=min(baseline_ms(2),t(end));
if ~(bmin<bmax), warning('baseline: invalid'); return; end
bmask=(t>=bmin&t<=bmax); if ~any(bmask), warning('baseline: empty'); return; end
for tr=1:size(data,3), bl=mean(data(:,bmask,tr),2,'omitnan'); data_bc(:,:,tr)=data(:,:,tr)-bl; end
end

function [mS,sS,mD,sD]=means_sems(Xs,Xd)
mS=mean(Xs,1,'omitnan'); sS=std(Xs,0,1,'omitnan')/sqrt(max(1,size(Xs,1)));
mD=mean(Xd,1,'omitnan'); sD=std(Xd,0,1,'omitnan')/sqrt(max(1,size(Xd,1)));
end

function p = pointwise_pvals(Xa,Xb,mode)
T=size(Xa,2); p=nan(1,T);
switch lower(mode)
    case 'trial_level'
        for tt=1:T
            xa=Xa(:,tt); xb=Xb(:,tt); xa=xa(~isnan(xa)); xb=xb(~isnan(xb));
            if numel(xa)>=2 && numel(xb)>=2, [~,p(tt)]=ttest2(xa,xb,'Vartype','unequal'); end
        end
    otherwise, error('Only trial_level supported.');
end
end

function [h, crit_p, adj_p] = fdr_bh_local(pvals, q)
if nargin<2||isempty(q), q=0.05; end
p = pvals(:); n=numel(p);
[sp,idx]=sort(p,'ascend'); thresh=(1:n)'/n*q;
h = sp<=thresh;
crit_p=NaN; if any(h), crit_p=max(sp(h)); end
wtd_p = n*sp./(1:n)'; adj_temp=cummin(wtd_p(end:-1:1)); adj_temp=adj_temp(end:-1:1);
adj_p=zeros(n,1); adj_p(idx)=min(1,adj_temp);
end

function draw_vertical_bands(ax,x_positions,ylims,rgb,alpha,width_ms)
if isempty(x_positions), return; end
hb=gobjects(0);
for k=1:numel(x_positions)
    x0=x_positions(k)-width_ms/2; x1=x_positions(k)+width_ms/2;
    hb(end+1)=patch(ax,[x0 x1 x1 x0],[ylims(1) ylims(1) ylims(2) ylims(2)],rgb, ...
                    'FaceAlpha',alpha,'EdgeColor','none','HandleVisibility','off');
end
try uistack(hb,'bottom'); catch, end
end

function w=pillar_width(ax,t,pw)
if ~isempty(pw), w=pw; else, w=(numel(t)>1)*4*median(diff(t)) + (numel(t)<=1)*4; end
end
