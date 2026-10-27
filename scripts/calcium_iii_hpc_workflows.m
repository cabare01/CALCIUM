function calcium_iii_hpc_workflows()
% CALCIUM_III_HPC_WORKFLOWS
% MATLAB batch entry point for workflows described in Calcium III Wavelet Methods.
% Set CALCIUM_WORKFLOW and the input environment variables before calling.
% Each run requires a new CALCIUM_OUTPUT_ROOT; existing directories are refused.

workflow = getenv('CALCIUM_WORKFLOW');
if isempty(workflow)
    error('Set CALCIUM_WORKFLOW to spectral, correlation, corr_summary, individual_maps, combine_maps, group_spatial, or boundary_qc.');
end
outRoot = getenv('CALCIUM_OUTPUT_ROOT');
if isempty(outRoot)
    error('Set CALCIUM_OUTPUT_ROOT to a new, non-existing directory.');
end
if isfolder(outRoot) || isfile(outRoot)
    error('Refusing to use existing output path: %s', outRoot);
end
mkdir(outRoot);
fprintf('Workflow: %s\nOutput: %s\n', workflow, outRoot);

switch lower(workflow)
    case 'spectral'
        runSpectral(outRoot);
    case 'correlation'
        runCorrelations(outRoot);
    case 'corr_summary'
        runCorrelationSummary(outRoot);
    case 'individual_maps'
        runIndividualMaps(outRoot);
    case 'combine_maps'
        runCombineMaps(outRoot);
    case 'group_spatial'
        runGroupSpatial(outRoot);
    case 'boundary_qc'
        runBoundaryQC(outRoot);
    otherwise
        error('Unknown CALCIUM_WORKFLOW: %s', workflow);
end
fprintf('Workflow finished: %s\n', workflow);
end

function runSpectral(outRoot)
% PDF spectral driver: 0.02-0.20 Hz, 80 bins, analytic Morlet CWT.
rawRoot = requiredEnv('CALCIUM_RAW_TRACE_ROOT');
pattern = getenv('CALCIUM_RAW_TRACE_PATTERN');
if isempty(pattern), pattern = '*.csv'; end
files = dir(fullfile(rawRoot, pattern));
if isempty(files), error('No raw trace files found under %s (%s).', rawRoot, pattern); end
embryoNames=cell(numel(files),1);
for k=1:numel(files), [~,b,~]=fileparts(files(k).name); embryoNames{k}=safeName(erase(b,' - Not_self_Merged')); end
if numel(unique(embryoNames))~=numel(embryoNames), error('Input files resolve to duplicate embryo output names; refusing possible overwrite.'); end
freqAxis = linspace(0.02, 0.20, 80);
makePlots = strcmpi(getenv('CALCIUM_MAKE_TRACE_PLOTS'),'true');
for fi = 1:numel(files)
    inFile = fullfile(files(fi).folder, files(fi).name);
    [~, base, ~] = fileparts(files(fi).name);
    embryo = erase(base, ' - Not_self_Merged');
    embryoDir = fullfile(outRoot, embryoNames{fi});
    mkdir(embryoDir); mkdir(fullfile(embryoDir,'tables'));
    if makePlots, mkdir(fullfile(embryoDir,'figures')); end
    fprintf('Processing %d/%d: %s\n', fi, numel(files), files(fi).name);
    T = readtable(inFile, 'VariableNamingRule','preserve');
    if ~ismember('ID', T.Properties.VariableNames)
        error('Missing ID column in %s', inFile);
    end
    ids = unique(string(T.ID), 'stable');
    resultRows = cell(numel(ids),1);
    lib = struct('CellID',{},'RunID',{},'SheetName',{},'SourceFile',{}, ...
        'MeanX',{},'MeanY',{},'Frequency',{},'Power',{});
    nOut = 0;
    for ci = 1:numel(ids)
        C = T(string(T.ID)==ids(ci),:);
        try
            d = preprocessTrace(C);
            [feat, spec] = cwtFeatures(d, freqAxis);
            nOut = nOut + 1;
            runID=string(datetime('now','Format','yyyy-MM-dd_HHmmss'));
            resultRows{nOut} = table(ids(ci), runID, ...
                string(files(fi).name), feat.DominantFreqHz, feat.SpectralCentroid, ...
                feat.SpectralSpread, feat.SpectralEntropy, feat.TotalBandPower, ...
                feat.PeakPowerFraction, (d.time(end)-d.time(1))/60, mean(d.X), mean(d.Y), ...
                'VariableNames', {'CellID','RunID','SourceFile','DominantFreqHz', ...
                'SpectralCentroid','SpectralSpread','SpectralEntropy','TotalBandPower', ...
                'PeakPowerFraction','Duration_min','MeanX','MeanY'});
            lib(nOut) = struct('CellID',char(ids(ci)),'RunID',char(runID), ...
                'SheetName','CSV','SourceFile',files(fi).name,'MeanX',mean(d.X), ...
                'MeanY',mean(d.Y),'Frequency',spec.Frequency,'Power',spec.MeanPower);
            if makePlots
                fig=figure('Visible','off','Color','w'); plot(d.time,d.dFF,'LineWidth',1.2); xlabel('Time (s)'); ylabel('DeltaF/F'); title([files(fi).name ' | Cell ' char(ids(ci))],'Interpreter','none'); grid on;
                exportgraphics(fig,fullfile(embryoDir,'figures',[safeName(embryo) '_Cell_' safeName(ids(ci)) '.png']),'Resolution',200); close(fig);
            end
        catch ME
            warning('Skipping embryo %s cell %s: %s', embryo, ids(ci), ME.message);
        end
    end
    if nOut == 0, error('No usable cell traces in %s', inFile); end
    writetable(vertcat(resultRows{1:nOut}), fullfile(embryoDir,'tables','singlecell_results_final.csv'));
    FrequencyAxis = freqAxis; SpectralLibrary = lib; %#ok<NASGU>
    save(fullfile(embryoDir,'SpectralLibrary.mat'),'SpectralLibrary','FrequencyAxis','-v7.3');
end
end

function d = preprocessTrace(T)
vars = T.Properties.VariableNames;
timeNames = {'Time','time','Timestamp','timestamp','t','Time_s','Time (s)'};
ix = find(ismember(vars,timeNames),1);
if isempty(ix), error('No recognized time column.'); end
timeName = vars{ix}; rawTime = T.(timeName);
if isdatetime(rawTime)
    t = seconds(rawTime - rawTime(1));
elseif isduration(rawTime)
    t = seconds(rawTime);
else
    t = double(rawTime);
    unit = lower(getenv('CALCIUM_TIME_UNIT'));
    if isempty(unit), unit = 's'; end
    if contains(lower(timeName),'ms') || strcmp(unit,'ms'), t=t/1000;
    elseif contains(lower(timeName),'min') || strcmp(unit,'min'), t=t*60;
    elseif ~ismember(unit,{'s','sec','seconds'}), error('Unsupported CALCIUM_TIME_UNIT: %s',unit); end
end
ix = find(strcmpi(vars,'MeanIntensity'),1);
if isempty(ix), vnames=string(vars); ix=find(contains(vnames,'mean','IgnoreCase',true) & contains(vnames,'intens','IgnoreCase',true),1); end
if isempty(ix), error('No MeanIntensity-like column.'); end
trace = double(T{:,ix});
ix= find(strcmpi(vars,'PositionX'),1); iy=find(strcmpi(vars,'PositionY'),1);
if isempty(ix)||isempty(iy), error('Missing PositionX/PositionY.'); end
x=double(T{:,ix}); y=double(T{:,iy});
ok=isfinite(t)&isfinite(trace)&isfinite(x)&isfinite(y); t=t(ok); trace=trace(ok); x=x(ok); y=y(ok);
if numel(t)<2, error('Trace has fewer than two valid observations.'); end
if any(diff(t)<0), error('Time values are not monotonic; refusing to reorder this trace.'); end
[t,ia]=unique(t(:),'stable'); trace=trace(ia); x=x(ia); y=y(ia);
dt=diff(t); dtMed=median(dt); if ~isfinite(dtMed)||dtMed<=0, error('Invalid median sample interval.'); end
gaps=find(dt>2*dtMed); starts=[1;gaps+1]; ends=[gaps;numel(t)];
tu=[]; fu=[]; xu=[]; yu=[]; Fs=1/dtMed;
for k=1:numel(starts)
    a=starts(k); b=ends(k); if b-a+1<2, continue; end
    n=max(2,round((t(b)-t(a))/dtMed)+1); tg=linspace(t(a),t(b),n)';
    if n<10, continue; end
    tu=[tu;tg]; fu=[fu;interp1(t(a:b),trace(a:b),tg,'linear')]; %#ok<AGROW>
    xu=[xu;interp1(t(a:b),x(a:b),tg,'linear')]; yu=[yu;interp1(t(a:b),y(a:b),tg,'linear')]; %#ok<AGROW>
end
if numel(tu)<10, error('No valid continuous segment has at least 10 resampled observations.'); end
F0=max(prctile(fu,10),eps); d=struct('time',tu,'dFF',(fu-F0)/F0,'Fs',Fs,'X',xu,'Y',yu);
end

function [f,spec] = cwtFeatures(d,axisHz)
f=struct('DominantFreqHz',NaN,'SpectralCentroid',NaN,'SpectralSpread',NaN, ...
    'SpectralEntropy',NaN,'TotalBandPower',NaN,'PeakPowerFraction',NaN);
spec=struct('Frequency',axisHz,'MeanPower',[]);
if numel(d.dFF)<10, return; end
[cfs,fr]=cwt(d.dFF,d.Fs,'amor'); pRaw=mean(abs(cfs).^2,2);
[fr,ord]=sort(fr); pRaw=pRaw(ord);
p=max(interp1(fr,pRaw,axisHz,'linear',0),0); total=sum(p); spec.MeanPower=p;
if total<=0, return; end
prob=p/(total+eps); [peak,idx]=max(p);
f.DominantFreqHz=axisHz(idx); f.SpectralCentroid=sum(axisHz.*p)/(total+eps);
f.SpectralSpread=sqrt(sum(((axisHz-f.SpectralCentroid).^2).*p)/(total+eps));
f.SpectralEntropy=-sum(prob.*log2(prob+eps)); f.TotalBandPower=total;
f.PeakPowerFraction=peak/(total+eps);
end

function runCorrelations(outRoot)
eventFile=requiredEnv('CALCIUM_EVENT_FILE'); masterFile=requiredEnv('CALCIUM_CELL_MASTER'); summaryFile=requiredEnv('CALCIUM_EMBRYO_SUMMARY');
E=readtable(eventFile,'VariableNamingRule','preserve'); W=readtable(masterFile,'VariableNamingRule','preserve'); M=readtable(summaryFile,'VariableNamingRule','preserve');
need(E,{'Name','ID'}); need(W,{'EmbryoID','CellID','Group','DominantFreqHz','SpectralCentroid','SpectralSpread','SpectralEntropy','TotalBandPower','PeakPowerFraction'});
eventPrefixes={'ActivePeaksPerHr_','Total_','CellWiseActiveIncidence_','OnlySpikes_','OnlyRosettes_','OnlyWaves_'};
evVars=E.Properties.VariableNames(startsWith(E.Properties.VariableNames,'ActivePeaksPerHr_'));
if isempty(evVars), error('Event file has no ActivePeaksPerHr_* columns.'); end
pdts=erase(evVars,'ActivePeaksPerHr_');
spec={'DominantFreqHz','SpectralCentroid','SpectralSpread','SpectralEntropy','TotalBandPower','PeakPowerFraction'};
out=table(); embryos=unique(string(M.EmbryoID));
for ei=1:numel(embryos)
    id=embryos(ei); evName=unique(string(E.Name)); match=strcmpi(erase(evName,'_'),erase(id,'_'));
    if sum(match)~=1, warning('Event embryo match not unique for %s; skipping.',id); continue; end
    eName=evName(find(match,1)); e=E(string(E.Name)==eName,:); w=W(string(W.EmbryoID)==id,:);
    e.CellID=string(e.ID); w.CellID=string(w.CellID); joined=innerjoin(e,w,'Keys','CellID');
    if height(joined)<3, warning('Too few merged cells for %s.',id); continue; end
    group=string(M.Group(find(string(M.EmbryoID)==id,1)));
    for pi=1:numel(pdts)
        vars=strings(0,1);
        for q=1:numel(eventPrefixes), v=string(eventPrefixes{q})+pdts(pi); if ismember(v,joined.Properties.VariableNames), vars(end+1)=v; end, end %#ok<AGROW>
        for q=1:numel(vars)
            for r=1:numel(vars)
                out=addCorr(out,id,group,pdts(pi),"Event-Event",vars(q),vars(r),joined.(vars(q)),joined.(vars(r)));
            end
            for r=1:numel(spec), out=addCorr(out,id,group,pdts(pi),"Event-Spectral",vars(q),string(spec{r}),joined.(vars(q)),joined.(spec{r})); end
        end
        for q=1:numel(spec), for r=1:numel(spec), out=addCorr(out,id,group,pdts(pi),"Spectral-Spectral",string(spec{q}),string(spec{r}),joined.(spec{q}),joined.(spec{r})); end, end
    end
end
if isempty(out), error('No correlation results generated.'); end
writetable(out,fullfile(outRoot,'All_Embryos_All_PDT_Correlations.csv'));
end

function out=addCorr(out,id,group,pdt,kind,aName,bName,x,y)
x=numericValues(x); y=numericValues(y); ok=isfinite(x)&isfinite(y); n=sum(ok); r=NaN; p=NaN;
if n>=3 && numel(unique(x(ok)))>1 && numel(unique(y(ok)))>1
    [r,p]=corr(x(ok),y(ok),'Type','Pearson');
end
row=table(string(id),group,string(pdt),kind,string(aName),string(bName),n,r,p, ...
    'VariableNames',{'Embryo','Group','PDT','CorrelationType','Feature1','Feature2','N','PearsonR','PValue'});
if isempty(out), out=row; else, out=[out;row]; end
end

function runCorrelationSummary(outRoot)
inFile=requiredEnv('CALCIUM_CORRELATION_FILE'); T=readtable(inFile,'VariableNamingRule','preserve');
need(T,{'Group','Feature1','Feature2','PearsonR'});
event={'ActivePeaksPerHr_','Total_','CellWiseActiveIncidence_','OnlySpikes_','OnlyRosettes_','OnlyWaves_'};
spec={'DominantFreqHz','SpectralCentroid','SpectralSpread','SpectralEntropy','TotalBandPower','PeakPowerFraction'};
features=[event spec]; groups=unique(string(T.Group),'stable');
for g=1:numel(groups)
    G=T(string(T.Group)==groups(g),:); n=numel(features); med=nan(n); lo=nan(n); hi=nan(n); sd=nan(n); count=zeros(n);
    for i=1:n, for j=i:n
        ok=string(G.Feature1)==features{i}&string(G.Feature2)==features{j}; v=double(G.PearsonR(ok)); v=v(isfinite(v));
        if ~isempty(v), med(i,j)=median(v);lo(i,j)=min(v);hi(i,j)=max(v);sd(i,j)=std(v,0);count(i,j)=numel(v); end
        med(j,i)=med(i,j);lo(j,i)=lo(i,j);hi(j,i)=hi(i,j);sd(j,i)=sd(i,j);count(j,i)=count(i,j);
    end,end
    A=array2table(med,'VariableNames',matlab.lang.makeValidName(features),'RowNames',features); writetable(A,fullfile(outRoot,safeName(groups(g))+"_Median_R.csv"),'WriteRowNames',true);
    A=array2table(lo,'VariableNames',matlab.lang.makeValidName(features),'RowNames',features); writetable(A,fullfile(outRoot,safeName(groups(g))+"_Min_R.csv"),'WriteRowNames',true);
    A=array2table(hi,'VariableNames',matlab.lang.makeValidName(features),'RowNames',features); writetable(A,fullfile(outRoot,safeName(groups(g))+"_Max_R.csv"),'WriteRowNames',true);
    A=array2table(sd,'VariableNames',matlab.lang.makeValidName(features),'RowNames',features); writetable(A,fullfile(outRoot,safeName(groups(g))+"_Std_R.csv"),'WriteRowNames',true);
    A=array2table(count,'VariableNames',matlab.lang.makeValidName(features),'RowNames',features); writetable(A,fullfile(outRoot,safeName(groups(g))+"_N_Correlations.csv"),'WriteRowNames',true);
    fig=figure('Visible','off','Color','w'); imagesc(med,[-1 1]); axis square; colorbar; colormap(parula); xticks(1:n);yticks(1:n);xticklabels(features);yticklabels(features);xtickangle(45); title(groups(g),'Interpreter','none');
    for i=1:n, for j=1:n, if isfinite(med(i,j)), text(j,i,sprintf('%.2f [%.2f, %.2f]',med(i,j),lo(i,j),hi(i,j)),'HorizontalAlignment','center','FontSize',5); end,end,end
    exportgraphics(fig,fullfile(outRoot,safeName(groups(g))+"_Median_Correlation_Matrix.png"),'Resolution',200); close(fig);
end
end

function runIndividualMaps(outRoot)
boundaryRoot=requiredEnv('CALCIUM_BOUNDARY_ROOT'); masterFile=requiredEnv('CALCIUM_CELL_MASTER');
features={'DominantFreqHz','SpectralCentroid','SpectralSpread','SpectralEntropy','TotalBandPower','PeakPowerFraction'};
B=dir(fullfile(boundaryRoot,'**','*-Boundary_SuperEmbryo.csv')); W=readtable(masterFile,'VariableNamingRule','preserve');
need(W,{'EmbryoID','CellID','Group','MeanX','MeanY',features{:}});
if isempty(B), error('No *-Boundary_SuperEmbryo.csv files found under %s.',boundaryRoot); end
AllMerged=cell(numel(B),1); AllX=[]; AllY=[]; featureValues=struct();
for f=1:numel(features), featureValues.(features{f})=[]; end
for k=1:numel(B)
    id=erase(B(k).name,'-Boundary_SuperEmbryo.csv'); b=readtable(fullfile(B(k).folder,B(k).name),'VariableNamingRule','preserve');
    if ~all(ismember({'Calcium_ID','Centroid_1','Centroid_2'},b.Properties.VariableNames)), warning('Boundary columns missing for %s.',id); continue; end
    wid=string(W.EmbryoID); w=W(wid==id | startsWith(wid,id+"_"),:); if isempty(w), warning('No cell rows matched %s.',id); continue; end
    b.Calcium_ID=string(b.Calcium_ID); w.CellID=string(w.CellID); m=innerjoin(b,w,'LeftKeys','Calcium_ID','RightKeys','CellID');
    if isempty(m), warning('No matched cells for %s.',id); continue; end
    AllMerged{k}=m; ok=isfinite(m.Centroid_1)&isfinite(m.Centroid_2); AllX=[AllX;double(m.Centroid_1(ok))]; AllY=[AllY;double(m.Centroid_2(ok))]; %#ok<AGROW>
    for f=1:numel(features), v=double(m.(features{f})); featureValues.(features{f})=[featureValues.(features{f});v(isfinite(v))]; end %#ok<AGROW>
end
if isempty(AllX)||isempty(AllY), error('No boundary/spectral rows with valid spatial coordinates matched.'); end
xm=.02*(max(AllX)-min(AllX)); ym=.02*(max(AllY)-min(AllY)); xlimGlobal=[min(AllX)-xm max(AllX)+xm]; ylimGlobal=[min(AllY)-ym max(AllY)+ym]; colorLimits=struct();
for f=1:numel(features), v=featureValues.(features{f}); if ~isempty(v), q=prctile(v,[1 99]); if q(1)==q(2), q=q+[-1 1]*max(abs(q(1))*0.05,1e-6); end; colorLimits.(features{f})=q; else, colorLimits.(features{f})=[0 1]; end, end
for k=1:numel(B)
    m=AllMerged{k}; if isempty(m), continue; end
    id=erase(B(k).name,'-Boundary_SuperEmbryo.csv'); folder=fullfile(outRoot,safeName(id)); if isfolder(folder), error('Duplicate embryo output folder would overwrite: %s',folder); end; mkdir(folder); writetable(m,fullfile(folder,[id '_Spectral_Mapped.csv']));
    for f=1:numel(features)
        if isempty(featureValues.(features{f})), continue; end
        fig=figure('Visible','off','Color','w'); scatter(m.Centroid_1,m.Centroid_2,25,m.(features{f}),'filled'); axis equal; xlim(xlimGlobal); ylim(ylimGlobal); clim(colorLimits.(features{f})); colorbar; xlabel('Template X');ylabel('Template Y');title([id ' — ' features{f}],'Interpreter','none');
        exportgraphics(fig,fullfile(folder,[id '_Mapped_' features{f} '.png']),'Resolution',300); close(fig);
    end
end
end

function runCombineMaps(outRoot)
mapRoot=requiredEnv('CALCIUM_MAPPED_ROOT'); summary=readtable(requiredEnv('CALCIUM_EMBRYO_SUMMARY'),'VariableNamingRule','preserve');
features={'DominantFreqHz','SpectralCentroid','SpectralSpread','SpectralEntropy','TotalBandPower','PeakPowerFraction'};
groups=unique(string(summary.Group),'stable');
for f=1:numel(features)
    pdf=fullfile(outRoot,[features{f} '_AllEmbryos.pdf']); first=true;
    for g=1:numel(groups)
        ids=string(summary.EmbryoID(string(summary.Group)==groups(g))); files=strings(0,1);
        for e=1:numel(ids)
            q=dir(fullfile(mapRoot,'**',[char(ids(e)) '_Mapped_' features{f} '.png'])); if ~isempty(q), files(end+1)=fullfile(q(1).folder,q(1).name); end %#ok<AGROW>
        end
        for st=1:16:numel(files)
            fig=figure('Visible','off','Color','w','Position',[100 100 1500 1100]); tiledlayout(4,4,'TileSpacing','compact','Padding','compact'); sgtitle([features{f} ' — ' char(groups(g))],'Interpreter','none');
            for i=st:min(st+15,numel(files)), nexttile; imshow(files(i)); title(erase(string(files(i)),"_Mapped_"+features{f}+".png"),'Interpreter','none','FontSize',8); end
            if first, exportgraphics(fig,pdf,'ContentType','image','Resolution',200); first=false; else, exportgraphics(fig,pdf,'ContentType','image','Resolution',200,'Append',true); end
            close(fig);
        end
    end
end
end

function runGroupSpatial(outRoot)
mappedRoot=requiredEnv('CALCIUM_MAPPED_ROOT'); assign=readtable(requiredEnv('CALCIUM_GROUP_ASSIGNMENTS'),'VariableNamingRule','preserve');
need(assign,{'Embryo','Group'}); features={'DominantFreqHz','SpectralEntropy','SpectralCentroid','SpectralSpread','TotalBandPower','PeakPowerFraction'};
cnpuSize=9; originX=120; originY=1; xRange=[120 1122]; yRange=[1 1190];
files=dir(fullfile(mappedRoot,'**','*_Spectral_Mapped.csv')); if isempty(files), error('No *_Spectral_Mapped.csv files found.'); end
AllCells=table();
for k=1:numel(files)
    T=readtable(fullfile(files(k).folder,files(k).name),'VariableNamingRule','preserve'); id=erase(files(k).name,'_Spectral_Mapped.csv');
    loc=find(strcmpi(string(assign.Embryo),id)); if numel(loc)~=1, warning('No unique group assignment for %s; skipping.',id); continue; end
    need(T,{'Centroid_1','Centroid_2',features{:}}); T.Embryo=repmat(string(id),height(T),1); T.Group=repmat(string(assign.Group(loc)),height(T),1); if isempty(AllCells), AllCells=T; else, AllCells=[AllCells;T]; end %#ok<AGROW>
end
if isempty(AllCells), error('No mapped rows had unique group assignments.'); end
x=double(AllCells.Centroid_1); y=double(AllCells.Centroid_2); ok=isfinite(x)&isfinite(y); AllCells=AllCells(ok,:); x=x(ok); y=y(ok);
AllCells.CNPU_X=floor((x-originX)/cnpuSize); AllCells.CNPU_Y=floor((y-originY)/cnpuSize); AllCells.CNPU_ID=string(AllCells.CNPU_X)+"_"+string(AllCells.CNPU_Y);
writetable(AllCells,fullfile(outRoot,'All_MappedCells_with_CNPU.csv'));
gx=0:ceil((xRange(2)-originX)/cnpuSize)-1; gy=0:ceil((yRange(2)-originY)/cnpuSize)-1; [xx,yy]=ndgrid(gx,gy); grid=table(xx(:),yy(:),'VariableNames',{'CNPU_X','CNPU_Y'});
grid.CNPU_StartX=originX+grid.CNPU_X*cnpuSize; grid.CNPU_EndX=grid.CNPU_StartX+cnpuSize;
grid.CNPU_StartY=originY+grid.CNPU_Y*cnpuSize; grid.CNPU_EndY=grid.CNPU_StartY+cnpuSize;
grid.CNPU_CenterX=originX+grid.CNPU_X*cnpuSize+cnpuSize/2; grid.CNPU_CenterY=originY+grid.CNPU_Y*cnpuSize+cnpuSize/2;
boundary13=readBoundaryXY(requiredEnv('CALCIUM_STAGE13_SOX2_BOUNDARY'));
boundary14=readBoundaryXY(requiredEnv('CALCIUM_STAGE14_SOX2_BOUNDARY'));
grid.Stage13Template=inpolygon(grid.CNPU_CenterX,grid.CNPU_CenterY,boundary13(:,1),boundary13(:,2));
grid.Stage14Template=inpolygon(grid.CNPU_CenterX,grid.CNPU_CenterY,boundary14(:,1),boundary14(:,2));
grid.AnySox2Template=grid.Stage13Template|grid.Stage14Template;
groups=unique(string(AllCells.Group),'stable'); summaryRows=table();
grid.Observed=false(height(grid),1); grid.NCells=zeros(height(grid),1);
[obsKey,~,obsIdx]=unique([AllCells.CNPU_X AllCells.CNPU_Y],'rows'); obsN=accumarray(obsIdx,1);
[isObs,obsLoc]=ismember([grid.CNPU_X grid.CNPU_Y],obsKey,'rows'); grid.Observed=isObs; grid.NCells(isObs)=obsN(obsLoc(isObs));
writetable(grid,fullfile(outRoot,'Complete_CNPU_Grid.csv'));
fig=figure('Visible','off','Color','w'); tiledlayout(1,3,'TileSpacing','compact');
nexttile; scatter(grid.CNPU_CenterX,grid.CNPU_CenterY,8,double(grid.Observed),'filled'); hold on; plot(boundary13(:,1),boundary13(:,2),':k'); plot(boundary14(:,1),boundary14(:,2),'--k'); hold off; axis equal; set(gca,'YDir','reverse'); title('Observed CNPUs and templates');
nexttile; scatter(grid.CNPU_CenterX,grid.CNPU_CenterY,8,double(grid.Stage13Template),'filled'); axis equal; set(gca,'YDir','reverse'); title('Stage 13 template (Stage 14 boundary)');
nexttile; scatter(grid.CNPU_CenterX,grid.CNPU_CenterY,8,double(grid.Stage14Template),'filled'); axis equal; set(gca,'YDir','reverse'); title('Stage 14 template (Stage 15 boundary)');
exportgraphics(fig,fullfile(outRoot,'CompleteCNPUGrid_QC.png'),'Resolution',250); close(fig);
for f=1:numel(features)
  feature=features{f}; for g=1:numel(groups)
    G=AllCells(string(AllCells.Group)==groups(g),:); embryos=unique(string(G.Embryo)); per=table();
    for e=1:numel(embryos)
      E=G(string(G.Embryo)==embryos(e),:); vals=double(E.(feature)); valid=isfinite(vals); E=E(valid,:); if isempty(E), continue; end
      [keys,~,ix]=unique([E.CNPU_X E.CNPU_Y],'rows'); med=splitapply(@median,double(E.(feature)),ix); n=splitapply(@numel,double(E.(feature)),ix);
      part=table(repmat(embryos(e),size(keys,1),1),keys(:,1),keys(:,2),med,n,'VariableNames',{'Embryo','CNPU_X','CNPU_Y','EmbryoMedian','NCells'}); if isempty(per), per=part; else, per=[per;part]; end %#ok<AGROW>
    end
    for i=1:height(grid)
      z=per(per.CNPU_X==grid.CNPU_X(i)&per.CNPU_Y==grid.CNPU_Y(i),:); if isempty(z), gm=NaN; gs=NaN; ne=0; nc=0; else, gm=median(z.EmbryoMedian); gs=std(z.EmbryoMedian,0); ne=height(z); nc=sum(z.NCells); end
      row=table(groups(g),string(feature),grid.CNPU_X(i),grid.CNPU_Y(i),grid.CNPU_CenterX(i),grid.CNPU_CenterY(i),gm,gs,ne,100*ne/numel(embryos),nc,'VariableNames',{'Group','Feature','CNPU_X','CNPU_Y','CNPU_CenterX','CNPU_CenterY','GroupMedian','GroupSD','NEmbryos','CoveragePct','NCells'}); if isempty(summaryRows), summaryRows=row; else, summaryRows=[summaryRows;row]; end %#ok<AGROW>
    end
  end
end
summaryRows.SmoothedMedian=nan(height(summaryRows),1); summaryRows.SmoothedZ=nan(height(summaryRows),1); halfWindow=3;
for f=1:numel(features), for g=1:numel(groups)
    loc=find(string(summaryRows.Feature)==features{f}&string(summaryRows.Group)==groups(g)); coords=[summaryRows.CNPU_X(loc),summaryRows.CNPU_Y(loc)]; values=summaryRows.GroupMedian(loc); sums=zeros(numel(loc),1); counts=zeros(numel(loc),1);
    for dx=-halfWindow:halfWindow, for dy=-halfWindow:halfWindow
        [found,at]=ismember([coords(:,1)+dx,coords(:,2)+dy],coords,'rows'); good=found; good(found)=isfinite(values(at(found))); sums(good)=sums(good)+values(at(good)); counts(good)=counts(good)+1;
    end,end
    sm=sums./counts; sm(counts==0)=NaN; valid=isfinite(sm); summaryRows.SmoothedMedian(loc)=sm; if any(valid), mu=mean(sm(valid)); sigma=std(sm(valid),0); if sigma>0, z=nan(size(sm)); z(valid)=(sm(valid)-mu)./sigma; summaryRows.SmoothedZ(loc)=z; end, end
end,end
writetable(summaryRows,fullfile(outRoot,'Group_CNPU_Summaries.csv'));
for f=1:numel(features)
  for metric={'GroupMedian','GroupSD','CoveragePct'}
    R=summaryRows(string(summaryRows.Feature)==features{f},:); fig=figure('Visible','off','Color','w'); tiledlayout(2,2,'TileSpacing','compact'); lim=[min(R.(metric{1}),[],'omitnan'),max(R.(metric{1}),[],'omitnan')]; if any(~isfinite(lim))||lim(1)==lim(2),lim=[0 1];end
    for g=1:numel(groups), nexttile; G=R(string(R.Group)==groups(g),:); scatter(G.CNPU_CenterX,G.CNPU_CenterY,25,G.(metric{1}),'filled'); hold on; if contains(lower(groups(g)),'13'), plot(boundary13(:,1),boundary13(:,2),':k','LineWidth',1.2); else, plot(boundary14(:,1),boundary14(:,2),':k','LineWidth',1.2); end; hold off; axis equal; xlim([min(grid.CNPU_StartX) max(grid.CNPU_EndX)]); ylim([min(grid.CNPU_StartY) max(grid.CNPU_EndY)]); set(gca,'YDir','reverse'); clim(lim); title(groups(g),'Interpreter','none'); colorbar; end
    sgtitle([features{f} ' — ' metric{1}],'Interpreter','none'); exportgraphics(fig,fullfile(outRoot,[features{f} '_' metric{1} '_GroupSpatialMaps.png']),'Resolution',250); close(fig);
  end
end
end

function runBoundaryQC(outRoot)
boundaryRoot=requiredEnv('CALCIUM_BOUNDARY_ROOT'); files=dir(fullfile(boundaryRoot,'**','*-Boundary_SuperEmbryo.csv')); out=table();
for k=1:numel(files)
    T=readtable(fullfile(files(k).folder,files(k).name),'VariableNamingRule','preserve'); vars=T.Properties.VariableNames;
    ix=find(strcmpi(vars,'Centroid_1'),1); iy=find(strcmpi(vars,'Centroid_2'),1); if isempty(ix)||isempty(iy), warning('No centroid columns in %s',files(k).name); continue; end
    bvars=vars(startsWith(vars,'Boundary')); nums=cellfun(@(s)sscanf(s,'Boundary_%d'),bvars); [~,ord]=sort(nums); bvars=bvars(ord); validCells=0;
    for r=1:height(T), b=table2array(T(r,bvars)); b=b(isfinite(b)); if numel(b)>=6 && mod(numel(b),2)==0, validCells=validCells+1; end, end
    row=table(string(files(k).name),height(T),validCells,sum(isfinite(T{:,ix})&isfinite(T{:,iy})),'VariableNames',{'BoundaryFile','Rows','CellsWithBoundary','CellsWithCentroid'}); if isempty(out), out=row; else, out=[out;row]; end %#ok<AGROW>
end
writetable(out,fullfile(outRoot,'Boundary_QC.csv'));
end

function need(T, names)
missing=names(~ismember(names,T.Properties.VariableNames)); if ~isempty(missing), error('Missing required columns: %s',strjoin(missing,', ')); end
end
function x=numericValues(x)
if isnumeric(x), x=double(x); else, x=str2double(string(x)); end
end
function value=requiredEnv(name)
value=getenv(name); if isempty(value), error('Required environment variable is unset: %s',name); end
if ~isfolder(value)&&~isfile(value), error('Configured input path does not exist: %s',value); end
end
function s=safeName(s)
s=regexprep(char(string(s)),'[^A-Za-z0-9_.-]','_');
end
function XY=readBoundaryXY(path)
T=readtable(path,'VariableNamingRule','preserve'); vars=T.Properties.VariableNames;
xc={'X','Centroid_1','Boundary_1','V1','Var1'}; yc={'Y','Centroid_2','Boundary_2','V2','Var2'}; ix=[]; iy=[];
for k=1:numel(xc), q=find(strcmpi(vars,xc{k}),1); if ~isempty(q), ix=q; break; end, end
for k=1:numel(yc), q=find(strcmpi(vars,yc{k}),1); if ~isempty(q), iy=q; break; end, end
if (isempty(ix)||isempty(iy))&&width(T)==2, ix=1; iy=2; end
if isempty(ix)||isempty(iy), error('Could not find X/Y columns in boundary file: %s',path); end
XY=[double(T{:,ix}),double(T{:,iy})]; XY=XY(all(isfinite(XY),2),:); if size(XY,1)<3, error('Boundary has fewer than 3 valid points: %s',path); end
if any(XY(1,:)~=XY(end,:)), XY(end+1,:)=XY(1,:); end
end
