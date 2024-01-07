str%
% merges DTW time segments data
%  and includes all other times series (e.g. RGR, solar forcing)
%
% kai wirtz (HZG) 2023
%
clear all; close all;
load_pars; % sets common parameters (scdir, cc, latlim, regs)
timelim = [2.8 10.2];

% parameters for band-pass filtering
tmov = 151;toff = 40; % smoothing parameters
tavg = 1.5;         % detrend

% ------ common time vector and data matrix
dt  = 0.01; di = 3;
time= timelim(1):dt:timelim(2);
dat = zeros(length(time),12) + NaN;
dat(:,1) = time; legdat{1} = 'time';

% ------------------------------------------------------------------
%  RGR Europe (methods)
load([scdir 'avg_rgr_all']);   % variables: leg  rgr_m  tirgr
ri = [6]; % 6 area-based
          % other methods can be added such as w normalized calbration or pooled
for i = 1:length(ri)
  j = ri(i);
  legdat{1+i} = leg{j};
  % re-grid on common time vector
  it = 1:min(find( time >= tirgr(end)));
  ts = interp1(tirgr,rgr_m(:,j),time(it),'linear','extrap');
  ts = movweighavg(time(it)*1E3,ts,tmov,toff); %
  dat(it,1+i) = ts;
end
legdat{3} = 'void';
i0=4;
% ------------------------------------------------------------------
%  RGR East Asia (China+)
file = [scdir 'mat/AllPop_EA_NoNorm_Bin100_all.mat'];
if exist(file)
  load(file); %poptime = tm, ymv = ymv,trgr = tirgr,rgr = rgrv,nreg = nregions
  trgr = trgr*1E-3;
  trgr = flipud(trgr); rgr = flipud(rgr);
  % bring both rgr estimates on same timeline
  it = find(time >= (trgr(1)) & time <= (trgr(end)) );
  ts = interp1(trgr,rgr,time(it),'linear','extrap');
  ts = movweighavg(time(it)*1E3,ts,tmov,toff)*1E3; %
  dat(it,i0)  =  movweighavg(time(it)*1E3,ts,tmov,toff);  %2nd smooth
end
legdat{i0} = ['RGR East Asia'];% (detrend)
i0=i0+1;

% ------------------------------------------------------------------
% RGR South America
load(['data/SA_spd_rgr']); %,'sa_rtim','spd1','sa_rgr'
it = find(time >= sa_rtim(1) & time <= sa_rtim(end)); %sa_time','spd1','sa_rgr
ts = interp1(sa_rtim,sa_rgr,time(it),'linear','extrap');
ts = movweighavg(time(it)*1E3,ts,tmov,toff); %
dat(it,i0) = ts; % 5
legdat{i0} = 'RGR South America';
i0=i0+1;

% ------------------------------------------------------------------
% Total Solar Irradiance (TSI) from Steinhilber et al 2012
ts = load(['data/Steinhilber2012_Solar.dat']);%
ts_time = ts(:,1);
it = 1:min(find(time >= ts_time(end)));
ts_m = interp1(ts_time,ts(:,2),time(it),'linear','extrap');
[ut avg1500] = movavg(time(it),ts_m,1.);
ts_m = ts_m-avg1500;
ts_m = ts_m/nanstd(ts_m(find(time(it)<8.2)));
dat(it,i0) = -movweighavg(time(it)*1E3,ts_m,tmov,toff); % revert in sign!
legdat{i0} = '- TSI';  % 6
i0=i0+1;

% ---------------------------------------------------------------
%  Northern Irish bog data provided by Rowan McLaughlin
load(['data/bog_std']); % variables: ,'ti',ringw_sm; bogstd_space; bogstd_time;
it = find(time <= ti(end) & time >= ti(1) );
ts  = interp1(ti,ringw_sm,time(it),'linear','extrap');
ts0 = interp1(ti,bogstd_space,time(it),'linear','extrap');
ts1 = interp1(ti,bogstd_time,time(it),'linear','extrap');
ts2 = interp1(ti,flood_m,time(it),'linear','extrap');
dat(it,i0) = (ts-mean(ts))/std(ts);
dat(it,i0+1) = -(ts0-mean(ts0))/std(ts0);
%%dat(it,i0+2) = -(ts1-mean(ts1))/std(ts1);
%%dat(it,i0+3) = (ts2-nanmean(ts2))/nanstd(ts2);
legdat{i0} = ['bog tree ring width'];
legdat{i0+1} = ['tree growth homogeneity'];
%%legdat{i0+2} = ['-STD tree growth (time)'];
%%legdat{i0+3} = ['floods'];
i0=i0+2;

% ---------------------------------------------------------------
% merge dtw & pca results for two time-segments
ri = [7]; jj=1;
split = 2; dtw_Dist_crit = 90;
tol = 0.1;
% clear all fields
datt  = zeros(split,length(time),length(ri))+NaN;
datm  = zeros(split,length(time))+NaN;
datall= zeros(1,length(time))+NaN;
datm1 = datm;datall1 = datall;datm2 = datm;datall2 = datall;
pca   = zeros(5,length(time))+NaN;;
datma = zeros(5,split,length(time))+NaN;
tmov2 = 181;
ln = ri(jj);
% loop over time-segments
tj = 1;
for tmax = [6.2 9.5]
  file = sprintf('%sdtwpca/dtwpca2_%3.2f_%2.0f_%1.0f_%d.mat',scdir,tmax,dtw_Dist_crit,tol*100,ln);
  fprintf('tm = %1.1f %d\t loading dtw-pca from %s\n',tmax,ln,file);
  if exist(file)
    load(file);
    exp_vara(tj,:) = exp_var;
    it = max(find(time <= (t_change(1)))):min(find(time >= (t_change(end))));
  %  pca_cm  =  movweighavg(t_change*1E3,pca_change,tmov2,toff);
  %  pca_cm  =  movweighavg(t_change*1E3,pca_cm,tmov/2,toff);
    pca_cm  =  movweighavg(t_change*1E3,pca_change,tmov,toff);%130,20
    datt(tj,it,j)  =  interp1(t_change,pca_cm,time(it),'linear','extrap');
    for jp = 1:5
      pcat  =  movweighavg(time20*1E3,pcats{jp} ,tmov,toff);%130,20
      datma(jp,tj,it)  =  interp1(time20,pcat,time(it),'linear','extrap');
    end %jp
    datm(tj,it)  =  nanmean(datt(tj,it,:),3);
    itv{tj} = it;
  end % if exists
  tj = tj+1;
end % for
% 2nd loop over time-segments
for tj = 1:2
  datall(itv{tj}) = datm(tj,itv{tj});
  for jp = 1:5
    pca(jp,itv{tj}) = datma(jp,tj,itv{tj});
  end
end
ii = intersect(itv{1},itv{2}); nt = length(ii);
ff = (ii(end)-ii)/(nt-1);
datm(1,ii(end)) = 0; datm(2,ii(1)) = 0;
datall(ii) = ff.*datm(1,ii)+(1-ff).*datm(2,ii);
for jp = 1:5
  datma(jp,1,ii(end))  =  0;
  datma(jp,2,ii(1))  =  0;
  pca(jp,ii) = ff.*squeeze(datma(jp,1,ii))'+(1-ff).*squeeze(datma(jp,2,ii))';
end
[ut avg1500] = movavg(time,datall,1);
datall  =  datall-avg1500;
clim_stability  =  datall/nanstd(datall);
dat(:,i0)  =  clim_stability;
for jp = 1:2
  dat(:,i0+jp)  =   pca(jp,:)/nanstd(pca(jp,:));
  legdat{i0+jp} = ['climate PCA' num2str(jp)];
end
legdat{i0} = 'climate stability';
file = sprintf('%sdtwpca/dtwpca_proxydata_%3.2f_%2.0f_%1.0f_%d.mat',scdir,tmax,dtw_Dist_crit,tol*100,ln);
file = sprintf('%spca_0.mat',scdir);
fprintf('save PCs/dPCs data in %s\n',file)
save('-v6',file,'time','pca','exp_vara','clim_stability');

% --------------------------------------
% save combined data to file
file = sprintf('%starget_ts_0.mat', scdir);
fprintf('save all data in %s\n',file)
save('-v6',file,'dat','legdat','tmov','toff');
