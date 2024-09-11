% This script calculates Dynamic Time Warping (DTW) of paleoclimate proxy time-series.
% It processes time-series data, performs Principal Component Analysis (PCA),
% and saves the results, including DTW distances and PCA components, into a MATLAB binary file.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input Files:
% - paleoclim/InEur_27_920.mat: MATLAB file containing paleoclimate proxy time-series
% - paleoclim/proxydescription_488_0_11.00_1.70.mat: MATLAB file containing meta-info for proxies

% Output Files:
% - out/dtwpca/dtwpca2_<timeLimits(2)>_<dtw_Dist_crit>_<tol*100>_<loop>.mat: MATLAB binary file containing DTW and PCA results

% Variables:
% - dtw_Dist_crit: DTW distance criterion
% - tol, tol2, dtol: Time tolerances
% - out: Flag for graphical output
% - nmax: Maximum number of principal components in PCA
% - tfac: Time factor for adjustments
% - ndbin: Number of bins for distance calculations
% - ppcol: Color map for plotting
% - fs: Font size for plotting
% - lw: Line width for plotting
% - ncol, nrow: Number of columns and rows for plot display
% - dxp, dyp: Plot display dimensions
% - markerColors: Colors for markers in plots
% - dt0: Time step for time vectors
% - timeLimitsv: Time limits for segments
% - tsmax: Maximum time segment
% - file, evinfofile: Filenames for input data
% - timeLimits: Current time limits for processing
% - time2, time20: Time vectors
% - nt: Length of time vector
% - tmin, tmax: Minimum and maximum times for time-series
% - numprox: Number of proxies
% - lons, lats: Longitude and latitude of proxies
% - vli: Length of InEur vector
% - wps, match: Weight and match vectors
% - matchInd: Match index matrix
% - cDTW_Dist, cDTW_sign, cPears, cLocD: Matrices for DTW distances, signs, Pearson correlations, and local distances
% - iv, ei: Indices for iteration and event
% - datprox: Data matrix for proxies
% - i, tim, ti, val: Indices and values for time-series
% - scs, ta: Strings for proxy names and tags
% - tmax2, tmin2: Maximum and minimum times for proxies
% - it, tiu, ia, ic: Indices for time vectors
% - time1, proxi, time0, time1b: Time and proxy data vectors
% - x2: Interpolated and normalized proxy data
% - dpi: Index vector for proxies
% - pca: PCA model
% - we_pca, we_dtw: Weight vectors for PCA and DTW
% - loop, out2: Loop index and output flag
% - nf: Figure index for plotting
% - datproxw: Weighted data matrix for proxies
% - dtw_weigh, dtw_weighi: Weight vectors for DTW
% - np: Number of valid time-series
% - iy, ix: Indices for plot positions
% - cs, csd, cDTW_D, cDTW_k: Cell arrays for correlation, distance, DTW distances, and warping paths
% - dists: Distance vector
% - deltat, deltaw: Delta time and weight vectors
% - i2, j, dist: Indices and distance for compared time-series
% - tmax, tmin: Maximum and minimum times for compared time-series
% - itj: Indices for time vectors
% - dtw_Dist, r, wp, w, li, sign: DTW distance, correlation, warping path, length, and sign
% - tiu, ia, ic: Indices for unique time vectors
% - tio, val: Time and value vectors for interpolation
% - weigh: Weight for DTW
% - deltat, deltaw: Delta time and weight vectors
% - itj0, itj1: Indices for time vectors
% - cDTW_w, cDTW_sign, cPears, cLocD: Cell arrays for DTW warping paths, signs, Pearson correlations, and local distances
% - ti, x2: Time and proxy data vectors
% - ind, vs, vi: Indices and sorted values for DTW distances
% - mDTW, mPear: Mean DTW distances and Pearson correlations
% - ii: Indices for delta weights
% - del, delm, ff: Delta values and factors
% - delmax: Maximum delta values
% - dilt: Delta time vector
% - le: Line handles for plotting
% - dt: Delta time vector for sorting
% - nl: Length of delta time vector
% - pca0: Initial PCA model
% - exp_var: Explained variance for PCA
% - pwei: Weight vector for PCA
% - m, ii: Sorted weights and indices for PCA
% - x: Normalized weights for PCA
% - we_pca: Weight vector for PCA
% - mlPear: Mean Pearson correlations
% - pca_ts: PCA time-series
% - kine: Kinetic energy for PCA stability
% - pca_sm: Smoothed PCA time-series
% - pcats: Cell array for PCA time-series
% - dpca: Delta PCA time-series
% - pca_change: PCA change vector
% - t_change: Time change vector
% - file: Filename for output data

clear all;
addpath('~/tools/mlabtools/pca_toolbox_1.5')

% settings
dtw_Dist_crit = 90;  % approximately the average
tol = 0.1;           % time tolerance in ka (thus 100yr)
tol2= 2*tol; dtol=tol/2;
fprintf('dtw_Dist_crit=%1.0f tol=%1.3f\n',dtw_Dist_crit,tol)

out=0;    % 1: graphical output
nmax=5;   % max number of PCs in PCA
tfac=0.2; %0.5;
ndbin=50;

% plot settings
ppcol = prism(7);ppcol(3,:)=[];
fs    = 20; lw=1;
ncol  = 4; % plot display: number of columns
nrow  = 4; % plot display: number of rows
dxp=0.9/ncol;dyp=0.89/nrow;
markerColors(1,:)=[1 0.5 0]; markerColors(2,:)=[0 0.5 1];

% time segments
dt0=0.01;
timeLimitsv=[[2.8 6.2];[5.8 9.5]];
tsmax=3.;

%load_pars;
outputDirectory='out/';
edir ='paleoclim/';

% load paleoclimate proxy time-series and meta-info
file=sprintf('%sInEur_%d_%03.0f.mat',edir,27,920);
load(file);
evinfofile=sprintf('%sproxydescription_488_0_%03.2f_%03.2f.mat',edir,11.,1.7);
load(evinfofile); % into struct evinfo

% -------------------------------------------------------
% loop over time segments
for tl=1:2
  timeLimits=timeLimitsv(tl,:);
  % time vectors
  time2=(timeLimits(1)-tol):dt0:(timeLimits(2)+tol);
  time20=(timeLimits(1)):dt0:(timeLimits(2));
  nt=length(time2)-1;

  % time range of TS
  tmin=evinfo.t_min(InEur); tmax=evinfo.t_max(InEur);
  numprox=length(InEur);

  % geographical position
  lons=evinfo.Longitude(InEur); lats=evinfo.Latitude(InEur);vli=length(InEur);

  % -------------------------------------------------------
  % initialize and clear variables/vectors
  wps       = zeros(numprox,1); match=wps;
  matchInd  = zeros(1+2,numprox);
  cDTW_Dist=zeros(numprox,numprox);
  cDTW_sign = cDTW_Dist; cPears=cDTW_Dist;
  cLocD     = cDTW_Dist;
  iv=1; ei=1;
  clear datprox;

  % -------------------------------------------------------
  % loop over proxy TS

  for ie=1:numprox
    i=InEur(ie);
    % dates and value of time-series (TS)
    tim=evinfo.time{i};  ti=tim;
    val=evinfo.value{i};

    % create name tag
    scs=evinfo.Proxy{i};
    scs(regexp(scs,'[$]'))=[];%%scs(regexp(scs,'\'))=[]; %'
    ta{ie}=sprintf('%s (%s)',evinfo.Plotname{i},scs);

    % check for upper overlap with time segment
    if(max(ti)>time20(end)-tol2 & max(ti)<time20(end))
      %%fprintf('%d %1.3f %1.3f\t',ie,min(ti),max(ti))
      ti = reshape(ti,1,length(ti));
      val= reshape(val,1,length(val));

      % extrapolation: small extension at the boundary
      if(ti(2)<ti(1)) % direction of time vector
        ti=[time20(end) ti]; val=[val(end) val];
      else
        ti=[ti time20(end)]; val=[val val(end)];
      end
      %%fprintf('%1.3f %1.3f\n',min(ti),max(ti))
    end

    % set time span of TS
    tmax2{ie}=max(ti)+tol;tmin2{ie}=min(ti)-tol;

    % sort dates
    it=find(time2>=tmin2{ie} & time2<=tmax2{ie});
    [tiu,ia,ic]=unique(ti);
    time1{ie}=time2(it);

    % interpolate on sorted time grid
    proxi{ie}=interp1(tiu,val(ia),time1{ie},'linear','extrap');

    % initialize date vectors
    time0{ie} = time1{ie};
    time1b{ie}= time1{ie};

    % check for overlap with time segment
    if(min(ti)-tol<time20(1) & max(ti)+tol>time20(end))
      % interpolate on sorted time grid
       x2=interp1(time1{ie},proxi{ie},time20,'linear','extrap');
      % normalize
       x2=x2-mean(x2); x2=x2/std(x2);

       datprox(iv,:)=x2;
       dpi(iv)=ie;
       iv=iv+1;
    else
       fprintf('skip %d\t%1.3f %1.3f\n',ie,min(ti),max(ti));
    end
  end

  % -------------------------------------------------------
  % calc initial PCA
  np=iv-1; % number of valid TS

  pca = pca_model(datprox(1:np,:)',nmax,'auto');%'
  fprintf('\npca: %1.3f %1.3f\n',pca.cum_var(2),pca.cum_var(nmax));

  we_pca = ones(numprox,1); we_dtw = ones(numprox,1);
  % -------------------------------------------------------
  % loop over iterations
  for loop=1:8
    out2=(loop==-1)+(loop==-2);
  % initialize and clear variables/vectors
    ie0=0;nf=1;iv=1;
    clear datprox;
    clear datproxw;
    dtw_weigh=zeros(numprox,1); dtw_weighi=dtw_weigh;

    % -------------------------------------------------------
    % loop over primer TS
    for ie=1:numprox
      np=mod(ie0,ncol*nrow);
      if out
        if(np==0)
          if(loop==1)
            gcf(nf)=figure(nf);clf;
            set(gcf(nf),'position',[20+nf*15 1 1640 990],'Color','w','Visible', 'off');
          end
          set(0, 'CurrentFigure', gcf(nf))
          %figure(nf);
          set(gcf(nf),'Visible','off');
          nf=nf+1;
        end
      end
      ie0=ie0+1;
      % -------------------------------------------------------
      % plot proxy
      iy = floor(np/ncol); ix = mod(np,ncol);

      it=find(time1{ie}>=tmin2{ie} & time1{ie}<=tmax2{ie});
      if out & out2
        gca = subplot('Position',[0.05+ix*1.05*dxp 0.07+(iy*1.05)*dyp dxp dyp]);
        set(gca,'YLim',[-1 1]*tsmax,'XLim',timeLimits,'Box','on','fontsize',fs,'XDir','reverse');
      if (loop==1) % plot only few iterations
        hold on
        %%plot(tim,evinfo.value{i},'-','color',ones(3,1)*0.,'LineWidth',4);
        plot(time1{ie}(it),proxi{ie},'-','color',ones(3,1)*0.,'LineWidth',3);

        text(mean(timeLimits)*1.0,0.91*tsmax,ta{ie},'HorizontalAlignment','center','fontsize',fs,'FontWeight','b');
        text(3.9,0.91*tsmax,num2str(ie),'HorizontalAlignment','center','fontsize',fs,'FontWeight','b','Color','w','BackgroundColor','k');
        if iy>0, set(gca,'XTickLabel',[]); end
        end
      end
      i=InEur(ie);

      cs{ie}=[];csd{ie}=[]; cDTW_D{ie}=[];  cDTW_k{ie}=[];
      dists=zeros(ndbin,1);
      % ------  calculating match factor
      deltat = zeros(1,nt+1);  deltaw=deltat;
      % -------------------------------------------------------
      % loop over compared TS
      for i2=1:numprox
       if length(time1b{i2})>9
        j=InEur(i2);
        dist=cl_distance(lons(i2),lats(i2),lons(ie),lats(ie));

        ti=time1b{i2};%ti=evinfo.time{j};
        tmax=max(ti);tmax=min(tmax,tmax2{ie});
        tmin=min(ti);tmin=max(tmin,tmin2{ie});

        itj=find(time1{ie}>=tmin & time1{ie}<=tmax);
        dtw_Dist=-9E9; r=0;wp=-1;w=[]; li=length(itj);sign=1;
        if(i~=j  & li*dt0>1 ) %% & numev(i2)>0 & (dist<maxdist)
          x1=proxi{ie}(itj);
          [tiu,ia,ic]=unique(ti);
          tio=time1{ie}(itj); %time2(it(itj));
          val=proxi{i2};% evinfo.value{j};
          x2=interp1(tiu,val(ia),tio,'linear','extrap');

          [r,p1]=corrcoef(x1,x2);
          r=r(1,2);  p=p1(1,2);
          c=r*std(x2)/std(x1); %yreg=mean(x2)+c*(xreg-mean(x1));
          cs{ie}=[cs{ie} r*r];csd{ie}=[csd{ie} dist];
          % Dynamic Time Warping
          [dtw_Dist,D,k,w]=dtw(x1,x2);

          if isempty(strmatch(evinfo.Proxy{i},evinfo.Proxy{j}))
            [dtw_Dist2,D2,k2,w2]=dtw(x1,-x2);
            if(dtw_Dist2<dtw_Dist)
              dtw_Dist=dtw_Dist2; w=w2;
              sign=-1;
            end
          end

          % calculate maximal distortion in time
          for j1=1:2
            clear nn
            for i1=1:li
               nn(i1)=sum(w(:,j1)==i1);
            end
            nc(j1)=max(nn)-1;
          end

          cDTW_Dist(ie,i2)=dtw_Dist;
          % penalize short overlap and large distortion
          if isnan(dtw_Dist)
             fprintf('%d %d :%1.3f*%1.3f\n',ie,i2,dtw_Dist,nt/li*(1+sum((nc/100).^2)));
          end
          dtw_Dist=dtw_Dist*nt/li*(1+sum((nc/100).^2));

    %Dist is the unnormalized distance D is the accumulated distance k is the length of the warping path

          ti0=tio;
          tio(w(:,1))=ti0(w(:,2));
          dti=tio-ti0;
          weigh=exp(-(dtw_Dist/dtw_Dist_crit)^2);
          dtw_weigh(ie)=dtw_weigh(ie)+weigh;
          dtw_weighi(ie)=dtw_weighi(ie)+1;
          weigh=weigh*we_dtw(i2); %*we_pca(i2);

          deltat(it(itj)) = deltat(it(itj)) + dti*weigh;
          deltaw(it(itj)) = deltaw(it(itj)) + weigh;
          itj0(ie,i2)=itj(1);itj1(ie,i2)=itj(end);
          cDTW_w{ie,i2}=w;
          cDTW_sign(ie,i2)=sign;
          cPears(ie,i2)=r*r;
          cLocD(ie,i2)=wp;
        end %if i~=j  & li>9
       end
      end %for i2

      % prep PCA
      ti=time1{ie};
      if(min(ti)-tol<time20(1) & max(ti)+tol>time20(end))
         x2=interp1(ti,proxi{ie},time20,'linear','extrap');
         x2=x2-mean(x2); x2=x2/std(x2);

         datprox(iv,:) = x2;
         datproxw(iv,:)= x2*we_dtw(ie);
         dpi(iv)=ie;
         iv=iv+1;
      end
      ind=find(cDTW_Dist(ie,:)>0);
      [vs vi]=sort(cDTW_Dist(ie,ind));
      mDTW(ie)=mean(cDTW_Dist(ie,ind));
      mPear(ie)=mean(cPears(ie,ind));
      ii=find(deltaw>0);
      %if np>=0 & mod(ie,10)==0, fprintf('%d %d\t%1.0f %1.0f\tl=%d we=%1.3f\n',loop,ie,mDTW(ie),mPear(ie)*1E3,length(ii),we_dtw(ie)); end

      if ii,
        %time1=time2;%time1{ie}
        del=deltat(ii)./(1E-3+deltaw(ii));
        delm=max(abs(del));
        ff=1;
        if(tfac*delm>tol*2 & 0 )
          ff=tol*2/(tfac*delm);
        end

        delmax(ie,loop)=delm*ff;

        time1{ie}(ii)=time1{ie}(ii)+tfac*del*ff;
        %fprintf('%d del %1.3f %1.3f %1.3f\n',ie,mean(abs(del)),mean(abs(tfac*del*ff)),ff);

        dilt=time1{ie}-time0{ie};
        ii=find(dilt>tol2-dtol);
        if ii, time1{ie}(ii)=time0{ie}(ii)+tol2-dtol*exp(-(dilt(ii)-tol2+dtol)/dtol); end
        ii=find(dilt<-tol2+dtol);
        if ii, time1{ie}(ii)=time0{ie}(ii)-tol2+dtol*exp((dilt(ii)+tol2-dtol)/dtol); end
      end
      if out
        le(loop)=plot(time1{ie}(it),proxi{ie},'-','color',ppcol(1+mod(loop-1,5),:),'LineWidth',1+2*floor((loop-1)/5));
      end
      dt=find(time1{ie}(2:end)-time1{ie}(1:end-1)<0);
      while length(dt)>0
       time1{ie}(dt+1)=time1{ie}(dt)+1E-4;
       dt=find(time1{ie}(2:end)-time1{ie}(1:end-1)<0);
      end
      nl=length(dt);
      dt(dt<3)=[];dt(dt>length(time1{ie})-3)=[];
      if nl>0 & length(dt)>0, fprintf('%d %d %d\t%d %d\t%1.4f %1.4f %1.3f %1.3f\n',loop,ie0,nl,dt(1),dt(end),time1{ie}(dt(1)+(-1:2))); end

    end %for ie
    % -------------------------------------------------------
    % loop over  TS
    for ie=1:numprox
      time1b{ie}=time1{ie};
      we_dtw(ie)=dtw_weigh(ie)/(1E-3+dtw_weighi(ie));
    end
    we_dtw=we_dtw/mean(we_dtw);

    % ------ calc PCA without weights ------------------
    np=iv-1;
    pca0 = pca_model(datprox(1:np,:)',nmax,'none');%'

    % ------ calc PCA with weights ----------------------
    pca = pca_model(datproxw(1:np,:)',nmax,'none');%' 'auto'

    exp_var=pca.exp_var;
    pwei=zeros(np,1);
    for i=1:np
     for j=1:nmax
      pwei(i)=pwei(i)+exp_var(j)*abs(pca.L(i,j));
     end
    end
    datproxw=pwei.*datproxw;

    % ------ calc PCA with weights ----------------------
    pca = pca_model(datproxw(1:np,:)',nmax,'none');%' 'auto'
    exp_var=pca.exp_var;

    [m ii]=sort(-pwei);
    x=-(m-mean(m))/std(m);
    fprintf('\nstd: %1.3f\n',-std(m)/mean(m));

    for i=1:np
      we_pca(dpi(ii(i)))= 1/(1+exp(-2*x(i)));
      fprintf('%d %d:%1.2f %1.0f  ',dpi(ii(i)),InEur(dpi(ii(i))),-m(i)*100,we_pca(dpi(ii(i)))*100);
    end
    fprintf('\nie=%d\t pca 0: %1.3f %1.3f\t',ie,pca0.cum_var(2),pca0.cum_var(nmax));
    fprintf('\t pca: %1.3f %1.3f\n',pca.cum_var(2),pca.cum_var(nmax));
    fprintf('dtw: %1.3f\n',mean(mDTW));
    mlPear(loop)=mean(mPear(ie));

    % ------ calc PCA with weights ----------------------
    pca_ts=pca.T;

    % calc stability from temporal change in PCs
    kine=0;
    for pi=1:nmax
       pca_sm = movweighavg(time20*1E3,squeeze(pca_ts(:,pi)),220,40); % smoothing
       pcats{pi}=pca_sm;
       dpca = diff(pca_sm);
       kine=kine+exp_var(pi)*abs(dpca);
    end
    pca_change=-(kine-mean(kine))/std(kine);
    t_change=(time20(2)-time20(1))/2+ time20(1:end-1);
    file=sprintf('%sdtwpca/dtwpca2_%3.2f_%2.0f_%1.0f_%d.mat',outputDirectory,timeLimits(2),dtw_Dist_crit,tol*100,loop);
    save(file,'mlPear','mDTW','proxi','time1','time20','pcats','cDTW_Dist','cPears','InEur','tol','exp_var','pca_change','t_change');

  end %for loop

  if out
    %legend(le,num2str([1:loop]'))%'
    set(0,'DefaultFigureVisible','on')
    for np=1:nf-1
      set(0, 'CurrentFigure', gcf(np))
      figure(np)
      set(gcf(np),'Visible','on');
      file=sprintf('%splots/proxy_dtwp_%03.2f_%3.2f_%d.png',outputDirectory,threshold,timeLimits(2),np);
      if np==nf-1,fprintf('saving graph to %s\n',file);end
      set(gcf(np),'PaperPositionMode','auto');%,'InvertHardCopy','off'
      print(gcf(np),'-dpng','-r300', file);
    end
  end
end %tl
