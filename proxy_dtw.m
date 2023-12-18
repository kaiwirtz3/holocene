#
# calculates Dynamic Time Warping  (DTW) of paleoclimate proxy time-series
#
# kai wirtz (HZG) 2023
#
clear all;
addpath('~/tools/mlabtools/pca_toolbox_1.5')

% settings
rinit='0';
dtw_Dist_crit=90;  %% 50
tol=0.1;
fprintf('dtw_Dist_crit=%1.0f tol=%1.3f\n',dtw_Dist_crit,tol)
out=0;
nmax=5;  tfac=0.5;
ndbin=50; %showcor=0;showextra=0; res=0.3;
tol2=2*tol; dtol=tol/2;
ppcol=prism(7);ppcol(3,:)=[];
fs=20; lw=1;
% time segments
timelimv=[[2.8 6.2];[5.8 9.5]];
tsmax=3.;
ncol=4;  nrow=4; % plot display: number of rows and columns
dt0=0.01; tfac=0.2;

%load_pars;
scdir='out/';
edir='paleoclim/';

% loop over time segments
for tl=1:2
  timelim=timelimv(tl,:);
  % time vectors
  time2=(timelim(1)-tol):dt0:(timelim(2)+tol);
  time20=(timelim(1)):dt0:(timelim(2));
  nt=length(time2)-1;
  dxp=0.9/ncol;dyp=0.89/nrow; fs=20;
  mcol(1,:)=[1 0.5 0]; mcol(2,:)=[0 0.5 1];

  % load paleoclimate proxy time-series and meta-info
  file=sprintf('%sInEur_%d_%03.0f.mat',edir,27,920);
  load(file);
  evinfofile=sprintf('%sproxydescription_488_0_%03.2f_%03.2f.mat',edir,11.,1.7);
  load(evinfofile); % into struct evinfo

  tmin=evinfo.t_min(InEur); tmax=evinfo.t_max(InEur);
  iin=find(tmin<timelim(1)+tol & tmax>timelim(2)-tol);
  InEur2=InEur;
  numprox=length(InEur);

  % geographical position
  lons=evinfo.Longitude(InEur); lats=evinfo.Latitude(InEur);vli=length(InEur);

  % -------------------------------------------------------
  ei=1;
  wps=zeros(numprox,1);match=wps;matchInd=zeros(1+2,numprox);;%c50=wps;
  cDTW_Dist=zeros(numprox,numprox);
  cDTW_sign=cDTW_Dist; cPears=cDTW_Dist;
  cLocD=cDTW_Dist;
  iv=1;clear datprox;

  % loop over proxy TS
  for ie=1:numprox
    % -------------------------------------------------------
    i=InEur(ie);
    tim=evinfo.time{i};
    scs=evinfo.Proxy{i};
    scs(regexp(scs,'[$]'))=[];%%scs(regexp(scs,'\'))=[]; %'
    ta{ie}=sprintf('%s (%s)',evinfo.Plotname{i},scs);
    ti=tim;
    val=evinfo.value{i};

    if(max(ti)>time20(end)-tol2 & max(ti)<time20(end))
      fprintf('%d %1.3f %1.3f\t',ie,min(ti),max(ti))
      ti=reshape(ti,1,length(ti));
      val=reshape(val,1,length(val));
      if(ti(2)<ti(1))
        ti=[time20(end) ti]; val=[val(end) val];
      else
        ti=[ti time20(end)]; val=[val val(end)];
      end
      fprintf('%1.3f %1.3f\n',min(ti),max(ti))
    end

    tmax2{ie}=max(ti)+1*tol;tmin2{ie}=min(ti)-1*tol;
    it=find(time2>=tmin2{ie} & time2<=tmax2{ie});
    [tiu,ia,ic]=unique(ti);
    time1{ie}=time2(it);

    proxi{ie}=interp1(tiu,val(ia),time1{ie},'linear','extrap');
    %time1{ie}=ti;
    time0{ie}=time1{ie};
    time1b{ie}=time1{ie};

    if(min(ti)-tol<time20(1) & max(ti)+tol>time20(end))
       x2=interp1(time1{ie},proxi{ie},time20,'linear','extrap');
       x2=x2-mean(x2); x2=x2/std(x2);
       datprox(iv,:)=x2;
       dpi(iv)=ie;
       iv=iv+1;
    else
       fprintf('skip %d\t%1.3f %1.3f\n',ie,min(ti),max(ti));
    end
  end

  % calc PCA
  np=iv-1;
  pca = pca_model(datprox(1:np,:)',nmax,'auto');%'
  fprintf('\npca: %1.3f %1.3f\n',pca.cum_var(2),pca.cum_var(nmax));

  we_pca=ones(numprox,1); we_dtw=ones(numprox,1);
  for loop=1:8
      out2=(loop==-1)+(loop==-2);

    ie0=0;nf=1;iv=1;
    clear datprox;
    clear datproxw;
    dtw_weigh=zeros(numprox,1); dtw_weighi=dtw_weigh;

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
        set(gca,'YLim',[-1 1]*tsmax,'XLim',timelim,'Box','on','fontsize',fs,'XDir','reverse');
      if (loop==1)
        hold on
        %%plot(tim,evinfo.value{i},'-','color',ones(3,1)*0.,'LineWidth',4);
        plot(time1{ie}(it),proxi{ie},'-','color',ones(3,1)*0.,'LineWidth',3);

        text(mean(timelim)*1.0,0.91*tsmax,ta{ie},'HorizontalAlignment','center','fontsize',fs,'FontWeight','b');
        text(3.9,0.91*tsmax,num2str(ie),'HorizontalAlignment','center','fontsize',fs,'FontWeight','b','Color','w','BackgroundColor','k');
        if iy>0, set(gca,'XTickLabel',[]); end
        end
      end
      i=InEur(ie);

      cs{ie}=[];csd{ie}=[]; cDTW_D{ie}=[];  cDTW_k{ie}=[];
      dists=zeros(ndbin,1);
      % ------  calculating match factor
      deltat = zeros(1,nt+1); deltaw=deltat;
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
      if np>=0 & mod(ie,10)==0, fprintf('%d %d\t%1.0f %1.0f\tl=%d we=%1.3f\n',loop,ie,mDTW(ie),mPear(ie)*1E3,length(ii),we_dtw(ie)); end

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
      i1=0;
      for ij=[] %length(ind)
        i2=ind(vi(ij));
        j=InEur(i2);
        ii=itj0(ie,i2):itj1(ie,i2);
        ti=time2(it(ii));ti0=ti;
        tia=evinfo.time{j};
        val=evinfo.value{j}*cDTW_sign(ie,i2);
        x2=interp1(tia,val,ti,'linear','extrap');
        ta=sprintf('%s(%d %d) %1.0f',evinfo.Plotname{j},j,cDTW_sign(ie,i2),cDTW_Dist(ie,i2));
        if out & out2
          col=ppcol(1+i1,:);
          text(5.6,-2+i1*0.8,ta,'fontsize',fs,'color',col,'FontWeight','b');
          plot(tia,val,'-','color',col,'LineWidth',2);

          t2a=ti0+0.5*dti;      t2b=ti0-0.5*dti;
          plot(ti,dti*2,'-','color','c','LineWidth',2);
          plot(t2a,proxi{ie}(ii),'-','color',ones(3,1)*0.5,'LineWidth',3-2*i1);
          plot(t2b,x2,'-','color',[0.8 0.4 0.2],'LineWidth',3-2*i1);
        end
        i1=i1+1;
      end %for ij
      %%text(timelim(2)+1.3-loop*0.8,0.72*tsmax,[num2str(mDTW(ie),'%3.f') ' ' num2str(mPear(ie)*1E3,'%3.0f')],'fontsize',fs-5,'FontWeight','b');
    end %for ie

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
    kine=0;
    for pi=1:nmax
      %pca = spline(time20,squeeze(pca_ts(:,pi)),tip2(itp));
       pca_sm = movweighavg(time20*1E3,squeeze(pca_ts(:,pi)),220,40);%130,20
       pcats{pi}=pca_sm;
       dpca = diff(pca_sm);
       kine=kine+exp_var(pi)*abs(dpca);
    end
    pca_change=-(kine-mean(kine))/std(kine);
    t_change=(time20(2)-time20(1))/2+ time20(1:end-1);
    file=sprintf('%sdtwpca/dtwpca2_%3.2f_%2.0f_%1.0f_%d.mat',scdir,timelim(2),dtw_Dist_crit,tol*100,loop);
    save(file,'mlPear','mDTW','proxi','time1','time20','pcats','cDTW_Dist','cPears','InEur2','tol','exp_var','pca_change','t_change');
    %%mlPear
    % tfac=tfac*0.8;
     %tfac=min(tfac,0.9);
  end %for loop
  if out
    %legend(le,num2str([1:loop]'))%'
    set(0,'DefaultFigureVisible','on')
    for np=1:nf-1
      set(0, 'CurrentFigure', gcf(np))
      figure(np)
      set(gcf(np),'Visible','on');
      file=sprintf('%splots/proxy_dtwp2_%03.2f_%3.2f_%d.png',scdir,threshold,timelim(2),np);
      if np==nf-1,fprintf('saving graph to %s\n',file);end
      set(gcf(np),'PaperPositionMode','auto');%,'InvertHardCopy','off'
      print(gcf(np),'-dpng','-r300', file);
    end
  end
end %tl
