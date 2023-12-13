%
% region kriging on a grid based on cluster points
%
% kai wirtz (hereon) 2023
%function gp = make_grid(offset,dll)
close all
clear all
addpath('~/tools/m_map')

%load_pars; % sets common parameters (scdir, cc, latlim, regs)
breaks=3000:400:9800;

latlim=[34 71]; lonlim=[-12 37]; % entire Europe
scdir='p3k14c/';
MaxOcc=4;
OnlySeaMask=0;
if exist('offset')
  offset=str2num(offset); npart=64;
else
  offset=1; npart=1;
end
if exist('dll')
  dlon=str2num(dll); dlat=dlon;
else
  dlon=0.05;  dlat=0.05;
%  dlon=0.2;  dlat=0.2;
end

lightgray=repmat(0.9,1,3);
%coastfilename=['europe.mat'];
coastfilename=['europe-norm.mat'];
%coastfilename='sns_coast.mat';

sc=load(coastfilename);
xcoast=sc.ncst(:,1);
ycoast=sc.ncst(:,2);
angl=(0:0.04:1)*2*3.1415;
sn=[[1 1];[1 -1];[-1 1];[-1 -1];];

%long=-5:dlon:10; latg=40:dlat:53;
long=lonlim(1):dlon:lonlim(end); latg=latlim(1):dlat:latlim(end);
nx=length(long);
ny=length(latg);

radmax = round(1.8/dlon);
for ix=0:radmax-1
  for iy=0:ix
    rad=sqrt(ix*ix+iy*iy);
    if rad<=radmax+0.01
      % weigh(1+ix,1+iy)=1./(1.+10*rad*dlon);
      % weigh(1+ix,1+iy)=exp(-(1.2*rad*dlon)^2);
       weigh(1+ix,1+iy) = exp(-(7*rad*dlon));
       weigh(1+iy,1+ix)=weigh(1+ix,1+iy);
    else
       weigh(1+ix,1+iy)=0;weigh(1+iy,1+ix)=0;
    end
  end
end
weigh(1,1)=1;
[iwx,iwy]=find(weigh>0);
%%size(xcoast)
if 0
  value=zeros(nx,ny)-9999;
  xp=ceil(nx/npart);
  for ix=1+(offset-1)*xp:min(nx,offset*xp)
    for iy=1:ny
     if inpolygon(long(ix),latg(iy),xcoast,ycoast)
       value(ix,iy)=0;
     end
    end
    fprintf('%d: %d/%d\n',offset,ix,min(nx,offset*xp));
  end
  save(['seamask_' num2str(dlon) '_' num2str(offset) '.mat'],'value','long','latg');
else
  %load seamask2_0.05.mat %seamask_High
  load seamask_norm_0.05.mat %seamask_High
  %values(1,:,:)=value;
end
%figure(12);
%imagesc(value);

if OnlySeaMask==0

%landfrac=length(find(value>=0))/length(value(:));
load([scdir 'mat/C14_p3k14c_europe_neo']);%'lonsn','latsn','C14agesn','C14SDsn','SiteIDsn','datIDsn'
nt=length(breaks);
ncol=round(sqrt(nt))+1;  nrow=ceil(nt/ncol);
dxp=0.97/ncol;dyp=0.97/nrow;

tii=0;
for ti=breaks
  if 1
    values=zeros(MaxOcc,nx,ny); %-9999 all ocean
    regs=zeros(MaxOcc,nx,ny); % all ocean
    %%[lo,la] = meshgrid(long,latg);
    %m_usercoast(coastfilename,'patch',lightgray,'LineStyle','none');
  %  load([scdir 'C14_dASIS']);
    load([scdir 'mat/clusti3_' num2str(ti) '_120']);
    lons=lonsn;lats=latsn;
    %C14ages1=C14agesn;C14SDs1=C14SDsn;SiteIDs1=SiteIDsn; datIDs1=datIDsn;
  %    clustdat=load([scdir 'cluster.dat']);
  %    regionlon=clustdat(:,4); regionlat=clustdat(:,5);
    clustn=unique(clusti);
    clustn=clustn(~isnan(clustn));
    clustn=clustn(clustn~=0);
    ncolor=length(clustn);
    fprintf('%d\n',ti);
    for i=1:ncolor
      ii=find(clusti==i);
      lon=lonsn(ii); lat=latsn(ii);
      regionlon(i)=mean(lon); regionlat(i)=mean(lat);
      fprintf('%2d:%1.1f %1.1f\t',i,regionlon(i),regionlat(i));
    end
    fprintf('\n');
    make_grid_regions
  else
    load([scdir 'mat/regiongrid3_' num2str(ti)]);
  end

  gcf=figure(2);
  set(gcf,'position',[15 20 850 750],'Color','w','Visible','on');
  clf;
  gca = subplot('Position',[0.01 0.01 0.98 0.98]);
  imagesc(flipud(reg')); %'
  for i=1:ncolor
    x=(regionlon(i)-lonlim(1))/(lonlim(2)-lonlim(1));
    y=(regionlat(i)-latlim(1))/(latlim(2)-latlim(1));
    annotation('textbox',0.01+[x y 0.05 0.05]*0.98,'string',[num2str(i) ':' num2str(round(area(i)*1E-4))],'FontWeight','bold','FontSize',24,'EdgeColor','none','Color','w');%,'FontName','Arial'
  end
  annotation('textbox',[0.15 0.82 0.12 0.05],'string',num2str(ti),'FontWeight','bold','FontSize',32,'EdgeColor','none','Color','w');%,'FontName','Arial'
  outfilename=[scdir 'plots/grid_regnum_' num2str(dlon) '_' num2str(ti) '.png'];
  print('-dpng','-r300', outfilename);

  gcf=figure(1);
  set(gcf,'position',[10 15 1450 1050],'Color','w','Visible','on');

  iy  = floor(tii/ncol);  ix  = mod(tii,ncol);
  gca = subplot('Position',[0.01+ix*1*dxp 0.01+iy*dyp dxp*0.98 dyp*0.98]);
  imagesc(flipud(val')); %'
  set(gca,'XTick',[],'YTick',[]);
  text(100,100,num2str(ti),'FontWeight','bold','FontSize',22,'Color','w');
  tii=tii+1;
  nreg(tii)=ncolor;
end %ti

% save areas to file
save([scdir 'area_' num2str(dlon) '_' num2str(breaks(1)) '_' num2str(breaks(end)) '.mat'],'areav','nreg');

% save plot to file
set(gcf,'PaperPositionMode','auto','InvertHardCopy','off','Visible','on');
outfilename=[scdir 'plots/grid' num2str(dlon) '.png'];
print('-dpng','-r300', outfilename);

end
