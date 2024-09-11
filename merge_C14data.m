% This script merges C14 databases, creates a grid, and plots sites.
% It processes data from two sources (dASIS+ and p3k14c), merges them,
% and generates plots for different continental domains.
% The script also saves the processed data and plots to files.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz  <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input files:
% - c14mat/C14_dASIS.mat
% - c14mat/p3k14c.mat

% Output files:
% - c14mat/C14_<continent>_<mode>_<patch>.mat
% - <outputDirectory>/plots/map_<continent><mode>.png

% Variables:
% - lats, lons: Latitude and longitude of sites
% - C14ages, C14SDs: C14 ages and their standard deviations
% - SiteIDs, datIDs: Site and data identifiers
% - clust: Array defining the boundaries of continental domains
% - nx, ny: Grid resolution parameters
% - mode: Configuration mode for data processing
% - gcf, gca: Figure and axis handles for plotting
% - cc: Continent name
% - coastfilename: Filename for high-resolution coastline data
% - pclust: Vector for all dates and sub-continental patches
% - bo: Boundary coordinates for the current cluster
% - ii, ii0, ii1: Indices for site selection and processing
% - lon, lat, C14age, C14SD, SiteID, datID: Processed data for saving

close all; %clear all
addpath('~/tools/m_map'); %
load_pars; % sets common parameters (outputDirectory, cc, latitudeLimits, regs)

% ----------  division lines for continental domain
% left, bottom, right, top
cont='';clust=[[-12 33 37 71];];contname={'europe'};
%%cont='EA_';clust=[[69 20 149 58];]; contname={'EAsia'};%
%cont='IreN_';clust=[[-10 54 -5.5 55.3];];
%cont='NA_';clust=[[-125 25 -60 51];];contname={'NAmerica'};
%cont='SA_';clust=[[-82 -56 -34 12];];contname={'SAmerica'};
%cont='Af_';clust=[[-19 -36 42 36];];contname={'Africa'};
%cont='Au_';clust=[[112 -40 155 -9];];contname={'Australia'};

% load own compilation (dASIS+)
load(['c14mat/C14_dASIS']);
lats0=lats; lons0=lons; C14ages0=C14age; C14SDs0=C14SDs; SiteIDs0=SiteIDs; datIDs0=datIDs;

% load p3k14c dates
load(['c14mat/p3k14c']);

% grid resolution
nx=8; ny=8;

% ----------  loop over different data-base configurations
for mode=1:1+3*strcmp(contname,'europe') % modes 2-4 only for Europe

% open figure
  gcf=figure(mode); clf;
  set(gcf,'position',[1+mode*20 1 880 900],'Color','w','Visible','on');
  gca=subplot('Position',[0.0 0.01 0.99 0.98]);
  hold on;
  dlo=0.; %add margin: dlo=1;
  m_proj('mercator','lat',[clust(2)+dlo clust(4)-dlo],'long',[clust(1)+dlo clust(3)-dlo]); %albers equal-area ,'rect','on'mercator
  m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'linest','none','fontsize',18,'backcolor',[0.8 0.8 0.9]);%'box','fancy',
  %% Coast line
  cc = cell2mat(contname);
  coastfilename=['data/' cc '-high.mat'];
  % if no high resolution coastline is found, create one; this takes time ...
  if ~exist(coastfilename,'file')
      fprintf('creating high resolution coastline (needs patience)\n')
      m_gshhs('hc','save',coastfilename); %'hc' high resolution coastline
      fprintf('... ready \n')
  end
  lightgray=repmat(1,3,1);
  m_usercoast(coastfilename,'patch',lightgray','LineStyle','none'); %'

  %%ncolor=length(cc)-1; % number of sub-continental patches
  ncolor=1; % number of sub-continental patches
  cmap=jet(nx*ny+1);

  % create vectors for all dates and sub-continental patches
  pclust = zeros(length(C14age),1);

  bo=clust(1,:);
  switch(mode) % data-base configurations
    case 1  % only p3k14c
      lons1=lons;lats1=lats;C14ages1=C14ages;C14SDs1=C14SDs;SiteIDs1=SiteIDs; datIDs1=datIDs;

    case 2  % only dASIS
      lons1=lons0;lats1=lats0;C14ages1=C14ages0;C14SDs1=C14SDs0;SiteIDs1=SiteIDs0; datIDs1=datIDs0;

    case 3  % merged p3k14c+dASIS

      % identify individual sites according to geocoordinates at 3 digit accuracy

      i2=find( lons>=bo(1) & lats>=bo(2) &  lons<bo(3) & lats<bo(4) );
      a=num2str(lats(i2),'%05.3f');
      b=num2str(lons(i2),'%05.3f');
      ab=[a b]; site1=unique(cellstr(ab));

      i0=find( lons0>=bo(1) & lats0>=bo(2) &  lons0<bo(3) & lats0<bo(4) );
      a=num2str(lats0(i0),'%05.3f');
      b=num2str(lons0(i0),'%05.3f');
      ab0=[a b]; site0=unique(cellstr(ab0));

      % find identical sites
      [site1d,ia] = setdiff(site1,site0);
      ii=[];
      for i=1:length(site1d)
        in=strmatch(site1d{i},ab,'exact');
        ii=[ii' in']'; %'
        if mod(i,500)==0
          fprintf('%5d/%5d\n',i,length(site1d))
        end
      end

      % select only new entries
      ii=i2(ii);
      lons1=lons(ii);lats1=lats(ii);C14ages1=C14ages(ii);C14SDs1=C14SDs(ii);SiteIDs1=SiteIDs(ii); datIDs1=datIDs(ii);

      % add new p3k14c to dASIS
      lonsn=[lons0(i0)' lons(ii)']'; latsn=[lats0(i0)' lats(ii)']'; C14agesn=[C14ages0(i0)' C14ages(ii)']';
      C14SDsn=[C14SDs0(i0)' C14SDs(ii)']';  datIDsn=[datIDs0(i0) datIDs(ii)];
      SiteIDsn=1:length(lonsn);
      fname=['c14mat/C14_' cc];
      fprintf('saving %02d dates into %s\n',length(lonsn),fname);
      save(fname,'lonsn','latsn','C14agesn','C14SDsn','SiteIDsn','datIDsn');

    case 4 % ------ loop over patches

      lons1=lonsn;lats1=latsn;C14ages1=C14agesn;C14SDs1=C14SDsn;SiteIDs1=SiteIDsn; datIDs1=datIDsn;

      % adjust European grid coordinates to topography
      yoff=zeros(1,nx)+1.5;
      xoff=1.;
  %%    yoff(1:round(nx/2.5))=1;
      yoff(1)=3; yoff(2)=3.5; yoff(3:4)=4.5;
      cl0=clust(1,:);
      ncolor=nx*ny;
      dx=(cl0(3)-cl0(1)-2.36)/(nx);
      dy=(cl0(4)-cl0(2)-1.2)/(ny);
      i=1;
      % create boundary boxes
      for ix=0:nx-1
        for iy=0:ny-1
          clust(i,:)=[xoff+cl0(1)+ix*dx  cl0(2)+iy*dy+yoff(1+ix) xoff+cl0(1)+(ix+1)*dx  cl0(2)+(iy+1)*dy+yoff(1+ix)];
          i=i+1;
        end
      end
    end % switch

    ii0=[];ii1=[];
    nc=400;
    % loop over domains: 1 for pooled continent
    for i=1:ncolor %length(cl)
      n=0;
      bo=clust(i,:);
      ii=find( lons1>=bo(1) & lats1>=bo(2) &  lons1<bo(3) & lats1<bo(4) );
      fprintf('%02d: %04d\t+ %d  %d? +%d ->\n',i,length(ii),length(ii0),mod(i,nx)==0,length(ii1));
      n=length(ii);
      if mod(i,ny)==0
         ii=[ii' ii0' ii1']';ii1=[];ii0=[];
         if (length(ii)<nc & i<ncolor)
           ii1=ii; ii=[];
         end
      else
         ii=[ii' ii0']';  %'
         ii0=[];
         if (length(ii)<nc)
           ii0=ii; ii=[];
         end
      end
      if i==ncolor, ii=[ii' ii0']';%'
      end
      fprintf('\t%04d (%d %d)\n',length(ii),length(ii0),length(ii1));

      lon=lons1(ii);lat=lats1(ii);C14age=C14ages1(ii);C14SD=C14SDs1(ii);SiteID=SiteIDs1(ii); datID=datIDs1(ii);
      %correct_ID
      file=['c14mat/C14_' cc '_' num2str(mode) '_' num2str(i)];
      if mode==4,
        file=['c14mat/C14_' cc '_' num2str(i)];
      else
        file=['c14mat/C14_' cc '-' num2str(mode)];
      end

      save(file,'lon','lat','C14age','C14SD','SiteID','datID','bo');
      %fprintf('saving C14_%s_%d_%d\n',cc,mode,i);
      m_plot(lon,lat,'o','Markersize',4,'MarkerFaceColor',cmap(i,:),'Color','none');
      m_text(mean(lon)-1.5,mean(lat),[num2str(i) ':' num2str(length(C14age),'%d')],'Fontsize',20,'Color','k','Fontweight','bold');
    end %for i

  % save plot to file
  set(gcf,'PaperPositionMode','auto', 'InvertHardCopy','off','Visible','on');
  outfilename=[outputDirectory 'plots/map_' cont num2str(mode) '.png'];
  print('-dpng','-r300', outfilename); % '-r600' for reproduction quality
end % mode
