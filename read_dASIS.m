close all;
load_pars; % sets common parameters (scdir, cc, latlim, regs)

% load own compilation (dASIS+)
fname='C14_dASIS'; off=3; IDi=2; Dbi=8; lli=4;
% load p3k14c dates
%%fname='p3k14c';  off=0; IDi=1; Dbi=17; lli=7;

% read C14 dates from XLS file
  [num,txt,raw] = xlsread(['c14mat/' fname '.xlsx']);  %"labnr"	"c14age"	"c14std" "lat"	"lon"

  % distribute data into named vectors
  C14age=num(:,1+off);
  C14SD=num(:,2+off);
  SiteID=txt(2:end,IDi);
  Dbase=txt(2:end,Dbi);
  lon=num(:,lli+off); lat=num(:,lli+1+off);

  % correct few wrong lon/lat value formats
  ii=find(abs(lon)>70);
  if ii, lon(ii)=lon(ii)*1E-3; end
  ii=find(abs(lat)>90);
  if ii, lat(ii)=lat(ii)*1E-3; end

  % identify and remove NaNs
  ii=find(~isnan(C14age) & ~isnan(C14SD) & ~isnan(lon) & ~isnan(lat) & ~contains(Dbase,'OUT'));
  fprintf('%s: %d -> %d (no NaNs, OUT)\n',fname,length(lat),size(ii));
  lat=lat(ii);lon=lon(ii);C14age=C14age(ii);C14SD=C14SD(ii); SiteID=SiteID(ii);

  % uses index vector instead of SiteID strings (because of UTF-8 problems)
  SiteID=strtrim(SiteID);
  [situ,ia,ic]=unique(SiteID);
  SiteID=num2str(ic);
  fprintf('from %d sites\n',length(ia));

  lats=lat; lons=lon; C14ages=C14age; C14SDs=C14SD; SiteIDs=SiteID;
  datIDs=1:length(C14age); % create index vector

  % save data to matlab binary
  save(['c14mat/' fname],'C14age','lons','lats','C14SDs','SiteIDs','datIDs');
