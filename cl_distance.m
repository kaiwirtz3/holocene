function distance=cl_distance(lon,lat,lon0,lat0)
% distance=cl_distance(lon,lat,lon0,lat0)

% requires m_map library

if nargin<2 error('At least two input arguments are required'); end
if ~exist('lon0','var') lon0=lon; end
if ~exist('lat0','var') lat0=lat; end

if length(size(lon))~=2 error('lon/lat must have at most 2 dimensions'); end
if length(size(lat))~=2 error('lon/lat must have at most 2 dimensions'); end

if numel(lon0)==1 && numel(lat0)==1
  lon0=repmat(lon0,size(lon,1),size(lon,2));
  lat0=repmat(lat0,size(lat,1),size(lat,2));
else
  error('Multi-element reference lon/lat not implemented yet');
end

n=numel(lat);
lat=reshape(lat,n,1);
lat0=reshape(lat0,n,1);
lon=reshape(lon,n,1);
lon0=reshape(lon0,n,1);

dlat=[lat lat0];
dlon=[lon lon0];

dlat=reshape(dlat',2*n,1);
dlon=reshape(dlon',2*n,1);
dists=m_lldist(dlon,dlat);
dists=dists(1:2:end);

if nargout>0 distance=dists; end

return;
end
