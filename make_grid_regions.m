%
% kriging: distribution of site positions -> grid occupation
%
% kai wirtz (hereon) Dec 2023
%
vc=1; minc=100;

% loop over clusters/regions
for i=1:ncolor
  if exist('clustdat')
     ii=find(clusti==clustn(clustdat(i,2)));
  else
     ii=find(clusti==i);
  end
  ii(find(ii>length(lats)))=[];

% loop over sites in cluster
  for ip=1:length(ii)

  % check for singularity points
  % after clustering, few sites can be located within or close to 'neighbor' clusters

  % distance to nearest neighbor site of the same cluster should be minc (=100km)
    mind=9E9;
    jp=1;
    % only screen  north and east regions for singularity points
    if(regionlon(i)>30 | regionlat(i)>59)
      while mind>minc & jp<=length(ii)
        cd=cl_distance(lons(ii(ip)),lats(ii(ip)),lons(ii(jp)),lats(ii(jp)));
        if cd<mind & cd>2 & ip~=jp
          mind=cd;
        end;
        jp=jp+1;
      end
    end

    if jp<=length(ii)

    % calc grid index position from site geo-location
      ix= 1+floor( (lons(ii(ip))-long(1))/dlon);
      iy= 1+floor( (lats(ii(ip))-latg(1))/dlat);

    % valid grid index ?
      if ix<=nx & iy<=ny & ix>0 & iy>0
        iso=0;

        % occupation of grid cell by other clusters
        for nv=1:MaxOcc
           iso=iso+(regs(nv,ix,iy)~=i | values(nv,ix,iy)<vc);
        end

        % full occupation by others or small contribution by own cluster:
        if iso==4 & value(ix,iy)>-1E-3
         for dd=1:4
          sea=0;occ=0;

          %check for marine straits
          for jxi=1+(sn(dd,1)<0):radmax*0.7
            jx=jxi-1;
            for jyi=1+(sn(dd,2)<0):radmax*0.7
              jy=jyi-1;
              x1=ix+jx*sn(dd,1);y1=iy+jy*sn(dd,2);
              if  x1<=nx &  x1>0 & y1<=ny &  y1>0
               if value(x1,y1)<0
                  sea=sea+1./(jx*jx+jy*jy);
               end %if

              end %if
            end % for jyi
          end % for jxi

          % loop over neighboring grid cells - x
          for jxi=1+(sn(dd,1)<0):radmax %*1/(0*sea+1)
            jx=jxi-1;
            % loop over neighboring grid cells - y
            for jyi=1+(sn(dd,2)<0):radmax
              jy=jyi-1;
              if weigh(jxi,jyi)/(sea*sea*sea+1)>5E-4 % inside radius
                 x1=ix+jx*sn(dd,1);
                 y1=iy+jy*sn(dd,2);
                 if  x1<=nx &  x1>0 & y1<=ny &  y1>0  % check for domaiin
                  %   fprintf('%d\t %1.3f %1.3f\t%d  %d %d\t %1.4f ->',ip,lons(ii(ip)),lats(ii(ip)), dd,ix+jx,iy+jy,values(1,ix+jx,iy+jy))
                    if(value(x1,y1)>=0)
                      iso=1;nv=1;
                      while iso==1 & nv<=MaxOcc
                      % occupied cell by 1st
                        iso=(regs(nv,x1,y1)~=i & regs(nv,x1,y1)>0);
                        nv=nv+1;
                      end
                      if iso==0 & nv==MaxOcc+1 % all 4 layers occupied: find lowest layer
                        [mv mi]=min(values(:,x1,y1));
                        if mv<weigh(jxi,jyi)
                          regs(mi,x1,y1)=i;
                          values(mi,x1,y1)=weigh(jxi,jyi);
                        end
                      else % set a new entry
                        nv=nv-1;
                        values(nv,x1,y1)=values(nv,x1,y1)*(regs(nv,x1,y1)==i)+weigh(jxi,jyi);
                        regs(nv,x1,y1)=i;
                      end
                     %fprintf('%1.4f\n',values(1,ix+jx,iy+jy))
                   end %if(value(x1,y1)>-1)
                  end % if in domain
                end % for dd
             %fprintf('%d %d\t%1.3f\n',ix+jx,iy+jy,values(1,ix+jx,iy+jy))
          end % if weigh
        end % for jyi
       end % for jxi
      end % new
    %  fprintf('%d %d\t%1.3f\n',ix+jx,iy+jy,values(1,ix+jx,iy+jy))
     end
    else
     fprintf('rm %d %d\t %1.1f/%1.1f\t%1.1f/%1.1f\t%1.3f\n',i,ip,lons(ii(ip)),regionlon(i),lats(ii(ip)),regionlat(i),mind)
    end %if mind
  end % for ip
end %i

for i=0:-ncolor
    fprintf('%d %d %d\n',i,length(find(regs(1,:,:)==i)),length(find(regs(2,:,:)==i)));
end

ind1=find(value>=0);
vmax=4;
for nv=1:MaxOcc
  ind=find(values(nv,:,:)>vmax);
  values(nv,ind)=vmax;
end
%fprintf('MEAN %1.3f\tMAX %1.3f\n',mean(mean(values(1,ind1))),max(max(values(1,ind1))));

% normalize all value fields
sval=(sum(values,1))+1E-5;%squeeze
%fprintf('sval MEAN %1.3f\tMAX %1.3f\n',mean(mean(sval(ind1))),max(max(sval(ind1))));
ind=find(value<0);
for nv=1:MaxOcc
  values(nv,:,:)=values(nv,:,:)./sval;
  % mask sea areas
  values(nv,ind)=NaN; regs(nv,ind)=0;
end
fprintf('MEAN %1.3f\tMAX %1.3f\n',nanmean(nanmean(values(1,ind1))),max(max(values(1,ind1))));

val=squeeze(values(1,:,:));
reg=squeeze(regs(1,:,:));
ind=find(val<0);val(ind)=NaN;reg(ind)=0;

% sort: bring higher weight to front
for nv=2:MaxOcc
  ind=find(squeeze(values(nv,:,:))>val);
  if ind, reg(ind)=squeeze(regs(nv,ind)); val(ind)=values(nv,ind); end
end

%calc area of each region
clear area
for i=1:ncolor
  x=regionlon(i);
  y=regionlat(i);
  arf=cl_distance(x,y,x,y+dlat);
  arf=arf*cl_distance(x,y,x+dlon,y);
  ind=find(reg==i);
  area(i)=arf*length(ind);
  fprintf('%2d %1.0f\t%1.2f * %d\n',i,area(i),arf,length(ind));
 end
areav{tii+1}=area;

% store region information in binary file
save([scdir 'mat/regiongrid_' num2str(ti)],'val','reg','values','regs','ncolor','regionlon','regionlat','area');
