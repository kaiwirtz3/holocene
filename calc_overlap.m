%
% ------- calculate overlap

% coefficients for linear filter
tmaxo=8.15;
fi=0;tot=0;tweigh=0;
ii2=find(tip2<=tmaxo & tip2>=3);
crit=0.5*std(avgrde(ii2));
%%fprintf('crit=%1.1f\n',crit*1E3);
dtt=abs(tip2(2)-tip2(1));
delta=0.095;
grv  =1;
dim=round(delta/dtt);
%dim=3;crit=0.5*stdc;

% graphical settings
yl=get(gca,'Ylim');
eps=abs(yl(2)-yl(1))/700;
yli=[yl(2)-eps eps -eps yl(1)+eps];
fa =0.667; % 'FaceAlpha' factor 1

pcol=[0.6 0.2 0.];      % bar colour
bgcol=[225 200 150]/230;%

% loop over booms and busts
for i=1:2
  overlap =ones(1,size(dat,1));
  boombust=zeros(1,size(dat,1));
  sgp=(3-2*i); % sign of booms and busts
  iy1=2*(i-1)+1;

  % find positive or negative phases/peaks
  [pg pgi]=findpeaks(tip2,avgrde*sgp,stdc,9); %0.667

  % delete peaks outside time-window
  ii2=find(pg<3 | pg>tmaxo);
  %fprintf('delete %d peaks:\t',length(ii));
  if ii2, pg(ii2)=[]; pgi(ii2)=[]; end

  jj=length(pg);
  %%fprintf('%d\tpeaks:%d\n',i,jj);

% ------- loop over phases/peaks
  for j3=1:jj
    xm=pg(j3);  % time
    im=pgi(j3); % index
    weigh=zeros(1,size(dat,1)); % clear weighing vector

    % end of phase (smallest age kaBP)
    i0 = max(find(tip2<xm & avgrde*sgp < crit));
    di = round(min(1,abs(avgrde(im))/crit)*dim);
    if isempty(i0) i0=im; end
    if i0<=di, i0=1;
    else i0 = i0 - di;
    end

    % start of phase
    i1=min(find(tip2>xm & avgrde*sgp < crit));
    if isempty(i1), i1=im; end
    if i1>=length(tip2)-di, i1=length(tip2);;
    else  i1 = i1 + di;
    end

    % linear weighing change at edge
    weigh(i0:(im-1))=((i0:(im-1))-i0)/dim;
    weigh(im:i1)=(i1-(im:i1))/dim;
    i01=min(find(weigh>=1)); if isempty(i01) i01=im; end
    i10=max(find(weigh>=1)); if isempty(i10) i10=im; end
    weigh(weigh>1)=1;

    tweigh=tweigh+sum(weigh);
    isp=i0:i1; isp=isp'; %'  vector of indices
    nl=length(isp);

    boombust(min((isp)):max((isp)))=1; % mark phases
    %%fprintf('%d %d\t%1.1f\t%d - %d\t%1.1f %1.1f - %1.1f %1.1f\n',j3,sgp,xm,min(isp),max(isp),avgrde(i0)*1E3,avgrde(i0+1)*1E3,avgrde(i1-1)*1E3,avgrde(i1)*1E3);

    % add background bars
    add_bars
    % overlap : weighed sum of product
    df = nansum(abs(tsi.*weigh(isp1).*overlap(isp1)));
    fi = fi+df;
    dft= nansum(abs(tsi.*overlap(isp1))); %isp
    tot= tot+dft;

    overlap(isp)=0; % avoid double counting

  end %for j3

  %% dat(:,9+i) = boombust;
  %% pgv{i}=pg;
end % for i=1:2

ii2 = find(tip1<tmaxo & tip1>=3);
to2 = nansum(abs(ts(ii2)));
fr  = fi/(to2*0.92); % match factor of equal periodic signals

%fprintf('fi=%1.3f\tsum %1.2f %1.2f\tfr=%1.3f\n',fi*1E3,tot*1E3,to2*1E3,fr*100);
