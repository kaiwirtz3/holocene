%
% load logit model results and
%  store low-high-pass filtered probability difference to matrix
%
% kai wirtz (hereon) Dec 2023
%
load([scdir 'target_ts_0.mat']); %'dat','legdat'
nd0 = length(legdat); %size of old data matrix

% ------ load logit model results glmres_-TSI_clim_area.mat
vis={'-TSI_clim_area','-TSI_clim','area_area','3'};
varname={'TSI','stress-free tree growth','climate stability','RGR($t$-350a)'};
for nv=1:3
  % name of model variant
  tag=[vis{nv}]; %'8.8_'
  load([scdir 'glmres_' tag '.mat']); %2
  %skill=val,timres,prob,tprob,statcoeff,PearsTot,bestbb=bobu[bi],bestcor=corn[bc]

  sp=sum(prob,2); fac=sp(2)/sum(sp);
  tip1=dat(:,1)'; %' time vector
  it=find(tip1>=timres(1) & tip1<=timres(end)); %it=find(time>=timres(1) & time<=timres(end));

  % new index; extend data matrix
  j = nd0+nv;

  % loop over boom/bust/pseudo-RGR
  for i=3:3 % only pseudo-RGR
    if i<=2
     ts=prob(i,:);  % probability of boom or bust
    else
     % probability of boom minus probability of bust
     ts=fac*prob(1,:)-(1-fac)*prob(2,:);
    end
    ts=interp1(timres,ts,tip1(it),'linear','extrap');
    ts= movweighavg(tip1(it)*1E3,ts,tmov,toff); %
% refine description
    if i==3,
       dat(it,j)=ts;
       switch(nv)
       case 1
           adds='(climate+)';
       case 2
           adds='(climate)';
       case 3
           adds='(autoregressive lags=350a+210a)';
       end
       legdat{j}=['logit model ' num2str(4-nv) 'V ' adds];
    end %if i==3
  end %for i
end %for nv

% write out variables and indices
fid=fopen('legdat.dat','w');
j=floor(length(legdat)/2);
for i=1:j
  fprintf(fid,'%2d %25s\t\t',i,legdat{i});  fprintf(fid,'%2d %23s\n',j+i,legdat{j+i});
  fprintf('%2d %25s\t\t',i,legdat{i});  fprintf('%2d %23s\n',j+i,legdat{j+i});
end
fclose(fid);
