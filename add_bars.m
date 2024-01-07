% add background bars

% ------- RGR phases
xx=tip2(isp);
xx=reshape(xx,1,length(xx));

% ------- left edge
ii2=(i0:i01)-(i0-2);
nli=length(ii2);
x=[xx(ii2) fliplr(xx(ii2))];
y=[yli(iy1)*ones(1,2) yli(iy1+1)*ones(1,2)];
ii20=ii2(end);
for i3=1:nli-1
  col=bgcol*grv; fac=(i3-1)/(nli-1);
  patch('XData',[xx(ii2(i3:(i3+1))) fliplr(xx(ii2(i3:(i3+1))))],'YData',y,'Facecolor',col,'EdgeColor','none','FaceAlpha',fa*fac);
end

% ------- right edge
ii2=(i10:i1)-(i1-nl+1);
ii2(ii2<ii20)=[];
ii2(ii2<=0)=[];
nli=length(ii2);
for i3=1:nli-1
  col=bgcol*grv; fac=(nli-1-i3)/(nli-1);
  patch('XData',[xx(ii2(i3:(i3+1))) fliplr(xx(ii2(i3:(i3+1))))],'YData',y,'Facecolor',col,'EdgeColor','none','FaceAlpha',fa*fac);
end

% ------- central bar
if(i10>i01+1)
   ii2=((i01+1):(i10-1))-(i0-1);
   ii2(ii2<=0)=[];nli=length(ii2);
   patch([xx(ii2) fliplr(xx(ii2))],[yli(iy1)*ones(1,nli) yli(iy1+1)*ones(1,nli)],'-','Facecolor',bgcol*grv,'EdgeColor','none','FaceAlpha',fa);
end

% ------- overlap of RGR phases with 2nd variable
isp1=find(tip2>=min(tip2(isp)) & tip2<=max(tip2(isp)));
ts=reshape(ts,1,length(ts));
tsi=ts(isp1);

ix=find(sign(tsi)~=sgp);
tsi(ix)=0;
xx=tip2(isp1); %
xx=reshape(xx,1,length(xx));
nl1=length(isp1);

iic=4:(length(xx)-3);
nl1=length(iic);
patch([xx(iic) fliplr(xx(iic))],[tsi(iic) zeros(1,nl1)],'-','Facecolor',pcol,'EdgeColor','none','FaceAlpha',0.1);
