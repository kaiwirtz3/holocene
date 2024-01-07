function [ptimes,index]=findpeaks(tim,val,pc,tm)
if ~exist('tm') tm=7.6; end

% check for correct value and sign
ii=find( (abs(val)>(pc)) & sign(val)==1 );
i=1;
ptimes=[];index=[];
if ii,
  di=ii(2:end)-ii(1:end-1);
  inext=find(di>5);% distant
  i0=1;
  for i=1:length(inext)
    t0=tim(ii(i0));
    t1=tim(ii(inext(i)));
    ini=ii(i0):ii(inext(i));
    nini=length(ini);
    [m in]=max(abs(val(ini)));
    index =[index ini(in)];
    ptimes=[ptimes tim(ini(in))];
    ini2=[];
    if in<=0.39*nini
       ini2=ini(round(nini*0.55):end);
    end
    if in>=0.61*nini
       ini2=ini(1:round(nini*0.45));
    end
    if ini2,
      [m in]=max(abs(val(ini2)));
      index =[index ini2(in)];
      ptimes=[ptimes tim(ini2(in))];
    end
    i0=inext(i)+1;
  end
  t0=tim(ii(i0));t1=tim(ii(end));
  ptimes=[ptimes mean([t0 t1])];
  index =[index round(mean([ii(i0),ii(end)]))];
end

return
end
