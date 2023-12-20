%
% calculate moving average
%
% kai wirtz (hereon) Dec 2023
%
function [arg1 arg2]=movavg(times,values,window)
twin=window/2;
nut=length(times);
for it = 1:nut
  ind=find(times>times(it)-twin & times<times(it)+twin);
  mavg(it)=nanmean(values(ind));
end
%movavg=values-mavg;
if nargout>1
  arg2=mavg;
  arg1=times;
else
  arg1=mavg;
end

return
