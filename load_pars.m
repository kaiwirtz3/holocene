% sets common parameters

addpath('~/tools/m_map'); %addpath('../');

% geo-boundary of our study area
latlim=[34 71]; lonlim=[-12 37]; % entire Europe

% time period of interest kyrBP
timelim=[3 9.3];
dtb=200;
breaks=3000:(2*dtb):9800;

% name of data treatments in SPD generation
tags={'_NoNorm_Bin100','_Norm_Bin100'}; %'_Norm','_NoNorm',

%  output directory
scdir='out/';
%  work directory
sdir=pwd; sdir=[sdir '/'];

ccol=[223 233 130]/256;
mcol(1,:)=[1 0.5 0]; mcol(2,:)=[0 0.5 1];
abc='abcd';
