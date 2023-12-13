addpath('~/tools/m_map'); %addpath('../');

% geo-boundary of our study area
latlim=[34 71]; lonlim=[-12 37]; % entire Europe

% time period of interest kyrBP
timelim=[3.1 9.9];
maxradius=920; %%800;
dtb=200;
breaks=3000:(2*dtb):9800;

% name of data treatments in SPD generation
tags={'_Norm','_NoNorm','_NoNorm_Bin100','_Norm_Bin100'};

%  output directory
scdir='out/';
%  work directory
sdir=pwd; sdir=[sdir '/'];

ccol=[223 233 130]/256;
mcol(1,:)=[1 0.5 0]; mcol(2,:)=[0 0.5 1];
