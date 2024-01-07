function nsum = nansum(arg)
arg(isnan(arg)) = 0;
nsum = sum(arg);
