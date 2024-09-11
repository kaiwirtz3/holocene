%% SPDX-FileCopyrightText: 2012-2024 Helmholtz-Zentrum hereon GmbH
%% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
%% SPDX-License-Identifier: CC0-1.0
%
% NORMINV(p,mu,sigma) returns the inverse of the normal cumulative
%  density function with mean mu and standard deviation sigma,
%  evaluated at the probability values in p.

function mp=norminv(p,mu,sigma)
  mp=sigma.*(-sqrt(2)*erfcinv(2*p))+mu;
  return
