function [lab,V] = eigK_with_vec(x,K)
%-------------------------------------------------------------------------------------------------
% This file is a modified version of the "eigK.m" included in the SeDuMi package.
% We modify the interface to return the corresponding eigenvectors in addition to the eigenvalues.
%-------------------------------------------------------------------------------------------------

% Complete rewrite for SeDuMi 1.3 by Michael C. Grant
% Copyright (C) 2013 Michael C. Grant.
%
% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA

% The existence of rsdpN is code for 'is this the internal format?'
if isfield(K,'rsdpN')
    % The internal SeDuMi cone format.
    is_int = true;
    nf = 0;
    nl = K.l;
    nq = length(K.q);
    nr = 0;
    ns = length(K.s);
    nrsdp = K.rsdpN;
    N = nl+2*length(K.q)+sum(K.s);
else
    % The external SeDuMi cone format.
    is_int = false;
    N = 0;
    if isfield(K,'f'), nf = K.f; else nf = 0; end
    if isfield(K,'l'), nl = K.l; N = N + nl; else nl = 0; end
    if isfield(K,'q'), nq = length(K.q); N = N + 2 * nq; else nq = 0; end
    if isfield(K,'r'), nr = length(K.r); N = N + 2 * nr; else nr = 0; end
    if isfield(K,'s'), ns = length(K.s); N = N + sum(K.s); else ns = 0; end
    if isfield(K,'z'), ns = ns + length(K.z); N = N + sum(K.z); end
    nrsdp = ns;
end
li = 0;
xi = nf;
lab = zeros(N,1);
lab(li+1:li+nl) = x(xi+1:xi+nl);
xi = xi + nl;
li = li + nl;
if nq
    tmp = sqrt(0.5);
    if is_int
        % Internal version: all of the x0 values are at the front, and the
        % vectors are stacked in order after that.
        zi = xi;
        xi = xi + nq;
        for i = 1:nq
            kk = K.q(i) - 1;
            x0 = x(zi+i);
            lab(li+1:li+2) = tmp*(x0+[-1;+1]*norm(x(xi+1:xi+kk)));
            xi = xi + kk;
            li = li + 2;
        end
    else
        for i = 1:length(K.q)
            kk = K.q(i);
            x0 = x(xi+1);
            lab(li+1:li+2) = tmp*(x0+[-1;+1]*norm(x(xi+2:xi+kk)));
            xi = xi + kk;
            li = li + 2;
        end
    end
end
for i = 1:nr
    % Only the external format needs to be implemented here
    ki = K.r(i);
    x1 = xx(xi+1);
    x2 = xx(xi+2);
    lab(li+1:li+2) = 0.5*(x1+x2+[-1;+1]*norm([x1-x2;2*x(xi+3:xi+ki)]));
    xi = xi + ki;
    li = li + 2;
end

V = cell(ns,1);
for i = 1 : ns
    ki = K.s(i);
    qi = ki * ki;
    XX = x(xi+1:xi+qi); 
    xi = xi + qi;
    if i > nrsdp
        XX = XX + 1i*x(xi+1:xi+qi); 
        xi = xi + qi;
    end
    XX = reshape(XX,ki,ki);
    XX = XX + XX';
    [P,D] = eig( 0.5*XX );
    V{i,1} = P;
    lab(li+1:li+ki) = diag(D);
    li = li + ki;
end
