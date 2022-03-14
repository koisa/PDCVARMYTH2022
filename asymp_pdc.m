%% ASYMP_PDC
%        Compute partial directed coherence magnitude for all three
%        metrics -- Euclidean, diagonal and informational -- as well as asymptotic 
%        statistics from vector autoregressive (VAR) coefficients in the frequency
%        domain.
%% Syntax
%        c = ASYMP_PDC(u,A,pf,nFreqs,metric,alpha)
%
%% Input Arguments
%        u:      time series
%        A:      AR estimate matrix by MVAR routine
%        pf:     covariance matrix provided by MVAR routine
%        nFreqs: number of point in [0,fs/2] frequency scale
%        metric: 'euc'  -- Euclidean   ==> original PDC
%                'diag' -- diagonal    ==> gPDC (generalized)
%                'info' -- information ==> iPDC
%        alpha:  significance level
%                 if alpha is zero, statistical analysis not performed
%
%% Output Arguments 
%        c structure variable with following fields:
%              |-- .pdc2      - |PDC|^2 estimates
%              |-- .cpdc      - complex PDC
%              |-- .pvalues   - p-values associated to PDC2 estimates. 
%              |-- .th        - Threshold value with (1-avalue) significance level.
%              |-- .{ic1,ic2} - upper and lower (1 - alpha) confidence interval of |PDC|^2 estimates
%              |-- .metric    - metric used for PDC calculation 
%              |-- .alpha     - significance level
%              |-- .p         - VAR model order
%              |-- .patdenr   - 
%              |-- .patdfr    - degree of freedom 
%              |-- .SS        - power spectra
%              +-- .coh2      - squared spectral coherence
%      or
%   c.{pdc2,cpdc,pvalues,th,ic1,ic2,metric,alpha,p,patdenr,patdfr,SS,coh2}
%
%% Description
%   Compute any of three types of $|PDC|^2$ --- connectivity measure --- as well
%   as its  allied asymptotic statistics [2] measures for chosen metric option: 
%        * 'euc'  - original or Euclidean PDC as proposed in [1]; 
%        * 'diag' - generalized PDC; 
%        * 'info' - information PDC. 
%
%% References
%
% [1] Baccala LA and Sameshima K (2001). Partial directed coherence: a new
%     concept in neural structure determination. Biol Cybern 84:463--474.
%                   <https://doi.org/10.1007/PL00007990>
%
% [2] Takahashi DY, Baccala LA and Sameshima K (2007). Connectivity
%     inference between neural structures via partial directed coherence. 
%     J Appl Stat 34:1259--1273. 
%               <https://doi.org/10.1080/02664760701593065>
%
% See also PDC_TOT_P, MVAR, MCARNS,
%          | <pdc_tot_p.html> | <mvar.html> | <mcarns.html> |
%

% (C) Koichi Sameshima & Luiz A. Baccala, 2021. See file Readme.pdf in
% installation directory for licensing terms.
%
%%

function c = asymp_pdc(u,A,pf,nFreqs,metric,alpha)
if nargin < 6
    error('ASYMP_PDC requires six input arguments.')
end
[m,Ndata] = size(u);
if m > Ndata
    u = u.';
    [m,Ndata] = size(u);
end

[nChannels, ~ ,p] = size(A); 
if nChannels ~= m, 
   error('Numbers of channels for data set and its A coefficients do not match.')
end

Af = A_to_f(A, nFreqs);

flgVerbose = 0;

% Variables pre-alocation
pdc2 = zeros(nChannels,nChannels,nFreqs);

if alpha ~= 0
   th  = zeros(nChannels,nChannels,nFreqs);
   ic1 = zeros(nChannels,nChannels,nFreqs);
   ic2 = zeros(nChannels,nChannels,nFreqs);
   varass1 = zeros(nChannels,nChannels,nFreqs);
   varass2 = zeros(nChannels,nChannels,nFreqs);
   patdfr = zeros(nChannels,nChannels,nFreqs);
   patdenr = zeros(nChannels,nChannels,nFreqs);
   pvalues = zeros(nChannels,nChannels,nFreqs);
   if flgVerbose
      switch lower(metric)
         case {'euc'}
            disp('* Original PDC and asymptotic statistics')
         case {'diag'}
            disp('* Generalized PDC and asymptotic statistics')
         case {'info'}
            disp('* Information PDC and asymptotic statistics')
         otherwise
            error('Unknown metric.')
      end
   end
elseif flgVerbose
   switch lower(metric)
      case {'euc'}
         disp('* Original PDC estimation')
      case {'diag'}
         disp('* Generalized PDC estimation')
      case {'info'}
         disp('* Information PDC estimation')
      otherwise
         error('Unknown metric.')
   end
end

% ======================================================================

switch lower(metric)
   case {'euc'}           % for original PDC
      dpdc_dev = zeros(1,(nChannels*(nChannels + 1))/2);
      
   case {'diag'}          % for gPDC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2*nChannels), evar_d);
      pinv_evar_d_big = pinv(evar_d_big);

      %'derivative of vec(Ed-1) by vecE'
      de_deh = Dup(nChannels);
      debig_de = fdebig_de(nChannels);
      
      dedinv_dev = diagtom(vec(-pinv_evar_d_big*pinv_evar_d_big));
      dedinv_deh = dedinv_dev*debig_de*de_deh;
      
   case {'info'}          % for iPDC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2*nChannels), evar_d);
      pinv_evar_d_big = pinv(evar_d_big);
      
      evar_big = kron(eye(2*nChannels), pf);
      pinv_evar_big = sparse(pinv(evar_big));

      %'derivative of vec(Ed-1) by vecE'
      de_deh = Dup(nChannels);
      debig_de = fdebig_de(nChannels);
      
      dedinv_devd = sparse(diagtom(vec(-pinv_evar_d_big*pinv_evar_d_big)));
      dedinv_dehd = sparse(dedinv_devd*debig_de*de_deh);
      
      %dedinv_dev = sparse(-kron(inv_e.', inv_e));
      dedinv_dev = sparse(-kron(pinv_evar_big.', pinv_evar_big));
      dedinv_deh = sparse(dedinv_dev*debig_de*de_deh);

   otherwise
      error('Unknown metric.')
end

gamma = bigautocorr(u, p);
omega = kron(inv(gamma), pf);
omega_evar = 2*pinv(Dup(nChannels))*kron(pf, pf)*pinv(Dup(nChannels)).';

%if (is_octave)
icdf_norm_alpha = norminv(1 - alpha/2.0,0,1);
%else
%   icdf_norm_alpha = icdf('norm',1 - alpha/2.0,0,1);
%end

for ff = 1:nFreqs
   f = (ff - 1)/(2*nFreqs); %Corrected 7/25/2011, f starting at 0 rad/s.
   Ca = fCa(f, p, nChannels);
   
   a = Af(ff,:,:); a = a(:);    %Equivalent to a = vec(Af[ff, :, :])
   a = [real(a); imag(a)];      %a = cat(a.real, a.imag, 0)

   omega2 = Ca*omega*Ca.';
   L = fChol(omega2);

   for j = 1:nChannels
      Ij = fIj(j,nChannels);
      switch lower(metric)
         case {'euc'}           % for PDC
            pfe=eye(nChannels);
            Ije = Ij;
         case {'diag'}          % for gPDC
            Ije = Ij*pinv_evar_d_big;
            pfe=pf;
         case {'info'}          % for iPDC
            Ije = Ij*pinv_evar_big*Ij;
            pfe=pf;
         otherwise
            error('Unknown metric.')
      end
      
      for i = 1:nChannels
         Iij = fIij(i,j,nChannels);
         %For metric = diag or info case, include evar in the expression'
         switch lower(metric)
            case {'euc'}           % for PDC
               Iije = Iij;
            case {'diag'}          % for gPDC
               Iije = Iij*pinv_evar_d_big;
            case {'info'}          % for iPDC
               Iije = Iij*pinv_evar_d_big;
            otherwise
               error('Unknown metric.')
         end

         num = a.'*Iije*a;
         den = a.'*Ije*a;
         pdc2(i,j,ff) = num/den; % $|PDC_{ij}(ff)|^2$
          numd=Af(ff,i,j)/sqrt(pf(i,i));
          cpdc(i,j,ff)=numd/sqrt(den);
         
         % If alpha == 0, no statistical calculation for faster computation.
         if alpha ~= 0
            %'Add evar derivation'
            switch lower(metric)
               case {'euc'}
                  %nop
               case {'diag'}
                  %'derivative of num by vecE'
                  dnum_dev = kron((Iij*a).', a.')*dedinv_deh;
                  %'derivative by den by vecE'
                  dden_dev = kron((Ij*a).', a.')*dedinv_deh;
                  dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2);
               case {'info'}
                  %'derivative of num by vecE'
                  dnum_dev = kron((Iij*a).', a.')*dedinv_dehd;
                  %'derivative of den by vecE'
                  dden_dev = kron((Ij*a).', a.'*Ij)*dedinv_deh;
                  dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2);
               otherwise
                  error('Unknown metric.')
            end

            G1a = 2 * a.' * Iije/den - 2*num*a.'*Ije/(den^2);
            G1 = -G1a * Ca;
            varalpha = G1 * omega * G1.';
            varevar = dpdc_dev * omega_evar * dpdc_dev.';
            varass1(i,j,ff) = (varalpha + varevar)/Ndata;
            
            ic1(i,j,ff) = pdc2(i,j,ff) ...
                    - sqrt(varass1(i,j,ff))*icdf_norm_alpha;
            ic2(i,j,ff) = pdc2(i,j,ff) ...
                    + sqrt(varass1(i,j,ff))*icdf_norm_alpha;

            G2a = Iije/den;
                        
            d = fEig(real(L), real(G2a)); % real() 28May2013            

            patdf = (sum(d).^2)./sum(d.^2);
            patden = sum(d)./sum(d.^2);

            th(i,j,ff) = chi2inv((1 - alpha), patdf)/(patden*Ndata);
            pvalues(i,j,ff) = 1 - chi2cdf(pdc2(i,j,ff)*patden*Ndata, patdf);
 
            varass2(i,j,ff) = patdf/(patden*Ndata).^2;
            patdfr(i,j,ff) = patdf;
            patdenr(i,j,ff) = patden;
         else % alpha is zero, no asymptotics calculation
            % nop;
         end
      end
   end
end

% assigning values to the c structure.
if alpha ~= 0
   c.pdc2 = pdc2;
   c.cpdc=cpdc;
   c.th = th;
   c.ic1 = ic1;
   c.ic2 = ic2;
   c.metric = metric;
   c.alpha = alpha;
   c.p = p;
   c.pvalues = pvalues; % p-values associated to PDC/gPDC/iPDC
   c.patden = patdenr;
   c.patdf = patdfr;
   c.varass1 = varass1;
   c.varass2 = varass2;   
%  Statistically significant PDC on frequency scale
   pdc2_temp = ((abs(pdc2) - abs(th)) > 0).*pdc2 + ((abs(pdc2) - abs(th)) <= 0)*(-1);
   pdc2_temp(pdc2_temp < 0) = NaN;
   c.pdc2_th = pdc2_temp;
else % No statistics, just the |PDC^2| values/
   c.pdc2 = pdc2;
   c.cpdc=cpdc;
   c.metric = metric;
   c.alpha = 0;
   c.p = p;
   c.pvalues = [];   
   c.th = [];
   c.ic1 = [];
   c.ic2 = [];
   c.patdenr = [];
   c.patdfr = [];
   c.varass1 = [];
   c.varass2 = [];
end

% Power spectra and coherence assignments to c structure.
c.SS = ss_alg(A, pf, nFreqs);
c.coh2 = coh_alg(c.SS);
%
end

%==========================================================================
function gamma = bigautocorr(x, p)
%Autocorrelation. Data in row-wise orientation. From order 0 to p-1.
%Output: n x n blocks of autocorr of lags i. (Nuttall Strand matrix)'''
[n, nd] = size(x);

gamma = zeros(n*p, n*p);
for i = 1:p
   for j = 1:p
      gamma(((i - 1)*n + 1):i*n, ((j - 1)*n + 1):j*n) = xlag(x, i - 1)*(xlag(x,j - 1).')/nd;
   end
end
end

%==========================================================================
function c = xlag(x,tlag)
if tlag == 0
   c = x;
else
   c = zeros(size(x));
   c(:,(tlag + 1):end) = x(:,1:(end - tlag));
end
end

%==========================================================================
function d = fEig(L, G2)
%'''Returns the eigenvalues'''

%L = mat(cholesky(omega, lower=1))
D = L.'*G2*L;

%disp('fEig: eig or svd?')

d  = svd(D);
d1 = sort(d);
%
% the two biggest eigenvalues no matter which values (non negative by
% construction
%
d = d1(length(d) - 1:length(d));

if (size(d) > 2)
   disp('More than two Chi-squares in the sum:')
end
end

%==========================================================================
function c = fIij(i,j,n)
%'''Returns Iij of the formula'''
Iij = zeros(1,n^2);
Iij(n*(j - 1) + i) = 1;
Iij = diag(Iij);
c = kron(eye(2), Iij);

c = sparse(c);  % SPARSED
end

%==========================================================================
function c = fIj(j,n)
%'''Returns Ij of the formula'''
Ij = zeros(1,n);
Ij(j) = 1;
Ij = diag(Ij);
Ij = kron(Ij,eye(n));
c = kron(eye(2), Ij);

c = sparse(c);  % SPARSED
end

%==========================================================================
function d = fCa(f, p, n)
%'''Returns C* of the formula'''
C1 = cos(-2*pi*f*(1:p));
S1 = sin(-2*pi*f*(1:p));
C2 = [C1; S1];
d = kron(C2, eye(n^2));
end

%==========================================================================
function c = fdebig_de(n)
%'''Derivative of kron(I(2n), A) by A'''
%c = kron(TT(2*n, n), eye(n*2*n)) * kron(eye(n), kron(vec(eye(2*n)), eye(n)));
A = sparse(kron(TT(2*n, n), eye(n*2*n)));
B = sparse(kron(vec(eye(2*n)), eye(n)));
c = A * kron(eye(n), B);
c = sparse(c);
end

%==========================================================================
function c = vec(x)
%vec = lambda x: mat(x.ravel('F')).T
c = x(:);
end

%==========================================================================
function t = TT(a,b)
%''' TT(a,b)*vec(B) = vec(B.T), where B is (a x b).'''
t = zeros(a*b);
for i = 1:a
   for j = 1:b
      t((i - 1)*b + j,(j - 1)*a + i) = 1;
   end
end
t = sparse(t);
end

%==========================================================================
function L = fChol(omega)
% Try Cholesky factorization
try
   L = chol(omega)';
   % If there's a small negative eigenvalue, diagonalize
catch
   %   disp('linalgerror, probably IP = 1.')
   [v,d] = eig(omega);
   L = zeros(size(v));
   for i = 1:length(d)
      if d(i,i) < 0
         d(i,i) = eps;
      end
      L(:,i) = v(:,i)*sqrt(d(i,i));
   end
end
end

%==========================================================================
function c = diagtom(a)
b = sparse(a');
c = sparse(diag(b(:)));
end

%==========================================================================
function c = mdiag(a)
%  diagonal matrix
c = diag(diag(a)); % UNSPARSED
end

%==========================================================================
function d=Dup(n)
% '''D*vech(A) = vec(A), with symmetric A'''
d = zeros(n*n, (n*(n + 1))/2);
count = 1;
for j = 1:n
   for i = 1:n
      if i >= j
         d((j - 1)*n + i,count) = 1;
         count = count + 1;
      else
         d((j - 1)*n + i,:) = d((i - 1)*n + j,:);
      end
   end
end
end

%==========================================================================

%% Change Log:
%Corrected 7/25/2011 to match the frequency range with plotting
%routine, f = 0 was include in the frequency for loop:
%                                for ff = 1:nFreqs,
%                                   f = (ff - 1)/(2*nFreqs); % 
%                                          ^?^^
% Optimization 7/1/2015 \(^o^)/

% [EOF]
