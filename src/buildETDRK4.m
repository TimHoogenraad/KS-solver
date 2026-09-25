function KS = buildETDRK4(boundary, len, N, h, v2, v4, advection)
% Coefficients used by the latest periodic and Dirichlet long-run solvers.
if nargin < 5, v2 = 1; end
if nargin < 6, v4 = 1; end
if nargin < 7, advection = 0; end
if ~isscalar(advection) || ~isreal(advection) || ~isfinite(advection)
    error('advection must be a finite real scalar');
end

switch lower(char(boundary))
    case 'periodic'
        nfft = N;
        domain = len;
    case 'dirichlet'
        nfft = 2*(N+1);
        domain = 2*len;
    otherwise
        error('boundary must be periodic or dirichlet');
end

if N < 2 || N ~= floor(N) || mod(nfft,2) ~= 0
    error('N must give an even Fourier grid with at least two points');
end

k = [0:(nfft/2-1) 0 (-nfft/2+1):-1]' * (2*pi/domain);
L = v2*k.^2 - v4*k.^4;
E = exp(h*L);
E2 = exp(h*L/2);
M = 16;
r = exp(1i*pi*((1:M)-0.5)/M);
LR = h*L(:,ones(M,1)) + r(ones(nfft,1),:);

KS.k = k;
KS.E = E;
KS.E2 = E2;
KS.Q = h*real(mean((exp(LR/2)-1)./LR,2));
KS.f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2));
KS.f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3,2));
KS.f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2));
KS.g = -0.5*v2*1i*k;
KS.advection = advection;
end
