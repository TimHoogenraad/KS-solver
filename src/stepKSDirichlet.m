function [x,flow] = stepKSDirichlet(KS,len,N,init,Cs)


k  = KS.k;
E  = KS.E;
E2 = KS.E2;
Q  = KS.Q;
f1 = KS.f1;
f2 = KS.f2;
f3 = KS.f3;
g  = KS.g;

useSGS = (nargin == 5) && ~isempty(Cs);

% ----- interior Dirichlet grid (what you return) -----
dx = len/(N+1);
x  = (1:N)' * dx;

% ----- doubled periodic grid (what you evolve internally) -----
Np  = 2*(N+1);          % periodic points on [0,2L)
dxp = 2*len/Np;         % equals dx
% periodic x-grid not needed externally

% Initial condition on interior
if isa(init,'function_handle')
    u0 = init(x);
else
    u0 = init(:);
end
assert(numel(u0)==N, 'init must be length N (interior values) or a function handle.');

% Build odd-extended periodic field U on Np points:
% Indices: 1 corresponds to x=0, N+2 corresponds to x=L
U = zeros(Np,1);
U(2:N+1)       = u0;            % (0,L) interior
U(N+2)         = 0;             % x=L
U(N+3:end)     = -flipud(u0);   % (L,2L) mirror
U(1)           = 0;             % x=0

v = fft(U);


if useSGS
    [~, nu_full] = sgs_full(v,k,dxp,Cs);
    Nu(:,1) = nu_full(1:N+2);
end

% ---- main ETDRK4 loop  ----

    % 1
    U  = real(ifft(v));
    Nvl = g .* fft(U.^2);

    if useSGS
        [Nvs, nu_full] = sgs_full(v,k,dxp,Cs);
        Nv = Nvl - Nvs;
    else
        Nv = Nvl;
    end
    a = E2.*v + Q.*Nv;

    % 2
    Ua  = real(ifft(a));
    Nal = g .* fft(Ua.^2);
    if useSGS
        Nas = sgs_full(a,k,dxp,Cs);
        Na = Nal - Nas;
    else
        Na = Nal;
    end
    b = E2.*v + Q.*Na;

    % 3
    Ub  = real(ifft(b));
    Nbl = g .* fft(Ub.^2);
    if useSGS
        Nbs = sgs_full(b,k,dxp,Cs);
        Nb = Nbl - Nbs;
    else
        Nb = Nbl;
    end
    c = E2.*a + Q.*(2*Nb - Nv);

    % 4
    Uc  = real(ifft(c));
    Ncl = g .* fft(Uc.^2);
    if useSGS
        Ncs = sgs_full(c,k,dxp,Cs);
        Nc = Ncl - Ncs;
    else
        Nc = Ncl;
    end

    % update
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;

    % ---- enforce odd symmetry strongly (Dirichlet guarantee) ----
    v = enforceOdd(v,N);

    % save interior
    U = real(ifft(v));


    % stats on interior (same as solveKS.m but restricted)
    u1 = real(ifft( 1i*k .* v));
    u2 = real(ifft(-(k.^2).* v));
    u3 = real(ifft(-(1i*k.^3).*v));
    u4 = real(ifft((k.^4).* v));


if useSGS
    [~, nu_full] = sgs_full(v,k,dxp,Cs);
    flow = [[0;U(2:N+1);0],u1(1:N+2),u2(1:N+2),u3(1:N+2),u4(1:N+2),nu_full(1:N+2)];
else
    flow = [[0;U(2:N+1);0],u1(1:N+2),u2(1:N+2),u3(1:N+2),u4(1:N+2)];
end

end

% -----------------------------
% Enforce odd symmetry on [0,2L) given interior length N on (0,L)
% Indices: 1 -> x=0, N+2 -> x=L
% -----------------------------
function vhat = enforceOdd(vhat,N)
Np = 2*(N+1);
U  = real(ifft(vhat));

% enforce exact odd reflection
U(1)      = 0;
U(N+2)    = 0;
U(N+3:end)= -flipud(U(2:N+1));

vhat = fft(U);

end

% -----------------------------
% SGS (full periodic grid) — same structure as your solveKS.m sgs()
% Returns Ns (spectral) and nu(x) (physical) on full periodic grid
% -----------------------------
function varargout = sgs_full(v,k,delta,Cs)
d1 = 1i.*k;
d3 = 1i.*k.^3;

d = real(ifft(d3 .* v));    % u_xxx in physical space
s = real(ifft(d1 .* v));    % u_x   in physical space

nu = Cs.^2 .* delta.^2 .* abs(s);
Ns = d1 .* fft(nu .* d);    % spectral contribution

if nargout == 1
    varargout{1} = Ns;
else
    varargout{1} = Ns;
    varargout{2} = nu;
end
end
