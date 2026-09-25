function [x,flow] = stepKSPeriodic(KS,len,N,init,Cs)

% One timesetp for Periodic KS equation

% Spatial grid and initial condition:
x = len*(1:N)'/N;
u = init;
v = fft(u);
dx = x(1);


k  = KS.k;
E  = KS.E;
E2 = KS.E2;
Q  = KS.Q;
f1 = KS.f1;
f2 = KS.f2;
f3 = KS.f3;
g  = KS.g;
advection = KS.advection;


if (nargin == 5) % LES

    % 1
    [Nvs, ~] = sgs(v,k,dx,Cs);
    Nvl = transport(v,k,g,advection);
    Nv = Nvl - Nvs;
    a = E2.*v + Q.*Nv;
    % 2
    Nas = sgs(a,k,dx,Cs);
    Nal = transport(a,k,g,advection);
    Na = Nal - Nas;
    b = E2.*v + Q.*Na;
    % 3
    Nbs = sgs(b,k,dx,Cs);
    Nbl = transport(b,k,g,advection);
    Nb = Nbl - Nbs;
    stageC = E2.*a + Q.*(2*Nb-Nv);
    % 4
    Ncs = sgs(stageC,k,dx,Cs);
    Ncl = transport(stageC,k,g,advection);
    Nc = Ncl - Ncs;
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    %v(N/2+1) = 0;

    % Return to real space and save solution
    u = real(ifft(v));
    u1  = real(ifft( 1i*k .* v));
    u2  = real(ifft(-(k.^2) .* v));
    u3  = real(ifft( -(1i*k.^3).* v ));
    u4  = real(ifft( (k.^4) .* v));
    [~, nu] = sgs(v,k,dx,Cs);

    flow = [u,u1,u2,u3,u4,nu];



else % DNS

    Nv = transport(v,k,g,advection);
    a = E2.*v + Q.*Nv;
    Na = transport(a,k,g,advection);
    b = E2.*v + Q.*Na;
    Nb = transport(b,k,g,advection);
    stageC = E2.*a + Q.*(2*Nb-Nv);
    Nc = transport(stageC,k,g,advection);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    %v(N/2+1) = 0;
    % Return to real space and save solution
    u = real(ifft(v));
    u1  = real(ifft( 1i*k .* v));
    u2  = real(ifft(-(k.^2) .* v));
    u3  = real(ifft( -(1i*k.^3).* v ));
    u4  = real(ifft( (k.^4) .* v));

    flow = [u,u1,u2,u3,u4];

end

end

function N = transport(v,k,g,advection)
u = real(ifft(v));
N = g.*fft(u.^2) - advection.*(1i.*k).*v;
end

function [Ns,nu] = sgs(v,k,delta,Cs)

d1 = 1i.*k;
d3 = 1i.*k.^3;
d  = real(ifft(d3 .* v));
s  = real(ifft(d1 .* v));
nu = Cs.^2 .* delta.^2 .* abs(s);
Ns = d1.*fft(nu .* d);

end
