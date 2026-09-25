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


if (nargin == 5) % LES

    % 1
    [Nvs, ~] = sgs(v,k,dx,Cs);
    Nvl = g.*fft(real(ifft(v)).^2);
    Nv = Nvl - Nvs;
    a = E2.*v + Q.*Nv;
    % 2
    Nas = sgs(a,k,dx,Cs);
    Nal = g.*fft(real(ifft(a)).^2);
    Na = Nal - Nas;
    b = E2.*v + Q.*Na;
    % 3
    Nbs = sgs(b,k,dx,Cs);
    Nbl = g.*fft(real(ifft(b)).^2);
    Nb = Nbl - Nbs;
    c = E2.*a + Q.*(2*Nb-Nv);
    % 4
    Ncs = sgs(c,k,dx,Cs);
    Ncl = g.*fft(real(ifft(c)).^2);
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

    Nv = g.*fft(real(ifft(v)).^2);
    a = E2.*v + Q.*Nv;
    Na = g.*fft(real(ifft(a)).^2);
    b = E2.*v + Q.*Na;
    Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = g.*fft(real(ifft(c)).^2);
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

function [Ns,nu] = sgs(v,k,delta,Cs)

d1 = 1i.*k;
d3 = 1i.*k.^3;
d  = real(ifft(d3 .* v));
s  = real(ifft(d1 .* v));
nu = Cs.^2 .* delta.^2 .* abs(s);
Ns = d1.*fft(nu .* d);

end
