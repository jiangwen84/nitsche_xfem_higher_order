function[jcob,Pos,Weight] = get_4thOrderGaussQuad_info(quadNod,iPsi,iEta)

%% Parameters for fourth-order Gauss quadrature
r = sqrt(1.2);
gaussLocPos(1) = -sqrt((3 + 2*r)/7);     gaussWeight(1) = (1/2) - (6*r)^-1;
gaussLocPos(2) = -sqrt((3 - 2*r)/7);     gaussWeight(2) = (1/2) + (6*r)^-1;
gaussLocPos(3) =  sqrt((3 - 2*r)/7);     gaussWeight(3) = (1/2) + (6*r)^-1;
gaussLocPos(4) =  sqrt((3 + 2*r)/7);     gaussWeight(4) = (1/2) - (6*r)^-1;

%% Gauss weights for the considered Gauss point
Weight = struct('psi',gaussWeight(iPsi),'eta',gaussWeight(iEta));

%% Locale coordinates of the Gauss point
psi = gaussLocPos(iPsi);
eta = gaussLocPos(iEta);

%% Local shape functions evaluated at (psi,eta)
Nloc = zeros(4,1);
Nloc(1) = (1/4)*(1 - psi)*(1 - eta);
Nloc(2) = (1/4)*(1 + psi)*(1 - eta);
Nloc(3) = (1/4)*(1 + psi)*(1 + eta);
Nloc(4) = (1/4)*(1 - psi)*(1 + eta);

%% Derivatives of the local shape functions evaluated at (psi,eta,chi)
dNdpsi = (1/4)*[-(1-eta),(1-eta),(1+eta),-(1+eta)];
dNdeta = (1/4)*[-(1-psi),-(1+psi),(1+psi),(1-psi)];

%% Derivatives of the global coordinates w.r.t the local coordiantes
dxdpsi = dNdpsi*[quadNod.X]';
dxdeta = dNdeta*[quadNod.X]';
dydpsi = dNdpsi*[quadNod.Y]';
dydeta = dNdeta*[quadNod.Y]';

%% Jacobian of the local->global transformation
J = [dxdpsi,dydpsi
    dxdeta,dydeta];

jcob = abs(det(J));

%% Global position of the considered Gauss point
Pos = struct('X',0.0,'Y',0.0,'Z',0.0);
for i=1:4
    Pos.X = Pos.X + Nloc(i)*quadNod(i).X;
    Pos.Y = Pos.Y + Nloc(i)*quadNod(i).Y;
end
end