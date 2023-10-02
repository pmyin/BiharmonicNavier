function u = getsint(omega,tau,R,option)

u1 = omega/(4-4*pi/omega)*(tau*R)^(2-2*pi/omega);

quadorder = 10;

[GaussPnt, GaussWeight, Point]=buildGauss1d(quadorder);


maxIt = option.maxIt;
Ndim = 2^maxIt;
h= (R-tau*R)/Ndim;
LeftPoint = [tau*R:h:R-h]';
RighPoint = [tau*R+h:h:R]';

jacobian=h/(Point(2)-Point(1));

GlobalPnt=Local_to_Global1d(GaussPnt,[LeftPoint RighPoint]);

nQuad = size(GaussWeight,1)
u2loc = zeros(Ndim,1);
for i=1:nQuad
    u2loc(:)=u2loc(:)+geteta(GlobalPnt(:,i),tau,R).*geteta(GlobalPnt(:,i),tau,R).*GlobalPnt(:,i).^(1-2*pi/omega)*GaussWeight(i)*jacobian;
end

u2=sum(u2loc)*omega/2;

u =u1+u2;
