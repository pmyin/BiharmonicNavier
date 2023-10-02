function u=lapsingular(x,tau,R,omega)
    [th,r] = cart2pol(x(:,1),x(:,2));
    
    ct=find(th<0);
    th(ct)=th(ct)+2*pi;
    
    % the following two parameters should be changed together with the file
    % usingular.m
    R1=tau*2*R;
    R2=2*R;

%    theta=find(th>1*pi)
    
    b1 = find(r<=R1);
    b2 = find(r>=R2);


%     w10 = 2*r/(R2-R1)-(R2+R1)/(R2-R1);
%     c1 = -35/32;
%     c2 = 35/32;
%     c3 = -21/32;
%     c4 = 5/32;
%     w1 = 1/2 + c1*w10 + c2*w10.^3 + c3*w10.^5 + c4*w10.^7;
%     w1r = (c1 + 3*c2*w10.^2 + 5*c3*w10.^4 + 7*c4*w10.^6)*(-2/(R2-R1));
%     w1rr = (6*c2*w10 + 20*c3*w10.^3 + 42*c4*w10.^5)*(-2/(R2-R1))^2;

    w10 = 2*r/(R2-R1)-(R2+R1)/(R2-R1);
    c1 = -15/16;
    c2 = 5/8;
    c3 = -3/16;
    c4 = 0;
    w1 = 1/2 + c1*w10 + c2*w10.^3 + c3*w10.^5 + c4*w10.^7;
    w1r = (c1 + 3*c2*w10.^2 + 5*c3*w10.^4 + 7*c4*w10.^6)*(2/(R2-R1));
    w1rr = (6*c2*w10 + 20*c3*w10.^3 + 42*c4*w10.^5)*(2/(R2-R1))^2;

    w2 = r.^(-pi/omega);
    w2r = -(pi/omega)*r.^(-pi/omega-1);
    w2rr = -(pi/omega)*(-pi/omega-1)*r.^(-pi/omega-2);
    
    
    w3=sin(pi/omega*th);
    
    wtt= -(pi/omega)^2*w1.*w2.*w3;
    
    wr= w1r.*w2.*w3 + w1.*w2r.*w3;
    
    wrr= w1rr.*w2.*w3 +2*w1r.*w2r.*w3 + w1.*w2rr.*w3;
    
    u = wtt./(r.^2)+wr./r+wrr;

    u(b1,:)=0;
    u(b2,:)=0;
    
end

