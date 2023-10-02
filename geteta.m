function u=geteta(r,tau,R)
    
    % the following two parameters should be changed together with the file
    % lapsingular.m
    R1=tau*R;
    R2=R;
    
    
    b1 = find(r<=R1);
    b2 = find(r>=R2);
%    w1=exp((a*a)./(r.*(r-2*a))+1);


%     w10 = 2*r/(R2-R1)-(R2+R1)/(R2-R1);
%     c1 = -35/32;
%     c2 = 35/32;
%     c3 = -21/32;
%     c4 = 5/32;
%     w1 = 1/2 + c1*w10 + c2*w10.^3 + c3*w10.^5 + c4*w10.^7;

    w10 = 2*r/(R2-R1)-(R2+R1)/(R2-R1);
    c1 = -15/16;
    c2 = 5/8;
    c3 = -3/16;
    c4 = 0;
    w1 = 1/2 + c1*w10 + c2*w10.^3 + c3*w10.^5 + c4*w10.^7;

    w1(b1)=1;
    w1(b2)=0;
    
    u = w1;
end