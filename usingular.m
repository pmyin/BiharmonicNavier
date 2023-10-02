function u=usingular(x,tau,R,omega)
    [th,r] = cart2pol(x(:,1),x(:,2));
    
    ct=find(th<0);
    th(ct)=th(ct)+2*pi;
    
    % the following two parameters should be changed together with the file
    % lapsingular.m
    R1=tau*2*R;
    R2=2*R;
    
    
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

    w2=r.^(-pi/omega);
    
    w3=sin(th*pi/omega);

    b3 = find (r==0);
    w2(b3)=(1e-32).^(-pi/omega);
%    w2(b3)=0;
    w3(b3) = 1;
    if length(b3)>0
        disp('b3 not 0')
    end
    
%     if ( length(b3)>0 & length(x)==length(p) )
%         for i=1:length(b3)
%             for j=1:length(t)
%                 for k=1:length(t(1,:))
%                     if ( b3(i)==t(j,k) )
%                         if (k==1)
%                             k1 = 2; k2 = 3;
%                         end
%                         if (k==2)
%                             k1 = 1; k2 = 3;
%                         end
%                         if (k==3)
%                             k1 = 1; k2 = 2;
%                         end
%                         th(b3(i)) = (th(t(j,k1))+th(t(j,k2)))/2;
%                     end
%                 end
%             end
%         end
%     end
    
    u = w1.*w2.*w3;
end