function  [GaussPnt, GaussWeight, Point]=buildGauss1d(n)

Point=zeros(2,1);

Point(1,1)=-1.0;
Point(2,1)=1.0;

switch n
    
    case 3 
    GaussPnt=zeros(3, 1);
    GaussWeight=zeros(3, 1);
    GaussPnt(1) = 0;
    GaussPnt(2) = -sqrt(3/5);
    GaussPnt(3) = sqrt(3/5);
    GaussWeight(1) = 8/9;
    GaussWeight(2) = 5/9;
    GaussWeight(3) = 5/9;
    
    case 30 
    GaussPnt=zeros(3, 1);
    GaussWeight=zeros(3, 1);
    GaussPnt(1) = 0;
    GaussPnt(2) = -0.7745966692414834;
    GaussPnt(3) = 0.7745966692414834;
    GaussWeight(1) = 0.8888888888888888;
    GaussWeight(2) = 0.5555555555555556;
    GaussWeight(3) = 0.5555555555555556;
 
    case 4 
    GaussPnt=zeros(4, 1);
    GaussWeight=zeros(4, 1);
    GaussPnt(1) = -0.3399810435848563;
    GaussPnt(2) = 0.3399810435848563;
    GaussPnt(3) = -0.8611363115940526;
    GaussPnt(4) = 0.8611363115940526;
    GaussWeight(1) = 0.6521451548625461;
    GaussWeight(2) = 0.6521451548625461;
    GaussWeight(3) = 0.3478548451374538;
    GaussWeight(4) = 0.3478548451374538;

    
    case 5
    GaussPnt=zeros(5, 1);
    GaussWeight=zeros(5, 1);
    GaussPnt(1) = 0;
    GaussPnt(2) = -0.5384693101056831;
    GaussPnt(3) = 0.5384693101056831;
    GaussPnt(4) = -0.906179845938664;
    GaussPnt(5) = 0.906179845938664;
    GaussWeight(1) = 0.5688888888888889;
    GaussWeight(2) = 0.4786286704993665;
    GaussWeight(3) = 0.4786286704993665;
    GaussWeight(4) = 0.2369268850561891;
    GaussWeight(5) = 0.2369268850561891;    
    
    case  6
    GaussPnt=zeros(6, 1);
    GaussWeight=zeros(6, 1);
    GaussPnt(1) = 0.6612093864662645;
    GaussPnt(2) = -0.6612093864662645;
    GaussPnt(3) = -0.2386191860831969;
    GaussPnt(4) = 0.2386191860831969;
    GaussPnt(5) = -0.9324695142031521;
    GaussPnt(6) = 0.9324695142031521;
    GaussWeight(1) = 0.3607615730481386;
    GaussWeight(2) = 0.3607615730481386;
    GaussWeight(3) = 0.467913934572691;
    GaussWeight(4) = 0.467913934572691;
    GaussWeight(5) = 0.1713244923791704;
    GaussWeight(6) = 0.1713244923791704;    

    case  10
    GaussPnt=zeros(10, 1);
    GaussWeight=zeros(10, 1);
    GaussPnt(1) = -0.1488743389816312;
    GaussPnt(2) = 0.1488743389816312;
    GaussPnt(3) = -0.4333953941292472;
    GaussPnt(4) = 0.4333953941292472;
    GaussPnt(5) = -0.6794095682990244;
    GaussPnt(6) = 0.6794095682990244;
    GaussPnt(7) = -0.8650633666889845;
    GaussPnt(8) = 0.8650633666889845;
    GaussPnt(9) = -0.9739065285171717;
    GaussPnt(10) = 0.9739065285171717;
    GaussWeight(1) = 0.2955242247147529;
    GaussWeight(2) = 0.2955242247147529;
    GaussWeight(3) = 0.2692667193099963;
    GaussWeight(4) = 0.2692667193099963;
    GaussWeight(5) = 0.219086362515982;
    GaussWeight(6) = 0.219086362515982;
    GaussWeight(7) = 0.1494513491505806;
    GaussWeight(8) = 0.1494513491505806;
    GaussWeight(9) = 0.0666713443086881;
    GaussWeight(10) = 0.0666713443086881;    
    
    case  20
    GaussPnt=zeros(20, 1);
    GaussWeight=zeros(20, 1);
    GaussPnt(1) = -0.0765265211334973;
    GaussPnt(2) = 0.0765265211334973;
    GaussPnt(3) = -0.2277858511416451;
    GaussPnt(4) = 0.2277858511416451;
    GaussPnt(5) = -0.3737060887154195;
    GaussPnt(6) = 0.3737060887154195;
    GaussPnt(7) = -0.5108670019508271;
    GaussPnt(8) = 0.5108670019508271;
    GaussPnt(9) = -0.636053680726515;
    GaussPnt(10) = 0.636053680726515;
    GaussPnt(11) = -0.7463319064601508;
    GaussPnt(12) = 0.7463319064601508;
    GaussPnt(13) = -0.8391169718222188;
    GaussPnt(14) = 0.8391169718222188;
    GaussPnt(15) = -0.9122344282513259;
    GaussPnt(16) = 0.9122344282513259;
    GaussPnt(17) = -0.9639719272779138;
    GaussPnt(18) = 0.9639719272779138;
    GaussPnt(19) = -0.9931285991850949;
    GaussPnt(20) = 0.9931285991850949;
    GaussWeight(1) = 0.1527533871307258;
    GaussWeight(2) = 0.1527533871307258;
    GaussWeight(3) = 0.1491729864726037;
    GaussWeight(4) = 0.1491729864726037;
    GaussWeight(5) = 0.142096109318382;
    GaussWeight(6) = 0.142096109318382;
    GaussWeight(7) = 0.1316886384491766;
    GaussWeight(8) = 0.1316886384491766;
    GaussWeight(9) = 0.1181945319615184;
    GaussWeight(10) = 0.1181945319615184;
    GaussWeight(11) = 0.1019301198172404;
    GaussWeight(12) = 0.1019301198172404;
    GaussWeight(13) = 0.0832767415767048;
    GaussWeight(14) = 0.0832767415767048;
    GaussWeight(15) = 0.0626720483341091;
    GaussWeight(16) = 0.0626720483341091;
    GaussWeight(17) = 0.0406014298003869;
    GaussWeight(18) = 0.0406014298003869;
    GaussWeight(19) = 0.0176140071391521;
    GaussWeight(20) = 0.0176140071391521;

   
end

end
