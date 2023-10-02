function [soln,solnw,eqn,info,c_coef] = Poisson(node,elem,bdFlag,pde,option,usingint,tau,R,omega)
%% POISSON Poisson equation: P1 linear element.
%
%   u = Poisson(node,elem,bdFlag,pde) produces the linear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
% 
%   The mesh is given by node and elem and the boundary edge is by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%   
%   soln = Poisson(node,elem,bdFlag,pde,option) specifies the options.
%
%   In the output, soln structure contains
%     - soln.u: solution u
%     - soln.Du: gradient of u
%
%   In the input, option structures contains
%    - option.dquadorder: quadrature order for diffusion coefficients
%    - option.fquadorder: quadrature order for computing right hand side f
%    - option.solver
%      'direct': the built in direct solver \ (mldivide)
%      'mg':     multigrid-type solvers mg is used.
%      'amg':    algebraic multigrid method is used.
%      'none':   only assemble the matrix equation but not solve. 
%   The default setting is to use the direct solver for small size problems
%   and multigrid solvers for large size problems. For more options on the
%   multigrid solver mg, type help mg.
%
%   The function Poisson assembes the matrix equation AD*u = b and solves
%   it by the direct solver (small size <= 2e3) or the multigrid solver
%   (large size > 2e3). The Dirichlet boundary condition is built into the
%   matrix AD and the Neumann boundary condition is build into b.
%
%   The diffusion coefficient d is a scalar function or a column array with
%   the same length as the elem array. 
%
%   When only one type of boundary condition is imposed, the input argument
%   bdFlag can be skipped. The boundary condition is implicitly given in
%   the pde structure by specifying g_D or g_N only. See examples below.
%
%   [soln,eqn] = Poisson(node,elem,bdFlag,pde) returns also the equation
%   structure eqn, which includes: 
%     - eqn.AD:  modified stiffness matrix AD;
%     - eqn.b:   the right hand side. 
%     - eqn.Lap: non-modified stiffness matrix
%
%   The solution u = AD\b. The matrix eqn.Lap can be used to evulate the
%   bilinear form a(u,v) = u'*eqn.Lap*v, especially the enery norm of a finite
%   element function u is given by by sqrt(u'*eqn.Lap*u). 
%
%   [soln,eqn,info] = Poisson(node,elem,bdFlag,pde) returns also the
%   information on the assembeling and solver, which includes:
%     - info.assembleTime: time to assemble the matrix equation
%     - info.solverTime:   time to solve the matrix equation
%     - info.itStep:       number of iteration steps for the mg solver
%     - info.error:        l2 norm of the residual b - A*u
%     - info.flag:         flag for the mg solver.
%       flag = 0: converge within max iteration 
%       flag = 1: iterated maxIt times but did not converge
%       flag = 2: direct solver
%       flag = 3: no solve
%
%   Example
%     squarePoisson; Poissonfemrate;
%
%   Example
%     clear all
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:4
%       [node,elem] = uniformrefine(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = Poisson(node,elem,[],pde);
%     figure(1); 
%     showresult(node,elem,u);
%
%   See also Poisson3, femPoisson, squarePoisson, Poissonfemrate
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Preprocess
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
% important constants
N = size(node,1); 
NT = size(elem,1);
Ndof = N;
time = cputime;  % record assembling time

% [node1,elem1,ndFlag1,bdFlag1,HB1] = gradedrefine(node,elem,ndFlag,bdFlag);
% N1 = size(node1,1);  NT1 = size(elem1,1); 
% Ndof1 = N1;
% 
% [node2,elem2,ndFlag2,bdFlag2,HB2] = gradedrefine(node1,elem1,ndFlag1,bdFlag1);
% [node3,elem3,ndFlag3,bdFlag3,HB3] = gradedrefine(node2,elem2,ndFlag2,bdFlag2);

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);
for i = 1:3
    for j = i:3
        % $A_{ij}|_{\tau} = \int_{\tau}K\nabla \phi_i\cdot \nabla \phi_j dxdy$ 
        Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Aij; Aij],Ndof,Ndof);        
        end        
    end
end
clear K Aij

%% Assemble the right hand side
b = zeros(Ndof,1);
if isreal(pde.f) % f is a real number or vector and not a function
   switch length(pde.f)
       case NT  % f is piecewise constant
         bt = pde.f.*area/3;
         b = accumarray(elem(:),[bt; bt; bt],[Ndof 1]);
       case N   % f is piecewise linear
         bt = zeros(NT,3);
         bt(:,1) = area.*(2*pde.f(elem(:,1)) + pde.f(elem(:,2)) + pde.f(elem(:,3)))/12;
         bt(:,2) = area.*(2*pde.f(elem(:,2)) + pde.f(elem(:,3)) + pde.f(elem(:,1)))/12;
         bt(:,3) = area.*(2*pde.f(elem(:,3)) + pde.f(elem(:,1)) + pde.f(elem(:,2)))/12;
         b = accumarray(elem(:),bt(:),[Ndof 1]);
       case 1   % f is a scalar e.g. f = 1
         bt = pde.f*area/3;
         b = accumarray(elem(:),[bt; bt; bt],[Ndof 1]);
   end
end
if ~isempty(pde.f) && ~isreal(pde.f)  % f is a function 
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[Ndof 1]);
end
clear pxy bt

bw = b;

%% Set up boundary conditions
[AD,bw,w,freeNode,isPureNeumann] = getbd(node,elem,bdFlag,pde,A,bw);

%% Solve the system of linear equations for w
[w, info] = P1solver(AD,bw,w,freeNode,Ndof,elem,option);

%% Assemble the right hand side \Delta s^-
b = zeros(Ndof,1);
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = lapsingular(pxy,tau,R,omega);
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[Ndof 1]);

clear pxy bt

bs = b;

%% Set up boundary conditions
[AD,bs,s,freeNode,isPureNeumann] = getbd(node,elem,bdFlag,pde,A,bs);

%% Solve the system of linear equations for s
[s, info] = P1solver(AD,bs,s,freeNode,Ndof,elem,option);


%% Find the singular triangles
% selem = find(ndFlag(:,2)==1);
% nelemloc1 = find(elem(:,1)==selem | elem(:,2)==selem | elem(:,3)==selem);
% nelemloc = find(elem(:,1)==selem);
% nsingular = length(nelemloc)*0; % 1 means singular points realted integral will be based on 4 times points, 0 means regular points.
% if length(nelemloc1)~=length(nelemloc)
%     disp('singular point need to be the first point of the triangle');
%     stop
% end
% level = 1;

% tri0 = setdiff([1:NT]',[tri1; tri2; tri3]);
% level = [1:3]';

%% compute the coefficient c_h
% wxinorm=0;
% xi2norm=0;
% 
% for ek=1:length(tri3)
%     e = tri3(ek);
%     nodes=elem(e,:);
%     Pe=node(nodes,:);
%     
%             %[newlambda,Area,Arearatio] = subint(node,elem,node1,elem1,option.fquadorder,e);
%             %[newlambda,Area,Arearatio] = subintP1(node,elem,node1,elem1,node2,elem2,option.fquadorder,e);
%             [newlambda,Area,Arearatio] = subintP2(node,elem,node1,elem1,node2,elem2,node3,elem3,option.fquadorder,e,level(3));
% 
%             newphi3 = newlambda;
%             Ntri = length(Area);
%             xi2local = 0;
%             wxilocal = 0;
%             pxy = zeros(Ntri,2);
%             for pd = 1:nQuad
% 		    % quadrature points in the x-y coordinate
%             for ie = 1:Ntri
%                 pxy(ie,:) = newlambda(pd,1,ie)*Pe(1,:) ...
%             	 + newlambda(pd,2,ie)*Pe(2,:) ...
%             	 + newlambda(pd,3,ie)*Pe(3,:);
%             end
% 
% %             pxy = newlambda(pd,1,:).*Pe(1,:) ...
% %             	 + newlambda(pd,2,:).*Pe(2,:) ...
% %             	 + newlambda(pd,3,:).*Pe(3,:)
% 
%             
%             sj = usingular(pxy,tau,R,omega);
% 
%             W0 = zeros(Ntri,1);
%             Wj = zeros(Ntri,1);
%             for ei = 1:Ntri
%                 for ej = 1:3
%                     W0(ei) = W0(ei) + w(nodes(ej))*newphi3(pd,ej,ei);
%                     Wj(ei) = Wj(ei) + s(nodes(ej))*newphi3(pd,ej,ei);
%                 end
%                 
%                 xi2local=xi2local+Arearatio(ei)*weight(pd)*(Wj(ei)*Wj(ei)+2*Wj(ei)*sj(ei)+0*sj(ei)*sj(ei));
%                 wxilocal=wxilocal+Arearatio(ei)*weight(pd)*(W0(ei)*Wj(ei)+W0(ei)*sj(ei));     
%                 
%             end
% 
%             end
% 
%             
%             xi2norm = xi2norm + xi2local*area(e);
%             wxinorm = wxinorm + wxilocal*area(e);
%             
% end
% 
% for ek=1:length(tri2)
%     e = tri2(ek);
%     nodes=elem(e,:);
%     Pe=node(nodes,:);
%     
%             %[newlambda,Area,Arearatio] = subint(node,elem,node1,elem1,option.fquadorder,e);
%             %[newlambda,Area,Arearatio] = subintP1(node,elem,node1,elem1,node2,elem2,option.fquadorder,e);
%             [newlambda,Area,Arearatio] = subintP2(node,elem,node1,elem1,node2,elem2,node3,elem3,option.fquadorder,e,level(2));
% 
%             newphi2 = newlambda;
%             Ntri = length(Area);
%             xi2local = 0;
%             wxilocal = 0;
%             pxy = zeros(Ntri,2);
%             for pd = 1:nQuad
% 		    % quadrature points in the x-y coordinate
%             for ie = 1:Ntri
%                 pxy(ie,:) = newlambda(pd,1,ie)*Pe(1,:) ...
%             	 + newlambda(pd,2,ie)*Pe(2,:) ...
%             	 + newlambda(pd,3,ie)*Pe(3,:);
%             end
% 
% %             pxy = newlambda(pd,1,:).*Pe(1,:) ...
% %             	 + newlambda(pd,2,:).*Pe(2,:) ...
% %             	 + newlambda(pd,3,:).*Pe(3,:)
% 
%             
%             sj = usingular(pxy,tau,R,omega);
% 
%             W0 = zeros(Ntri,1);
%             Wj = zeros(Ntri,1);
%             for ei = 1:Ntri
%                 for ej = 1:3
%                     W0(ei) = W0(ei) + w(nodes(ej))*newphi2(pd,ej,ei);
%                     Wj(ei) = Wj(ei) + s(nodes(ej))*newphi2(pd,ej,ei);
%                 end
%                 
%                 xi2local=xi2local+Arearatio(ei)*weight(pd)*(Wj(ei)*Wj(ei)+2*Wj(ei)*sj(ei)+0*sj(ei)*sj(ei));
%                 wxilocal=wxilocal+Arearatio(ei)*weight(pd)*(W0(ei)*Wj(ei)+W0(ei)*sj(ei));     
%                 
%             end
% 
%             end
% 
%             
%             xi2norm = xi2norm + xi2local*area(e);
%             wxinorm = wxinorm + wxilocal*area(e);
%             
% end
% 
% for ek=1:length(tri1)
%     e = tri1(ek);
%     nodes=elem(e,:);
%     Pe=node(nodes,:);
%     
%             %[newlambda,Area,Arearatio] = subint(node,elem,node1,elem1,option.fquadorder,e);
%             %[newlambda,Area,Arearatio] = subintP1(node,elem,node1,elem1,node2,elem2,option.fquadorder,e);
%             [newlambda,Area,Arearatio] = subintP2(node,elem,node1,elem1,node2,elem2,node3,elem3,option.fquadorder,e,level(1));
% 
%             newphi1 = newlambda;
%             Ntri = length(Area);
%             xi2local = 0;
%             wxilocal = 0;
%             pxy = zeros(Ntri,2);
%             for pd = 1:nQuad
% 		    % quadrature points in the x-y coordinate
%             for ie = 1:Ntri
%                 pxy(ie,:) = newlambda(pd,1,ie)*Pe(1,:) ...
%             	 + newlambda(pd,2,ie)*Pe(2,:) ...
%             	 + newlambda(pd,3,ie)*Pe(3,:);
%             end
% 
% %             pxy = newlambda(pd,1,:).*Pe(1,:) ...
% %             	 + newlambda(pd,2,:).*Pe(2,:) ...
% %             	 + newlambda(pd,3,:).*Pe(3,:)
% 
%             
%             sj = usingular(pxy,tau,R,omega);
% 
%             W0 = zeros(Ntri,1);
%             Wj = zeros(Ntri,1);
%             for ei = 1:Ntri
%                 for ej = 1:3
%                     W0(ei) = W0(ei) + w(nodes(ej))*newphi1(pd,ej,ei);
%                     Wj(ei) = Wj(ei) + s(nodes(ej))*newphi1(pd,ej,ei);
%                 end
%                 
%                 xi2local=xi2local+Arearatio(ei)*weight(pd)*(Wj(ei)*Wj(ei)+2*Wj(ei)*sj(ei)+0*sj(ei)*sj(ei));
%                 wxilocal=wxilocal+Arearatio(ei)*weight(pd)*(W0(ei)*Wj(ei)+W0(ei)*sj(ei));     
%                 
%             end
% 
%             end
% 
%             
%             xi2norm = xi2norm + xi2local*area(e);
%             wxinorm = wxinorm + wxilocal*area(e);
%             
% end
% tic
% for ek = 1:length(tri0)
%     e = tri0(ek);
%             
%             nodes=elem(e,:);
%             Pe=node(nodes,:);
% 
%             for pd = 1:nQuad
%                 % quadrature points in the x-y coordinate
%                 pxy = lambda(pd,1)*Pe(1,:) ...
%                     + lambda(pd,2)*Pe(2,:) ...
%                     + lambda(pd,3)*Pe(3,:);
%                 sj = usingular(pxy,tau,R,omega);
%         
%                 W0 = w(nodes(1))*phi(pd,1)+w(nodes(2))*phi(pd,2)+w(nodes(3))*phi(pd,3);
%                 Wj = s(nodes(1))*phi(pd,1)+s(nodes(2))*phi(pd,2)+s(nodes(3))*phi(pd,3);
%         
% %                xi2norm=xi2norm+area(e)*weight(pd)*(Wj*Wj+2*Wj*sj+0*sj*sj);
%                 xi2norm=xi2norm+area(e)*weight(pd)*(Wj*Wj+2*Wj*sj);
%                 wxinorm=wxinorm+area(e)*weight(pd)*(W0*Wj+W0*sj);     
% 
%             end  
% 
% end
% toc

%% compute the coefficient c_h
wxinormloc=zeros(NT,1);
xi2normloc=zeros(NT,1);

    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);

    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		sj = usingular(pxy,tau,R,omega);
        W0 = zeros(NT,1);
        Wj = zeros(NT,1);
        for i = 1:3
            W0 = W0 + phi(p,i)*w(elem(:,i));
            Wj = Wj + phi(p,i)*s(elem(:,i));
        end
        
        xi2normloc = xi2normloc + weight(p)*( Wj.*Wj+2*Wj.*sj );
        wxinormloc = wxinormloc + weight(p)*( W0.*Wj+W0.*sj );
    end

    xi2normloc = xi2normloc.*area;
    wxinormloc = wxinormloc.*area;
    
    xi2norm = sum(xi2normloc);
    wxinorm = sum(wxinormloc);

    xi2norm = xi2norm + usingint;

c_coef = wxinorm/xi2norm;
c_wxi_xixi = [c_coef wxinorm xi2norm]
%c_coefref = 0.477003329105857

%% Assemble the right hand side w_h-c_h xi_h
b = zeros(Ndof,1);
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,3);
    
    
%     for ek = 1:length(tri3)
%         e=tri3(ek);
%     
%         %[newlambda,Area,Arearatio] = subint(node,elem,node1,elem1,option.fquadorder,e);
%         %[newlambda,Area,Arearatio] = subintP1(node,elem,node1,elem1,node2,elem2,option.fquadorder,e);
%         [newlambda,Area,Arearatio] = subintP2(node,elem,node1,elem1,node2,elem2,node3,elem3,option.fquadorder,e,level(3));
%         
%     	newphi3 = newlambda;
%         Ntri = length(Area);
%         
%             for pd = 1:nQuad
%                 pxy = zeros(Ntri,2);
% 		    % quadrature points in the x-y coordinate
%             for ie = 1:Ntri
%                 pxy(ie,:) = newlambda(pd,1,ie)*Pe(1,:) ...
%             	 + newlambda(pd,2,ie)*Pe(2,:) ...
%             	 + newlambda(pd,3,ie)*Pe(3,:);
%             end
%             
% %             pxy = newlambda(pd,1,:).*Pe(1,:) ...
% %             	 + newlambda(pd,2,:).*Pe(2,:) ...
% %             	 + newlambda(pd,3,:).*Pe(3,:)
% 
%             
%             sj = usingular(pxy,tau,R,omega);
% 
%             W0 = zeros(Ntri,1);
%             Wj = zeros(Ntri,1);
%             for ei = 1:Ntri
%                 for ej = 1:3
%                     W0(ei) = W0(ei) + w(nodes(ej))*newphi3(pd,ej,ei);
%                     Wj(ei) = Wj(ei) + s(nodes(ej))*newphi3(pd,ej,ei);
%                 end  
%                 
%         for i = 1:3
%             bt(e,i) = bt(e,i) + weight(pd)*Arearatio(ei)*newphi3(pd,i,ei)*(W0(ei)-1*c_coef*(Wj(ei)+sj(ei)));
%         end                
%                 
%             end
% 
%             end 
%         
%         
%     end
%     
%     for ek = 1:length(tri2)
%         e = tri2(ek);
%         %[newlambda,Area,Arearatio] = subint(node,elem,node1,elem1,option.fquadorder,e);
%         %[newlambda,Area,Arearatio] = subintP1(node,elem,node1,elem1,node2,elem2,option.fquadorder,e);
%         [newlambda,Area,Arearatio] = subintP2(node,elem,node1,elem1,node2,elem2,node3,elem3,option.fquadorder,e,level(2));
%         
%     	newphi2 = newlambda;
%         Ntri = length(Area);
%         
%             for pd = 1:nQuad
%                 pxy = zeros(Ntri,2);
% 		    % quadrature points in the x-y coordinate
%             for ie = 1:Ntri
%                 pxy(ie,:) = newlambda(pd,1,ie)*Pe(1,:) ...
%             	 + newlambda(pd,2,ie)*Pe(2,:) ...
%             	 + newlambda(pd,3,ie)*Pe(3,:);
%             end
%             
% %             pxy = newlambda(pd,1,:).*Pe(1,:) ...
% %             	 + newlambda(pd,2,:).*Pe(2,:) ...
% %             	 + newlambda(pd,3,:).*Pe(3,:)
% 
%             
%             sj = usingular(pxy,tau,R,omega);
% 
%             W0 = zeros(Ntri,1);
%             Wj = zeros(Ntri,1);
%             for ei = 1:Ntri
%                 for ej = 1:3
%                     W0(ei) = W0(ei) + w(nodes(ej))*newphi2(pd,ej,ei);
%                     Wj(ei) = Wj(ei) + s(nodes(ej))*newphi2(pd,ej,ei);
%                 end  
%                 
%         for i = 1:3
%             bt(e,i) = bt(e,i) + weight(pd)*Arearatio(ei)*newphi2(pd,i,ei)*(W0(ei)-1*c_coef*(Wj(ei)+sj(ei)));
%         end                
%                 
%             end
% 
%             end 
%         
%         
%     end
%     
%         for ek = 1:length(tri1)
%             e = tri1(ek);
%         %[newlambda,Area,Arearatio] = subint(node,elem,node1,elem1,option.fquadorder,e);
%         %[newlambda,Area,Arearatio] = subintP1(node,elem,node1,elem1,node2,elem2,option.fquadorder,e);
%         [newlambda,Area,Arearatio] = subintP2(node,elem,node1,elem1,node2,elem2,node3,elem3,option.fquadorder,e,level(1));
%         
%     	newphi1 = newlambda;
%         Ntri = length(Area);
%         
%             for pd = 1:nQuad
%                 pxy = zeros(Ntri,2);
% 		    % quadrature points in the x-y coordinate
%             for ie = 1:Ntri
%                 pxy(ie,:) = newlambda(pd,1,ie)*Pe(1,:) ...
%             	 + newlambda(pd,2,ie)*Pe(2,:) ...
%             	 + newlambda(pd,3,ie)*Pe(3,:);
%             end
%             
% %             pxy = newlambda(pd,1,:).*Pe(1,:) ...
% %             	 + newlambda(pd,2,:).*Pe(2,:) ...
% %             	 + newlambda(pd,3,:).*Pe(3,:)
% 
%             
%             sj = usingular(pxy,tau,R,omega);
% 
%             W0 = zeros(Ntri,1);
%             Wj = zeros(Ntri,1);
%             for ei = 1:Ntri
%                 for ej = 1:3
%                     W0(ei) = W0(ei) + w(nodes(ej))*newphi1(pd,ej,ei);
%                     Wj(ei) = Wj(ei) + s(nodes(ej))*newphi1(pd,ej,ei);
%                 end  
%                 
%         for i = 1:3
%             bt(e,i) = bt(e,i) + weight(pd)*Arearatio(ei)*newphi1(pd,i,ei)*(W0(ei)-1*c_coef*(Wj(ei)+sj(ei)));
%         end                
%                 
%             end
% 
%             end 
%         
%         
%     end
    

    
%    sj = zeros(NT,1);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		sj = usingular(pxy,tau,R,omega);
        W0 = zeros(NT,1);
        Wj = zeros(NT,1);
        for ei = 1:3
            W0 = W0 + w(elem(:,ei))*phi(p,ei);
            Wj = Wj + s(elem(:,ei))*phi(p,ei);
        end
%        Fe(ii)=Fe(ii)+Area(e)*weight(pd)*phi(pd,ii)*(W0-1*c_coef*(Wj+sj));   

        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*(W0-1*c_coef*(Wj+sj));
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[Ndof 1]);

clear pxy bt

bu = b;

%% Set up boundary conditions
[AD,bu,u,freeNode,isPureNeumann] = getbd(node,elem,bdFlag,pde,A,bu);

%% Solve the system of linear equations for s

[u, info] = P1solver(AD,bu,u,freeNode,Ndof,elem,option);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Compute Du
dudx =  u(elem(:,1)).*Dphi(:,1,1) + u(elem(:,2)).*Dphi(:,1,2) ...
      + u(elem(:,3)).*Dphi(:,1,3);
dudy =  u(elem(:,1)).*Dphi(:,2,1) + u(elem(:,2)).*Dphi(:,2,2) ...
      + u(elem(:,3)).*Dphi(:,2,3);         
Du = [dudx, dudy];

%% Compute Dw
dwdx =  w(elem(:,1)).*Dphi(:,1,1) + w(elem(:,2)).*Dphi(:,1,2) ...
      + w(elem(:,3)).*Dphi(:,1,3);
dwdy =  w(elem(:,1)).*Dphi(:,2,1) + w(elem(:,2)).*Dphi(:,2,2) ...
      + w(elem(:,3)).*Dphi(:,2,3);         
Dw = [dwdx, dwdy];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'freeNode',freeNode,'Lap',A);
    info.assembleTime = assembleTime;
end

%% Output
if nargout == 1
    solnw = w;
else
    solnw = struct('u',w,'Du',Dw);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of Poisson
