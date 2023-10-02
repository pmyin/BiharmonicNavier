function [err,time,solver,eqn] = femPoisson(mesh,pde,option,varargin)
%% FEMPOISSON solve Poisson equation by various finite element methods
%
%   FEMPOISSON computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% See also Poisson, crack, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

node = mesh.node;          % nodes
elem = mesh.elem;          % triangles
bdFlag = mesh.bdFlag;

%[node,elem] = uniformrefine(node,elem);
%[node,elem,ndFlag,bdFlag,HB] = gradedrefine(node,elem,ndFlag,bdFlag);
%figure;
%showmesh(node,elem);
%findnode(node);
%findelem(node,elem,'all','index','FaceColor',[0.5 0.9 0.45]);
% 
% %[elem2dof,edge,bdDof] = dofP2(elem);
% findedgedof(node,edge);
% 
% display(elem2dof);
% display(bdDof);
% 
% stop

% showmesh(node,elem);
% findnode(node);
% findelem(node,elem,'all','index','FaceColor',[0.5 0.9 0.45]);

%% Parameters
nv = size(elem,2);     
dim = size(node,2);
option = femoption(option);
elemType = option.elemType;     
refType = option.refType;
maxIt = option.maxIt;       
maxN = option.maxN;
L0 = option.L0;

%% Generate an initial mesh by refine L0 times
% for k = 1:L0
%     if strcmp(refType,'red')
%         if nv == 4
%             [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
%         else
%             [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%         end        
%     elseif strcmp(refType,'bisect')
%         [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
%     end
% end


node0 = node;          % nodes
elem0 = elem;          % triangles
bdFlag0 = bdFlag;

%% Initialize err
errL2u = zeros(maxIt,1);
errL2w = zeros(maxIt,1);  
errH1u = zeros(maxIt,1);
errH1w = zeros(maxIt,1); 
rateL2u = zeros(maxIt-2,1); 
rateL2w = zeros(maxIt-2,1);
errC = zeros(maxIt,1); 
rateC = zeros(maxIt-2,1);
rateH1u = zeros(maxIt-2,1); rateH1w = zeros(maxIt-2,1); 
errLinfu = zeros(maxIt,1); errLinfw = zeros(maxIt,1);
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1); h = zeros(maxIt,1);

tau=1/8;
R=0.9; %half of R in the paper.
%omega = 2*pi-getangle(node(1,:),node(9,:),node(2,:))
omega = 1.5*pi;
rate = pi/omega

usingint = getsint(omega,tau,2*R,option)

%% Finite Element Method        
for k = 1:maxIt
    k
    % refine mesh
    t = cputime;
    
%    if strcmp(refType,'red')
%        if nv == 4
%            [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
%        else
%            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%        end        
%    elseif strcmp(refType,'bisect')
%        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
%    end

    NT=size(elem,1);
%    [node,elem,ndFlag,bdFlag,HB] = gradedrefine(node,elem,ndFlag,bdFlag);


if k==0
    showmesh(node,elem)
    a=0.25;
    b=0.3;
    snode=[b a; 1-b a];
    findnode(snode);
    findnode(node);
    title('Mesh after 4 mesh refinements with \kappa=0.2')
end

%     showmesh(node,elem);
%     findnode(node);
%     findelem(node,elem,'all','index','FaceColor',[0.5 0.9 0.45]);
%     stop

%    dl=HBdl(HB,dl);
%    dlelem = getdlelem(dlelem,Nt);

%% Find the singular triangles
% selem = find(ndFlag(:,2)==1);
% nelemloc1 = find(elem(:,1)==selem | elem(:,2)==selem | elem(:,3)==selem);
% nelemloc = find(elem(:,1)==selem);
% nsingular = length(nelemloc)*1; % 1 means singular points realted integral will be based on 4 times points, 0 means regular points.
% if length(nelemloc1)~=length(nelemloc)
%     disp('singular point need to be the first point of the triangle');
%     stop
% end
% tri3 = [1:nsingular]';
% 
% if k==1
%     tri2 = [];
%     tri1 = [];
% end
% 
% if k==2
%     tri1 = [];
%     tri2 = getri2(node0,elem0,ndFlag0,bdFlag0,tri3,k-1);
%     
% end
% 
% if k>2
%     tri2 = getri2(node0,elem0,ndFlag0,bdFlag0,tri3,k-1);
%     tri1 = getri1(node0,elem0,ndFlag0,bdFlag0,tri3,tri2,k-2);
% end
% 
% tri0 = [1:NT]';
% showmesh(node,elem(tri1,:));
% findnode(node);
% %findelem(node,elem,'all','index','FaceColor',[0.5 0.9 0.45]);
% title('Initial mesh for graded mesh')
% if k==3
% stop
% end

% tri1
% tri2
% tri3
% NT = size(elem,1)
% 
% if k==3
%     stop
% end

%    writemesh(node,elem,k);

    meshTime(k) = cputime - t;
    tic
    switch elemType
        case 'P1'     % piecewise linear function P1 element

            [soln,solnw,eqn,info,c_coef] = Poisson(node,elem,bdFlag,pde,option,usingint,tau,R,omega);
                  
         
    end

    toc
    
    if k<maxIt
        [node,elem,bdFlag,HB] = uniformrefine(node,elem,bdFlag);
    end
    
end

% b=1/cos(pi/8);
% a=tan(pi/8)+b;
figure
showsolution(node,elem,soln.u,2);
title('P^1 approximation after 7 mesh refinements')
%axis(2*[-a/sqrt(2)-sqrt(1/2) b -b a/sqrt(2)+sqrt(1/2)]);
axis(2*[-1 1 -1 1]);
colorbar;
grid off;

