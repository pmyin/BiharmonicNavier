function [u, info] = P1solver(AD,b,u,freeNode,Ndof,elem,option)

%% Solve the system of linear equations for u
if isempty(freeNode), info.solverTime=0; info.itStep=0;info.stopErr=0;info.flag=0; return; end
% Set up solver type

if isempty(option) || ~isfield(option,'solver')  || isfield(option,'mgoption')   % no option.solver
    if Ndof <= 2e10  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        t = cputime;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'mg'
        if ~isfield(option,'mgoption')   % no option.mgoption
            option.mgoption.x0 = u;
            option.mgoption.solver = 'CG';
        end
        [u,info] = mg(AD,b,elem,option.mgoption);
    case 'amg'
        if ~isfield(option,'amgoption')  % no option.amgoption
            option.amgoption.x0 = u;
            option.amgoption.solver = 'CG';
        end
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option.amgoption);                 
end