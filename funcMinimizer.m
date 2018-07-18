function [x,fval] = funcMinimizer(funfcn,x0)
%funcMinimizer tried to find a local minimizer 'x0' of function 'funfcn.'
% 'funfcn' is a function handle, and accepts input 'x0' and returns scalar 
% function value 'fval' evaluated at 'x', which can be scalar, vector or 
% matrix.

Params.max_func = 800;
Params.max_iterations = 800;
Params.tolerance_func = 1.0000e-4;
Params.tolerance_x = 1.0000e-4;
varargin = {};
Params.no_of_x = numel(x0);

% Initialize parameters
Params.rho = 1;
Params.chi = 2;
Params.psi = 0.5;
Params.sigma = 0.5;
fubar.ones_n_of_x = ones(1,Params.no_of_x);
fubar.two_2_n_p_1 = 2:Params.no_of_x+1;
fubar.one_2_n = 1:Params.no_of_x;

% Set first guesses.
fubar.x_in = x0(:); 
fubar.vee = zeros(Params.no_of_x,Params.no_of_x+1); fubar.fv = ...
    zeros(1,Params.no_of_x+1);
fubar.vee(:,1) = fubar.x_in; % put first guesses in
x(:) = fubar.x_in;    % Change x for funfcn
fubar.fv(:,1) = funfcn(x,varargin{:});
itercount = 0;

Params.usual_delta = 0.05; % 5 percent 
Params.zero_term_delta = 0.00025; % 0.025 percent for zero elements of x
for j = 1:Params.no_of_x
    fubar.y = fubar.x_in;
    if fubar.y(j) ~= 0
        fubar.y(j) = (1 + Params.usual_delta)*fubar.y(j);
    else
        fubar.y(j) = Params.zero_term_delta;
    end
    fubar.vee(:,j+1) = fubar.y;
    x(:) = fubar.y;
    fubar.f = funfcn(x,varargin{:});
    fubar.fv(1,j+1) = fubar.f;
end

% sort so vee(1,:) has the lowest function value
[fubar.fv,j] = sort(fubar.fv);
fubar.vee = fubar.vee(:,j);

Params.how = 'initial simplex';
itercount = itercount + 1;
Params.func_evals = Params.no_of_x+1;

while Params.func_evals < Params.max_func && itercount < ...
        Params.max_iterations
    if max(abs(fubar.fv(1)-fubar.fv(fubar.two_2_n_p_1))) <= ...
            max(Params.tolerance_func,10*eps(fubar.fv(1))) && ...
            max(max(abs(fubar.vee(:,fubar.two_2_n_p_1)-...
        fubar.vee(:,fubar.ones_n_of_x)))) <= ...
            max(Params.tolerance_x,10*eps(max(fubar.vee(:,1))))
        break
    end
    % Compute the reflection point
    % fubar.xbar = average of the n (NOT n+1) best points
    fubar.xbar = sum(fubar.vee(:,fubar.one_2_n), 2)/Params.no_of_x;
    fubar.xr = (1 + Params.rho)*fubar.xbar - Params.rho*fubar.vee(:,end);
    x(:) = fubar.xr;
    fubar.fxr = funfcn(x,varargin{:});
    Params.func_evals = Params.func_evals+1;
    if fubar.fxr < fubar.fv(:,1)
        % Calculate the expansion point
        fubar.xe = (1 + Params.rho*Params.chi)*fubar.xbar - Params.rho *...
            Params.chi * fubar.vee(:,end);
        x(:) = fubar.xe; fubar.fxe = funfcn(x,varargin{:});
        Params.func_evals = Params.func_evals+1;
        if fubar.fxe < fubar.fxr
            fubar.vee(:,end) = fubar.xe;
            fubar.fv(:,end) = fubar.fxe;
            Params.how = 'expand';
        else
            fubar.vee(:,end) = fubar.xr;
            fubar.fv(:,end) = fubar.fxr;
            Params.how = 'reflect';
        end
    else % fv(:,1) <= fubar.fxr
        if fubar.fxr < fubar.fv(:,Params.no_of_x)
            fubar.vee(:,end) = fubar.xr;
            fubar.fv(:,end) = fubar.fxr;
            Params.how = 'reflect';
        else % fubar.fxr >= fv(:,n)
            % Perform contraction
            if fubar.fxr < fubar.fv(:,end)
                % Perform an outside contraction
                fubar.xc = (1 + Params.psi * Params.rho) * ...
                    fubar.xbar - Params.psi * Params.rho * ...
                    fubar.vee(:,end);
                x(:) = fubar.xc;
                fubar.fxc = funfcn(x,varargin{:});
                Params.func_evals = Params.func_evals + 1;
                
                if fubar.fxc <= fubar.fxr
                    fubar.vee(:,end) = fubar.xc;
                    fubar.fv(:,end) = fubar.fxc;
                    Params.how = 'contract outside';
                else
                    % perform a shrink
                    Params.how = 'shrink';
                end
            else
                % Perform an inside contraction
                fubar.xcc = (1-Params.psi) * fubar.xbar + ...
                    Params.psi * fubar.vee(:,end);
                x(:) = fubar.xcc; fubar.fxcc = funfcn(x,varargin{:});
                Params.func_evals = Params.func_evals + 1;
                
                if fubar.fxcc < fubar.fv(:,end)
                    fubar.vee(:,end) = fubar.xcc;
                    fubar.fv(:,end) = fubar.fxcc;
                    Params.how = 'contract inside';
                else
                    % perform a shrink
                    Params.how = 'shrink';
                end
            end
            if strcmp(Params.how,'shrink')
                for j = fubar.two_2_n_p_1
                    fubar.vee(:,j) = fubar.vee(:,1) + Params.sigma * ...
                        (fubar.vee(:,j) - fubar.vee(:,1));
                    x(:) = fubar.vee(:,j); 
                    fubar.fv(:,j) = funfcn(x,varargin{:});
                end
                Params.func_evals = Params.func_evals + Params.no_of_x;
            end
        end
    end
    [fubar.fv,j] = sort(fubar.fv);
    fubar.vee = fubar.vee(:,j);
    itercount = itercount + 1;
end   % while

x(:) = fubar.vee(:,1);
fval = fubar.fv(:,1);
end
