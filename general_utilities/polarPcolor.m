function [varargout] = polarPcolor(R,theta,Z,varargin)
% [h,c] = polarPcolor1(R,theta,Z,varargin) is a pseudocolor plot of matrix 
% Z for a vector radius R and a vector angle theta. 
% The elements of Z specify the color in each cell of the 
% plot. The goal is to apply pcolor function with a polar grid, which 
% provides a better visualization than a cartesian grid.
%
%% Syntax
% 
% [h,c] = polarPcolor(R,theta,Z)
% [h,c] = polarPcolor(R,theta,Z,'Ncircles',10)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5,'colBar',0) 
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5,'labelR','r (km)')
% 
% INPUT
%	* R :
%        - type: float
%        - size: [1 x Nrr ] where Nrr = numel(R).
%        - dimension: radial distance.
%	* theta : 
%        - type: float
%        - size: [1 x Ntheta ] where Ntheta = numel(theta).
%        - dimension: azimuth or elevation angle (deg).
%        - N.B.: The zero is defined with respect to the North.
%	* Z : 
%        - type: float
%        - size: [Ntheta x Nrr]
%        - dimension: user's defined .
%	* varargin:
%        - Ncircles: number  of circles for the grid definition.
%        - Nspokes: number of spokes for the grid definition.
%        - colBar: display the colorbar or not.
%        - labelR: title for radial axis.
%        - RtickLabel: Tick label for the radial axis.
% 
% 
% OUTPUT
% h: returns a handle to a SURFACE object.
% c: returns a handle to a COLORBAR object.
%
%% Examples 
% R = linspace(3,10,100);
% theta = linspace(0,180,360);
% Z = linspace(0,10,360)'*linspace(0,10,100);
% figure
% polarPcolor(R,theta,Z,'Ncircles',3)
%
%% Adapted from polarPcolor created by:
% Etienne Cheynet
% https://se.mathworks.com/matlabcentral/fileexchange/49040-pcolor-in-polar-coordinates

p.Ncircles = nan;
p.Origin = 'auto';
p.Nspokes = 8;
p.labelR = '';
p.RtickLabel = [];
p.colBar = 1;
p.Rscale = 'linear';
p.rMin = 0;
p.rMax = max(R);
p.thetaMin = 0;
p.thetaMax = 360;
p.rStep = nan;

% Was the 'parameters' struct supplied?
if ~isempty(varargin)
    % Check for overrides of the defaults
    p = parsePropertyValuePairs(p, varargin);
end

if strcmpi(p.Origin,'auto') && strcmpi(p.Rscale,'log'),
    warning(' The origin cannot be set to 0 if R is expressed on a logarithmic axis. The value ''Rmin'' is used instead')
    Origin = 'Rmin';
end
if ~isempty(p.RtickLabel),
    if numel(p.RtickLabel)~=p.Ncircles,
        error(' The radial ticklabel must be equal to Ncircles');
    end
    if any(cellfun(@ischar,p.RtickLabel)==0),
        error(' The radial ticklabel must be a cell array of characters');
    end
end
%% Preliminary checks
% case where dimension is reversed
Nrr = numel(R);
Noo = numel(theta);
if isequal(size(Z),[Noo,Nrr]) && Noo~=Nrr,
  Z=transpose(Z);
end
% case where dimension of Z is not compatible with theta and R
if ~isequal(size(Z),[Nrr,Noo])
    fprintf('\n')
    fprintf([ 'Size of Z is : [',num2str(size(Z)),'] \n']);
    fprintf([ 'Size of R is : [',num2str(size(R)),'] \n']);
    fprintf([ 'Size of theta is : [',num2str(size(theta)),'] \n\n']);
    error(' dimension of Z does not agree with dimension of R and Theta')
end
if ~strcmpi(p.Origin,'auto'),  % if the origin is not automatically centered at rMin, then origin is fixed at 0
    p.Origin = 'Rmin';
end

%% data plot    
% Definition of the mesh
cax = newplot;
p.Rrange = p.rMax - p.rMin; % get the range for the radius
[rNorm] = getRnorm(R,p); % getRnorm is a nested function
YY = (rNorm)'*cosd(theta);
XX = (rNorm)'*sind(theta);
h = pcolor(XX,YY,Z,'parent',cax);
% disp([max(R/p.Rrange),max(rNorm)])
shading flat
set(cax,'dataaspectratio',[1 1 1]);axis off;
if ~ishold(cax);
    % make a radial grid
    hold(cax,'on')
    % Draw circles and spokes
    createSpokes(p);
    createCircles(p)
end
%% PLot colorbar if specified
if p.colBar==1,
    c=colorbar('location','WestOutside');
    caxis([quantile(Z(:),0.01),quantile(Z(:),0.99)])

else
    c = [];
end
pause(.2)
axpos = get(gca,'Position');
if ~isempty(c)
    cpos = c.Position;
    cpos(3) = .01;
    cpos(1) = cpos(1) + .02; 
    c.Position = cpos;
end
pause(.2)
set(gca,'Position',axpos)
%% Outputs
nargoutchk(0,2)
if nargout==1,
    varargout{1}=h;
elseif nargout==2,
    varargout{1}=h;
    varargout{2}=c;
end

    %%--createSpokes--%%
    function createSpokes(p)
        
        spokeMesh = linspace(p.thetaMin,p.thetaMax,p.Nspokes);
        if ~isnan(p.rStep)
	    circleMesh = [0:p.rStep:p.rMax p.rMax];
        else
            circleMesh = linspace(p.rMin,p.rMax,p.Ncircles); 
        end
        contourD = abs((circleMesh - circleMesh(1))/p.Rrange+R(1)/p.Rrange);
        
        cost = cosd(90-spokeMesh); % the zero angle is aligned with North
        sint = sind(90-spokeMesh); % the zero angle is aligned with North
        for kk = 1:p.Nspokes
            
            X = cost(kk)*contourD;
            Y = sint(kk)*contourD;
            
            if strcmpi(p.Origin,'Rmin'),
                X(1)=0;
                Y(1)=0;
            end
            plot(X,Y,'color',[0.35,0.35,0.35],'linewidth',0.5,...
                'handlevisibility','off');
            % plot graduations of angles
            % avoid superimposition of 0 and 360
            if and(p.thetaMin==0,p.thetaMax == 360),
                if spokeMesh(kk)<360,
                    
                    text(1.15.*contourD(end).*cost(kk),...
                        1.15.*contourD(end).*sint(kk),...
                        [num2str(spokeMesh(kk),3),char(176)],...
                        'horiz', 'center', 'vert', 'middle');
                end
            else
                text(1.15.*contourD(end).*cost(kk),...
                    1.15.*contourD(end).*sint(kk),...
                    [num2str(spokeMesh(kk),3),char(176)],...
                    'horiz', 'center', 'vert', 'middle');
            end
            
        end
    end

    %%-- createCircles --%%
	    function createCircles(p)

        
        if ~isnan(p.rStep)
            contourD = [p.rMin:p.rStep:p.rMax p.rMax];
            contourD = contourD./p.rMax;
        else
            if strcmpi(p.Origin,'Rmin'),  % if the origin is set at rMin
                contourD = linspace(0,1+R(1)/p.Rrange,p.Ncircles);
            else % if the origin is automatically centered at 0
                contourD = linspace(0,1,p.Ncircles)+R(1)/p.Rrange;
            end
        end
        if strcmpi(p.Rscale,'linear')||strcmpi(p.Rscale,'lin')
           if ~isnan(p.rStep)
   	      tickMesh =  0:p.rStep:p.rMax;
           else
              tickMesh = linspace(p.rMin,p.rMax,p.Ncircles);
           end
        elseif strcmpi(p.Rscale,'log')||strcmpi(p.Rscale,'logarithmic'),
            tickMesh  = logspace(log10(p.rMin),log10(p.rMax),p.Ncircles);
        else
            error('''Rscale'' must be ''log'' or ''linear'' ');
        end

        
        % define the grid in polar coordinates
        angleGrid = linspace(90-p.thetaMin,90-p.thetaMax,100);
        xGrid = cosd(angleGrid);
        yGrid = sind(angleGrid);
        
        spokeMesh = linspace(p.thetaMin,p.thetaMax,p.Nspokes);
        
        % plot circles
        for kk=1:length(contourD)
            X = xGrid*contourD(kk);
            Y = yGrid*contourD(kk);
            plot(X,Y,'color',[0.35,0.35,0.35],'linewidth',.5);
        end
        % radius tick label
        for kk=1:length(tickMesh)
            
            position = 0.51.*(spokeMesh(min(p.Nspokes,round(p.Ncircles/2)))+...
                spokeMesh(min(p.Nspokes,1+round(p.Ncircles/2))));
            
            if isempty(p.RtickLabel),
                rtick = num2str(tickMesh(kk),2);
            else
                rtick = p.RtickLabel(kk);
            end
            if abs(round(position)) ==90,
                % radial graduations
                text((contourD(kk)).*cosd(90-position),...
                    (0.1+contourD(kk)).*sind(86-position),...
                    rtick,'verticalalignment','BaseLine',...
                    'horizontalAlignment', 'center',...
                    'handlevisibility','off','parent',cax);
                % annotate spokes
                text(contourD(end).*0.6.*cosd(90-position),...
                    0.07+contourD(end).*0.6.*sind(90-position),...
                    [labelR],'verticalalignment','bottom',...
                    'horizontalAlignment', 'right',...
                    'handlevisibility','off','parent',cax);
            else
                % radial graduations
               text((contourD(kk)).*cosd(90-position),...
                    (contourD(kk)).*sind(90-position),...
                    rtick,'verticalalignment','BaseLine',...
                    'horizontalAlignment', 'right',...
                    'handlevisibility','off','parent',cax);
                % annotate spokes
                text(contourD(end).*0.6.*cosd(90-position),...
                    contourD(end).*0.6.*sind(90-position),...
                    [p.labelR],'verticalalignment','bottom',...
                    'horizontalAlignment', 'right',...
                    'handlevisibility','off','parent',cax);
            end
        end
        
    end
    
    %%--getRNorm--%%
    function [rNorm] = getRnorm(R,p)
        if strcmpi(p.Rscale,'linear')||strcmpi(p.Rscale,'lin'),
            if strcmpi(p.Origin,'Rmin'),
                rNorm = R-R(1);
                rNorm = (rNorm)/max(rNorm)*max(R/p.Rrange);
            else
                rNorm = R/p.Rrange;
            end
        elseif strcmpi(p.Rscale,'log')||strcmpi(p.Rscale,'logarithmic'),
            if rMin<=0,
                error(' The radial vector cannot be lower or equal to 0 if the logarithmic scale is used');
            end
            rNorm = log10(R); %normalized radius [0,1]
            rNorm =rNorm-rNorm(1);
            rNorm = (rNorm)/max(rNorm)*max(R/p.Rrange);
        else
            error('''Rscale'' must be ''log'' or ''linear'' ');
        end
    end

    %%-- parsePropertyValuePairs --%%
    %DErrico, John (2006). Parsing property/value pairs for function input
    %(http://www.mathworks.com/matlabcentral/fileexchange/9082-parse-pv-pairs),
    %MATLAB Central File Exchange. Retrieved Oct 7, 2015.

    function params = parsePropertyValuePairs(params, pv_pairs)
        %PARSEPROPERTYVALUEPAIRS: parses sets of property value pairs
        % usage: params = parse_pv_pairs(default_parameters, pv_pairs)
        %
        % arguments: (input)
        %  default_parameters - structure, with one field for every
        %   potential property/value pair. Each field will contain the
        %   default value for that property. If no default is supplied for
        %   a given property, then that field must be empty.
        %
        %  pv_array - cell array of property/value pairs.
        %   Case is ignored when comparing properties to the list of field
        %   names. Also, any unambiguous shortening of a field/property
        %   name is allowed.
        %
        % arguments: (output)
        %  params - parameter struct that reflects any updated
        %  property/value pairs in the pv_array.
        %
        % Example usage:
        % First, set default values for the parameters. Assume we
        % have four parameters that we wish to use optionally in
        % the function examplefun.
        %
        %  - 'viscosity', which will have a default value of 1
        %  - 'volume', which will default to 1
        %  - 'pie' - which will have default value 3.141592653589793
        %  - 'description' - a text field, left empty by default
        %
        % The first argument to examplefun is one which will always be
        % supplied.
        %
        %   function examplefun(dummyarg1,varargin)
        %   params.Viscosity = 1;
        %   params.Volume = 1;
        %   params.Pie = 3.141592653589793
        %
        %   params.Description = '';
        %   params=parse_pv_pairs(params,varargin);
        %   params
        %
        % Use examplefun, overriding the defaults for 'pie', 'viscosity'
        % and 'description'. The 'volume' parameter is left at its default.
        %
        %   examplefun(rand(10),'vis',10,'pie',3,'Description',
        %    'Hello world')
        %
        % params =
        %     Viscosity: 10
        %        Volume: 1
        %           Pie: 3
        %   Description: 'Hello world'
        %
        % Note that capitalization was ignored, the property 'viscosity'
        % was truncated as supplied. Also, note that the order the pairs
        % were supplied was arbitrary.
        
        npv = length(pv_pairs);
        n = npv/2;
        
        if n~=floor(n)
            error 'Property/value pairs must come in PAIRS.'
        end
        if n<=0
            % just return the defaults
            return
        end
        
        if ~isstruct(params)
            error 'No structure for defaults was supplied'
        end
        
        % there was at least one pv pair. process any supplied
        propnames = fieldnames(params);
        lpropnames = lower(propnames);
        for i_sf=1:n
            p_i = lower(pv_pairs{2*i_sf-1});
            v_i = pv_pairs{2*i_sf};
            
            ind = strmatch(p_i,lpropnames,'exact');
            if isempty(ind)
                ind = find(strncmp(p_i,lpropnames,length(p_i)));
                if isempty(ind)
                    error(['No matching property found for: ',...
                        pv_pairs{2*i_sf-1}])
                elseif length(ind)>1
                    error(['Ambiguous property name: ',pv_pairs{2*i_sf-1}])
                end
            end
            p_i = propnames{ind};
            
            % override the corresponding default in params
            params = setfield(params,p_i,v_i); %#ok
            
        end
    end
end

