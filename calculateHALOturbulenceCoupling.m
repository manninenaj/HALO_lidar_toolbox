function [ tkefield ] = calculateHALOturbulenceCoupling(epsilon,skewn,cloudmask,th,sun_rise_set,time)
%calculateHALOturbulenceCoupling creates a mask indicating is the turbulence
%in contact with surface, cloud, both, or neither of them.
%
%0: no signal
%1: non turbulent
%2: connected with surface
%3: connectedt with cloud
%4: in cloud
%5: unconnected

% rm nans
skewn(isnan(epsilon)) = nan;

% Skewness field
%0: no signal
%1: pos skewness connected with surface
%2: neg skewness
%3: in cloud
%4: neg skewness connected with cloud
%5: pos skewness unconnected

% Initialize
Skewnfield = zeros(size(epsilon));
Skewnfield(skewn > th.vert_velo) = 1;
Skewnfield(skewn <= th.vert_velo) = 2;
Skewnfield(cloudmask) = 3; % cloud

% Find negative skewness connected with cloud
for ii = 1:size(Skewnfield,1)
    if any(cloudmask(ii,:),2) % cloud in profile?
        for jj = size(Skewnfield,2):-1:3
            % cloud top
            if Skewnfield(ii,jj)==3 % if cloud
                jj2 = jj;
                % cloud bottom
                while jj2>4 && Skewnfield(ii,jj2)==3
                    jj2 = jj2-1;
                end
                % if neg. skewness
                if jj2>4 && (Skewnfield(ii,jj2)==2 || Skewnfield(ii,jj2)==0) 
                    jj3=jj2;
                    % while neg. skewness
                    while jj3>3 && (Skewnfield(ii,jj3)==2 || Skewnfield(ii,jj3)==0) 
                        jj3 = jj3-1;
                        % if pos. but next range gate is negative, i.e.
                        % allow one pos. in between neg. gates
                        if Skewnfield(ii,jj3)~=2 && Skewnfield(ii,jj3-1)==2
                            jj3 = jj3-1;
                        end
                    end
                    % From cloud bottom to when skewness changes to pos.
                    Skewnfield(ii,jj2:-1:jj3+1)=4;
                end
                
            end
        end
    end
end

% % % Find negative skewness connected to surface
% % jj7 = 4;
% % for ii7 = 1:size(Skewnfield,1)
% %     if  Skewnfield(ii7,jj7) == 2
% %         while Skewnfield(ii7,jj7) == 2 && jj7 < size(Skewnfield,2)-1
% %             jj7 = jj7 + 1;
% %         end
% %         Skewnfield(ii7,4:jj7) = 6;
% %         jj7 = 4;
% %     end
% % end

% Find skewness not in contact with surface
for ii = 1:size(Skewnfield,1)
    izero = find(Skewnfield(ii,4:end)~=1,1,'first');
    ione = find(Skewnfield(ii,4:end)==1);
    if ~isempty(ione) && ~isempty(izero)
        ione(ione<izero) = []; % if '1' higher than '0', remove
        Skewnfield(ii,3+ione) = 5;
    end
end


%0: no signal
%1: non turbulent
%2: epsilon > 10^-5 & connected with surface
%3: connected with cloud
%4: in cloud
%5: unconnected
%6: epsilon > 10^⁻4 & connected to surface (convective)
tkefield = zeros(size(epsilon));
tkefield(~isnan(epsilon)) = 1;
tkefield(epsilon > th.epsilon) = 2; % hi turbulence
tkefield(epsilon > th.epsilon & Skewnfield == 4) = 3; % connected w/ cloud
tkefield(cloudmask) = 4; % in cloud
tkefield(tkefield == 3 & repmat(~any(cloudmask,2),1,size(tkefield,2))) = 5; % if cloud driven but no clouds in profile
% tkefield(:,1:3) = nan; % ignore

jj = 3;
for ii = 1:size(tkefield,1)
    if time(ii) >= sun_rise_set(1) && time(ii) <= sun_rise_set(2)
        switch any(cloudmask(ii,:))
            case 1 % cloud within boundary layer, use also skewness
                while not(isnan(epsilon(ii,jj))) && epsilon(ii,jj) > th.epsilon_hi && ...
                       Skewnfield(ii,jj) ~= 4 && Skewnfield(ii,jj) ~= 3 
                    jj = jj + 1;
                end
            case 0 % no clouds within boundary layer, use only epsilon
                while epsilon(ii,jj) > th.epsilon_hi
                    jj = jj + 1;
                end
        end
        if jj>4
            tkefield(ii,4:jj) = 6;
        end
        jj = 4;
    end
end

% Find epsilon not in contact with surface
for ii = 1:size(tkefield,1)
    izero = find(tkefield(ii,4:end)~=2 & tkefield(ii,4:end)~=6,1,'first');
    ione = find(tkefield(ii,4:end)==2 | tkefield(ii,4:end)==6);
    if ~isempty(ione)
        ione(ione<izero) = []; % if '1' higher than '0', remove
        tkefield(ii,3+ione) = 5;
    end
end

% Fix false "in contact with cloud" in between "non turbulent"
for ii = 1:size(tkefield,1)
    for jj = 3:size(tkefield,2)-1
        if  tkefield(ii,jj) == 3
            jj2 = jj + 1;
            while jj2 < size(tkefield,2) && tkefield(ii,jj2) == 3
                jj2 = jj2 + 1;
            end
            % if cloud not above cloud driven within 100 m
            if ~any(ismember([tkefield(ii,jj2),tkefield(ii,jj2+1),tkefield(ii,jj2+2),tkefield(ii,jj2+3)],4))
%                 if tkefield(ii,jj2) == 0
                    tkefield(ii,jj:jj2) = 5; % "not in contact with either"
%                 end
            end
        end
    end
end


tkefield(isnan(epsilon) & tkefield ~= 4) = 0;
tkefield(isnan(tkefield)) = 0; % no signal
tkefield(:,1:3) = nan;
end

