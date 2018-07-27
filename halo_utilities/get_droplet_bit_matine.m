function droplet_bit = get_droplet_bit_matine(height, beta, temperature)

% get_droplet_bit_matine.m is modified from original Cloudnet function
% get_droplet_bit.m
% updated by AH since Feb 2017 (agreed with EOC)

beta0 = beta;
beta0(find(~isfinite(beta0))) = 0;

droplet_bit = zeros(size(beta));
dheight = median(diff(height));
nrays = size(beta,2);
ngates = length(height);
% size(beta)
threshold_beta = 2e-5;
min_jump = 10; % Factor of 20
jump_distance = 250; % m
jump_gates = ceil(jump_distance / dheight);
final_gate = ngates-jump_gates;
search_gates_below = ceil(100/dheight);
search_gates_above = ceil(300/dheight);
diff_factor = 0.25;
min_temperature = 273.16 - 40; % No droplets colder than this

disp('Calculating location of liquid water cloud droplets')
for ii = 1:nrays
    start_gate = 2;
    profile = beta0(:,ii);
    while start_gate <= final_gate
        ihighbeta = min(find(profile(start_gate:(length(profile)-search_gates_above)) > threshold_beta)) ...
            + start_gate - 1;
        % Peak beta might be a little higher: find the highest beta in
        % this vicinity - but would this be the best approach?
        
        if ~isempty(ihighbeta)
            if min(profile(ihighbeta+[1:jump_gates]) ./ profile(ihighbeta)) < 1./min_jump
                % Candidate liquid water
                % Find the base
                search_start = max(1, ihighbeta-search_gates_below);
                diff_profile = diff(profile(search_start:ihighbeta));
                max_diff = max(diff_profile);
                base_liquid = min(find(diff_profile > max_diff*diff_factor))+search_start;
                % Find the top
                % First see if profile goes to zero
                top_liquid = min(find(~profile(ihighbeta+[1:search_gates_above])))+ihighbeta-1;
                if isempty(top_liquid)
                    diff_profile = diff(profile(ihighbeta+[0:search_gates_above]));
                    max_diff = max(-diff_profile);
                    top_liquid = max(find(-diff_profile > max_diff*diff_factor))+ihighbeta-1;
                end
                % Set to wet
                droplet_bit(base_liquid:top_liquid, ii) = 1;
            end
        end
        start_gate = ihighbeta+1;
    end
end

droplet_bit(find(temperature < min_temperature)) = 0;