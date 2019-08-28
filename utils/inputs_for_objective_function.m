function [geo_map, d, Ginv_elev, channel_indexes, channel_elevations, n_K, compare_elements,lithology,outlet_indexes] = inputs_for_objective_function(theta)
% Modified by Felipe May 24, 2017

[G_elev,d,channel_indexes,outlet_indexes] = makeGForConcavity(theta);

% Expand G_elev to account for penalty constraints:

n_with_penalty = length(G_elev(:,1)) + 2;
G_elev_p = sparse(n_with_penalty,n_with_penalty);
G_elev_p(1:n_with_penalty-2,1:n_with_penalty-2) = G_elev;
G_elev_p(end-1:end,end-1:end) = speye(2,2);

Ginv_elev = inv(G_elev_p);

channel_elevations = getElevations();

[lithology, compare_elements] = getLithologies();
current_geotag_value = lithology(1);
all_geotag_values = [lithology(1)];
interval_start = 1;
number_of_entries = 1;

for i=1:length(lithology)
    if(lithology(i) ~= current_geotag_value)
        j = find(all_geotag_values == current_geotag_value);
        geo_map(number_of_entries,:) = [current_geotag_value, interval_start, i-1, j(1)];
        number_of_entries = number_of_entries + 1;
        j = find(all_geotag_values == lithology(i));
        if(length(j) == 0)
            all_geotag_values = [all_geotag_values, lithology(i)];
        end
        current_geotag_value = lithology(i);
        interval_start = i;
    end
end
j = find(all_geotag_values == current_geotag_value);
geo_map(number_of_entries,:) = [current_geotag_value, interval_start, length(lithology), j(1)];

n_K = length(all_geotag_values);