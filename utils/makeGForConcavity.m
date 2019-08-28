function [Go,d,channel_indexes,outlet_indexes] = makeGForConcavity(c)

% Read channel points:

[from_ind, from_val, to_ind, to_val, area] = textread('input/G_channels.txt', '%f, %f, %f, %f, %f');

% Read outlet points:

[outlet_index, elevation] = textread('input/G_outlets.txt','%f, %f');

% Allocate sparse matrix:

max_index = max([from_ind;to_ind;outlet_index])+1;
G = sparse(max_index,max_index);
d = zeros(length(outlet_index),1);

% Build sparse matrix for channel points:

channel_indexes = zeros(length(from_ind),1);

for i=1:length(from_ind)
   G(i,from_ind(i)+1) = from_val(i) .* area(i).^c;
   G(i,to_ind(i)+1) = to_val(i) .* area(i).^c;
   channel_indexes(i) = from_ind(i)+1;
end

% Build sparse matrix for outlets:

outlet_indexes = zeros(length(outlet_index), 1);

for i=1:length(outlet_index)
    G(i+length(from_ind),outlet_index(i)+1) = 1.0;
    d(i) = elevation(i);
    outlet_indexes(i) = outlet_index(i)+1;
end

% Reorder columns of G to treat non-outlets first, outlets second:

Go = sparse(max_index,max_index);
Go(:, 1:length(channel_indexes)) = G(:,channel_indexes);
Go(:, length(channel_indexes)+1:max_index) = G(:,outlet_indexes);