function [lithology_chan, compare_indexes] = getLithologies()
%Modified by Felipe May 2017

forbidden_codes = [0, 1, 7, 8, 9, 10, 11, 12, 13, 14, 38];

[i_chan, lithology_chan] = textread('input/G_channel_lithology.txt','%f, %f');

% compare_indexes = i_chan(find(~ismember(lithology_chan, forbidden_codes)))+1;

i_chan_newvect = 1:length(i_chan);
compare_indexes = i_chan_newvect(find(~ismember(lithology_chan, forbidden_codes)))';




