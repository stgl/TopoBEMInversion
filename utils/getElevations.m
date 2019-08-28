function [elevations,i] = getElevations()

[i, elevations] = textread('input/G_channel_elevations.txt','%f, %f');
