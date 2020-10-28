theta=0.4;
[geo_map, e_outlets, Ginv_elev, channel_indexes, e_chan, n_K, ...
    ind_chan_misfit, litho_chan, outlet_indexes] = ...
    inputs_for_objective_function(theta);

[x_channels, y_channels]=getCoordinates;
x_outlets=x_channels(outlet_indexes); y_outlets=y_channels(outlet_indexes);

%%
x_channels=x_channels(channel_indexes); 
y_channels=y_channels(channel_indexes);


%% Save outlets x,y,z
fid = fopen('outlets_xyz_UTM10N.txt', 'wt');
fprintf(fid, '%4.2f   %4.2f   %4.2f\n',[x_outlets,y_outlets,e_outlets]' ); 
fclose(fid);

%% Save channels x,y,z
fid = fopen('channels_xyz_UTM10N.txt', 'wt');
fprintf(fid, '%4.2f   %4.2f   %4.2f\n',[x_channels,y_channels,e_chan]' ); 
fclose(fid);

%%
a=[1;2;3]; b=[4;5;6];
fid = fopen('test.txt', 'wt');
fprintf(fid, '%4.4f   %4.4f\n',[a,b]' ); 
fclose(fid);