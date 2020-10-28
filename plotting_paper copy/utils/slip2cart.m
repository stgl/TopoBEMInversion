function [slip_xyz]=slip2cart(c,v,slip_str,slip_dip,slip_nor)


% c= faults.c.*1e3;
% v= faults.v;
% slip_str= str_MCMC;
% slip_dip= dip_MCMC;
% slip_nor= nor_MCMC;

% [momentrate,elemarea]=MomentRate(faults.c.*1e3,faults.v,slip_total_MCMC.*1e-3,ShearMod);
% function [momentrate,elemarea]=MomentRate(c,v,slip,ShearMod)


slip_xyz=zeros(length(v),3);


for i=1:length(v)
    
    i1=v(i,1);
    i2=v(i,2);
    i3=v(i,3);
    
    vec1=[c(i1,1)-c(i2,1),c(i1,2)-c(i2,2),c(i1,3)-c(i2,3)];
    vec2=[c(i1,1)-c(i3,1),c(i1,2)-c(i3,2),c(i1,3)-c(i3,3)];
    
    n= cross(vec1,vec2);   %find normal vector to triangle
    un = n./(norm(n));   %compute unit vector
    if un(3) < 0
        un=un.*-1;
    end
    
    %     fdd(i,:) = [ddx(i,:) ddy(i,:) ddz(i,:)]; %dip strike normal displacement discontinuty components
    
    ns = cross(n,[0 0 -1]);
    uns = ns./(norm(ns));%vector parallel to strike
    
    nd = cross(n,ns);
    und = nd./(norm(nd));%vector parallel to dip


    ss = slip_str(i)*uns; % strike displacements
    ds = slip_dip(i)*und; % dip displacements
    op = slip_nor(i)*un; % opening
    slip_xyz(i,:) = ss+ds+op;%total displacement discontinuity
       
    

        
end

end