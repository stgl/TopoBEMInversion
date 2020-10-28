function [momentrate,elemarea]=MomentRate(c,v,slip,ShearMod)

% function centroids=PatchCentroid(c,v)

% centroids=zeros(size(v));
elemarea=zeros(length(v),1);
momentrate=zeros(length(v),1);

for i=1:length(v)
    
    i1=v(i,1);
    i2=v(i,2);
    i3=v(i,3);
    
    vec1=[c(i1,1)-c(i2,1),c(i1,2)-c(i2,2),c(i1,3)-c(i2,3)];
    vec2=[c(i1,1)-c(i3,1),c(i1,2)-c(i3,2),c(i1,3)-c(i3,3)];
    
    elemarea(i)=(norm(cross(vec1,vec2)))/2;
    
    momentrate(i)=ShearMod*elemarea(i)*slip(i);
    

        
end

end





