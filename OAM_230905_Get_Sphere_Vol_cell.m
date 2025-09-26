function [Spherical_vol_cell] = OAM_230905_Get_Sphere_Vol_cell(ccell)
%% uses a binary mask of a cell or cellular structures as a matrix to calculate the volumen of a sphere with the same equivalent diameter
% ccellA=(bwlabel(ccell));
% ed=regionprops(ccellA,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
% Spherical_vol_cell=0.523*(ed.EquivDiameter)^3;

mask_cyt1=(bwlabel(ccell));% figure;imagesc(mask_cyt1)

Spherical_vol_cell=zeros(1,1);
for it=1:max(mask_cyt1(:))
    Ic=(mask_cyt1==it);
ed=regionprops(Ic,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
Spherical_vol_cell=Spherical_vol_cell+(0.523*(ed.EquivDiameter)^3);
end

