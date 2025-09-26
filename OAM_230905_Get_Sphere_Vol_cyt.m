function [Spherical_vol_cyt] = OAM_230905_Get_Sphere_Vol_cyt(mask_cyt)
%% uses a binary mask of a cell or cellular structures as a matrix to calculate the volumen of a sphere with the same equivalent diameter
% mask_cyt=mask_cyt_Orange;%figure;imagesc(mask1)

mask_cyt1=(bwlabel(mask_cyt));% figure;imagesc(mask_cyt1)


Spherical_vol_cyt=zeros(1,1);
for it=1:max(mask_cyt1(:))
    Ic=(mask_cyt1==it);
ed=regionprops(Ic,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
Spherical_vol_cyt=Spherical_vol_cyt+(0.523*(ed.EquivDiameter)^3);
end

