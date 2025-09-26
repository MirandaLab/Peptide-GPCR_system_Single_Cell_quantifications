function [Spherical_vol_nuc] = OAM_230905_Get_Sphere_Vol_nuc(mask_nuc)
%% uses a binary mask of a cell or cellular structures as a matrix to calculate the volumen of a sphere with the same equivalent diameter
mask_nuc=(bwlabel(mask_nuc));
if max(mask_nuc)==0

    Spherical_vol_nuc=nan;

elseif max(mask_nuc(:))>1


            Spherical_vol_nuc=zeros(1,1);
            for it=1:max(mask_nuc(:))
                Ic=(mask_nuc==it);
            ed=regionprops(Ic,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
            Spherical_vol_nuc=Spherical_vol_nuc+(0.523*(ed.EquivDiameter)^3);
            end

else

ed=regionprops(mask_nuc,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
Spherical_vol_nuc=0.523*(ed.EquivDiameter)^3;
end
