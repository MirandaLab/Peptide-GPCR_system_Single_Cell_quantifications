
function [p_nuc]=OAM_230906_Gaussian_nuclear_fit(IG,peak_cutoff,x_size,y_size,ccell)


    p_nuc=zeros(size(IG)); % allocate
    %-----gaussian fit--------------------
% 
% h = fspecial('average', 5); % structuring element for filtering (kernel)
% IG_GPU = gpuArray(IG);
% IG_GPU = imfilter(IG_GPU,h,'replicate');% filtering for bringing up spatial features
% Ignew=gather(IG_GPU);

h = fspecial('average', 5);
Ignew = imfilter(IG,h,'replicate');%

    nuc_tmp=Ignew.*ccell; % figure;imagesc(nuc_tmp); caxis([130 136])
    
    Amp=max(max(nuc_tmp));
    [ A, B ]=find(nuc_tmp==Amp);

   if size(A,1)~=1 % in case we have more than 1 max intensity pixels

        A = A(2);
        B = B(2);
   
   end

if A>=(size(nuc_tmp,1)-5) || B>=(size(nuc_tmp,2)-5) % exlusion of edge cells

       p_nuc=nan;

else
%     s_f=regionprops(logical(nuc_tmp==Amp),'Centroid'); % figure; imagesc(logical(nuc_tmp==Amp)) 
    xcenter=B;
    ycenter=A;
    wind=5;
    xcoord=max(1,xcenter-wind):min(y_size,xcenter+wind);
    ycoord=max(1,ycenter-wind):min(x_size,ycenter+wind);
    [X, Y] = meshgrid(xcoord,ycoord);
    
    sX=size(X);
    sY=size(Y);
    ok_size=wind.*2+1;
    if sX(1,1)==ok_size && sX(1,2)==ok_size && sY(1,1)==ok_size && sY(1,2)==ok_size
        In_part=nuc_tmp(round(diag(Y)),round(diag(X))); % imagesc(In_part)
        edge_med=median([In_part(1,:) In_part(end,:) In_part(:,1)' In_part(:,end)']);
        x0=xcenter;
        y0=ycenter;
        mu=[x0,y0];
        best_fit_score=1e9;% first threshold 
        best_fit_no=[0 0];
        for jj=1:1:25

            for jj2=1:1:25
                limit_h=round(sqrt(jj.*jj2));
                for ii=-(limit_h-1):2:(limit_h-1) % why 2
                    Sigma = [jj -ii; -ii jj2];

                   %            XOX=gpuArray([X(:) Y(:)]);
                   %  F = mvnpdf(XOX,mu,Sigma); % plot(F)
                   % F=  gather(F);

  F = mvnpdf([X(:) Y(:)],mu,Sigma); % plot(F)

                    F = reshape(F,length(X(:,1)),length(Y(:,1))); % imagesc(F)
                    Fmax=max(max(F));
                    Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med); % Z correction, tru normalization
                    Z=Z.*(In_part>0); % express the resuls in Z
                    tmp_score=sum(sum(abs((Z-In_part)))); 
%                     
%                      imagesc(Z);title([num2str(jj) '  ' num2str(jj2) 'score = ' num2str(tmp_score)]);
% %                     
% %                     
% %                     
%                     pause
                    if tmp_score<best_fit_score
                        best_fit_no=Sigma;%[ii jj]; % second threshold: should we establish a sigma threshold?
                        best_fit_score=tmp_score;  % third threshold: should we establish a sigma threshold?
                    end
                end
            end
        end
       
        
        %----------------- get the best solution------------------------
        %Sigma = [best_fit_no(1,2) -best_fit_no(1,1); -best_fit_no(1,1) best_fit_no(1,2)];
        %F = mvnpdf([X(:) Y(:)],mu,Sigma);
        F = mvnpdf([X(:) Y(:)],mu,best_fit_no);% figure;plot(F)
        F = reshape(F,length(X(:,1)),length(Y(:,1))); % figure;imagesc(F)
        Fmax=max(max(F));
        Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        Z=Z.*(In_part>0); % figure;imagesc(Z)
        %---------------------------------------------------------------    
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(double(Amp).*peak_cutoff));
        %---------------------------------------------------------------------     
       

    end

 p_nuc=double(p_nuc); % % figure;imagesc(p_nuc)
        % figure;imagesc(~p_nuc.*(nuc_tmp))

end
    end %size