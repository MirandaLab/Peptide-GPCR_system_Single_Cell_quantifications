
%% this scripts contains all steps required for the nalysis of microfluidic data 


clear all

% Save files as individual tiff uint16, and change into RGB format with plain
% numbers as input

path_h='Path\';
filenames=dir(fullfile(path_h, '*.avi'));
f_names = { filenames.name };

for it0=1:size(f_names,2)
    %% GRAY SCALE
    name1=replace(f_names{it0},'.avi','');
save_path=[path_h '\Second\' name1 '\Pos0\input\'];
mkdir(save_path)
v = VideoReader([path_h f_names{it0}]);
ending='.png';
numFrames = v.NumFrames;
for it = 1:numFrames

    frame = readFrame(v); % figure;imagesc(frame)

    if it0~=4 
    frame=rgb2gray(frame);
   % frame=frame(:,:,2);
    frame=imflatfield(frame,20);
    % frame=imresize(frame,0.5,"nearest"); % figure;imagesc(I2)

    I0=OAM_240716_simple_0255_scaling(frame,255);
frame =cat(3,I0,I0,I0);
   
    else

 frame=imresize(frame,0.5,'nearest');

    end

    %filename = [save_path sprintf(['Ashbya_' '%04d_Green.png'], i)];
    filename = [save_path num2str(it-1) ending];
    imwrite(uint8(frame), filename,"png"); % figure;imagesc(uint8(frame))
end
end


% Run the Python script for RIFE
% Use yeast vision to interpolate the time series and track teh cells: 
%% https://pypi.org/project/yeastvision/
%%



%%%% parallel extraction of fluorescent features
parpool('local', 16) % trigger the workers
clear all 
tic
% first level
path_h0=['Path\']; % path the experiments folders
% remove all files (isdir property is 0)
fold_names0=dir(path_h0);
fold_names0 = fold_names0([fold_names0(:).isdir]);
fold_names0 = fold_names0(~ismember({fold_names0(:).name},{'.','..','Interpol', 'Stitch', 'Tracks'}));
fold_names0 = { fold_names0.name};

for exp = 1%numel(fold_names0):-1:4
exp
exp_foldrs = dir([path_h0 fold_names0{exp}]);
exp_foldrs = exp_foldrs([exp_foldrs(:).isdir]);% remove '.' and '..' 
exp_foldrs = exp_foldrs(~ismember({exp_foldrs(:).name},{'.','..','Tracks','Morpho_Extracts','FL_extracts','Interpol','*.mat'}));
exp_fold_name = { exp_foldrs.name };

% Threshold for calculating nuclear features based on intensity
I_mean_modifier=1.5;
% gaussian fit parameters
peak_cutoff=0.75;
% window around the cell
cell_margin = 10; 


for it0 = 1:size(exp_fold_name,2) % from 57:75 OAM211124_Rim11_Htb1_2 had no images in the folder % imagesc(all_ob)
% it0=1

%% load tracks
%track_names = dir([path_h0 fold_names0{exp} '\' exp_fold_name{it0} '\Tracks\']);
track_names = dir([path_h0 fold_names0{exp}]);

track_names = track_names(~ismember({track_names(:).name},{'.','..','Tracks','Segs','FL_extracts','Interpol','*.mat'}));
track_names = { track_names.name };

% to get the right track ending name in file
for ih=1:numel(track_names)
    if contains(track_names{ih},'ART_Track_')
        num1=extractAfter(track_names{ih},'ART_Track_');
        break
    end
end


%load([path_h0 fold_names0{exp} '\' exp_fold_name{it0} '\Tracks\' exp_fold_name{it0} '_ART_Track_nointer_' num1 ])

load([path_h0 fold_names0{exp} '\'  exp_fold_name{it0} '_ART_Track_' num1 ])


Ipath=[path_h0 fold_names0{exp} '\' exp_fold_name{it0} '\']; 
  
    % Mask2=Mask6;  
    Mask2=Mask_down;  % volshow(Mask_down)

    x_size = size(Mask2,1);
    y_size = size(Mask2,2);

file_n=dir(fullfile(Ipath, '*.tif'));
file_n2={file_n.name};

% no_obj=unique(cell2mat(Mask2));

%Name=cell(1,1);

% get the all the channel names
Name=cell(1,1);
A=1;
    for it01=1:20 % maximum including all masks
        if ~contains(file_n2{it01}, 'Ph')
        %Name{it0,1}=char(file_n2{it0}(1,14:end)); % 14 is the position in the file_name after img_000000000
       channelName=char(file_n2{it01}(1,14:end));
       Name{A,1}=channelName;
       A=A+1;
        end
    end

if isempty(Name{1})
break
end

channels=unique(Name);

%%%% allocate the variable names
ALLDATA{1,1}='Channel_name';% str
ALLDATA{2,1}='Cell_Size'; % vector
ALLDATA{3,1}='cell_vol';
ALLDATA{4,1}='max_nuc_int1';
ALLDATA{5,1}='mean_cell_Fl1';
ALLDATA{6,1}='Conc_T_cell_Fl1';
ALLDATA{7,1}='mem_area1';
ALLDATA{8,1}='nuc_area1';
ALLDATA{9,1}='cyt_area1';
ALLDATA{10,1}='mean_fl_mem1';
ALLDATA{11,1}='std_fl_mem1';
ALLDATA{12,1}='tot_Fl_mem1';
ALLDATA{13,1}='tot_Fl_cyt1';
ALLDATA{14,1}='tot_Fl_nuc1';
ALLDATA{15,1}='mean_int_per_area_C1';
ALLDATA{16,1}='mean_int_per_area_N1';
ALLDATA{17,1}='nuc_Vol1';
ALLDATA{18,1}='cyt_Vol1';
ALLDATA{19,1}='cyt_Vol_sub1';
ALLDATA{20,1}='FL_Conc_T1';
ALLDATA{21,1}='FL_Conc_C1';
ALLDATA{22,1}='FL_Conc_N1';
ALLDATA{23,1}='FL_mean_int_N_thr1';
ALLDATA{24,1}='FL_mean_int_C_thr1';

for it1=1:size(channels,1) % channels

file1=dir(fullfile(Ipath, ['*' channels{it1}]));
file2={file1.name};
% allocation, matrixes represent cells as rows and time points as columns 
cell_Vol1=zeros(no_obj,size(Mask2,3)); 
max_nuc_int1=zeros(no_obj,size(Mask2,3)); 
mean_cell_Fl1=zeros(no_obj,size(Mask2,3)); 
Conc_T_cell_Fl1=zeros(no_obj,size(Mask2,3)); 
mem_area1=zeros(no_obj,size(Mask2,3)); 
nuc_area1=zeros(no_obj,size(Mask2,3)); 
cyt_area1=zeros(no_obj,size(Mask2,3)); 
mean_fl_mem1=zeros(no_obj,size(Mask2,3)); 
std_fl_mem1=zeros(no_obj,size(Mask2,3)); 
tot_Fl_mem1=zeros(no_obj,size(Mask2,3)); 
tot_Fl_cyt1=zeros(no_obj,size(Mask2,3)); 
tot_Fl_nuc1=zeros(no_obj,size(Mask2,3)); 
mean_int_per_area_C1=zeros(no_obj,size(Mask2,3)); 
mean_int_per_area_N1=zeros(no_obj,size(Mask2,3)); 
nuc_Vol1=zeros(no_obj,size(Mask2,3)); 
cyt_Vol1=zeros(no_obj,size(Mask2,3)); 
cyt_Vol_sub1=zeros(no_obj,size(Mask2,3)); 
FL_Conc_T1=zeros(no_obj,size(Mask2,3)); 
FL_Conc_C1=zeros(no_obj,size(Mask2,3)); 
FL_Conc_N1=zeros(no_obj,size(Mask2,3)); 
FL_mean_int_N_thr1=zeros(no_obj,size(Mask2,3)); 
FL_mean_int_C_thr1=zeros(no_obj,size(Mask2,3)); 
all_back=zeros(1,size(Mask2,3));
Cell_Size1=zeros(no_obj,size(Mask2,3)); 

%parpool('local', 2)
for c_time=1:size(Mask2,3) % images
disp(c_time)
% cell_Vol=zeros(no_obj,1); 
% max_nuc_int=zeros(no_obj,1); 

%  load the masks as a grayscaled indexed image where each cell as a cingle
%  value
    %Lcells=all_obj.cells(:,:,c_time); % figure;imagesc(Lcells) % or Mask2{1,it2}
    Lcells = Mask2(:,:,c_time); % figure;imagesc(Lcells) find(Mask2{1,1}==20583)% 23387)

    if ~isempty(Lcells) && ~sum(Lcells(:))==0

     I = imread([Ipath '/' file2{c_time}]);
     I = double(I); % figure;imagesc(I)
     I = medfilt2(I,'symmetric'); % filtering 
     I = imresize(I,0.5,"nearest"); %%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     I = imflatfield(I, 25); %%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% if sum(ttype{prefix}(1,1:4)=='_TET') == 4 || sum(ttype{prefix}(1,1:4)=='_MAT') == 4
% 
% if size(Lcells,1)~=size(I,1) || size(Lcells,2)~=size(I,2) 
% Lcells=imresize(Lcells, size(I), "nearest");
% end 
% 
% end
 
% background correction
     Lcells = imresize(Lcells,[size(I,1)  size(I,2)],"nearest"); %%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
     bck = (I.*(~Lcells)); % figure;imagesc((I))
     backgr =   median(bck(bck~=0));% 
     % I = (I-0.5*backgr); % background correction
     I = (I-1*backgr); % background correction
     I = I +1000;
     all_back(1,c_time)=backgr;
     

%% Allocation for parfor loop sliced variables
          cell_Vol=zeros(no_obj,1);
          max_nuc_int=zeros(no_obj,1);
          mean_cell_Fl=zeros(no_obj,1);
          Conc_T_cell_Fl =zeros(no_obj,1);
          mem_area =zeros(no_obj,1);
          nuc_area =zeros(no_obj,1);
          cyt_area =zeros(no_obj,1);
          mean_fl_mem =zeros(no_obj,1);
          std_fl_mem =zeros(no_obj,1);
          tot_Fl_mem =zeros(no_obj,1);
          tot_Fl_cyt =zeros(no_obj,1);
          tot_Fl_nuc =zeros(no_obj,1);
          mean_int_per_area_C =zeros(no_obj,1);
          mean_int_per_area_N =zeros(no_obj,1);
          nuc_Vol =zeros(no_obj,1);
          cyt_Vol =zeros(no_obj,1);
          cyt_Vol_sub =zeros(no_obj,1);
          FL_Conc_T =zeros(no_obj,1);
          FL_Conc_C =zeros(no_obj,1);
          FL_Conc_N =zeros(no_obj,1);
          FL_mean_int_N_thr =zeros(no_obj,1);
          FL_mean_int_C_thr =zeros(no_obj,1);
          cell_size =zeros(no_obj,1);


% the par variable is teh single cell, the par loop completes the colum of
% extraction for the time point and passes the column vector to the
% allocated matrix for the experiment

  parfor cell_no =1:no_obj % cell_no = 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  % for cell_no = 1:no_obj
        
%        if ctime==1 || ctime>cell_exists(cell_no,2)
                ccell = double(Lcells == cell_no); %  figure;imagesc(ccell); figure;imagesc(Lcells);

               % figure;imagesc(I-200*bwmorph(Lcells, "remove"));

      if sum(ccell(:))~=0
                [x_cn, y_cn] = get_wind_coord1(ccell, cell_margin);% defines a window around the cell, with cellmargin
                ccell=ccell(y_cn, x_cn);%  figure;imagesc(ccell);
                cell_size(cell_no,1)=sum(ccell(:));
                cell_Vol(cell_no,1) = OAM_230905_Get_Sphere_Vol_cell(ccell);

                    I_cell = I(y_cn,x_cn); % figure;imagesc(I_cell)
                    put_I = ccell.*I_cell; % figure;imagesc(put_I)
                    max_nuc_int(cell_no,1) = max(put_I(:)); % figure;imagesc(put_I)

                    mean_cell_Fl(cell_no,1) = sum((put_I(:)))./sum(ccell(:)); % total intensity/number of cell pixels
                    Conc_T_cell_Fl(cell_no,1) = sum((put_I(:)))./cell_Vol(cell_no,1); % total intensity / approx. volume
                    % Obtain nucleus by Gaussian fit  imagesc(mask_nuc)
                    % this gaussian corrects for multiple equally intense brightest pixels 
                    %[mask_nuc] = OAM_231006_Gaussian_nuclear_fit(I_cell,peak_cutoff); 
                    [mask_nuc] = OAM_230906_Gaussian_nuclear_fit(I_cell,peak_cutoff,x_size,y_size,ccell);  
                 
                    % visualize:  figure;imagesc(mask_nuc)
                    mask_mem = bwmorph(ccell, 'remove'); % figure;imagesc(mask_mem)
                    mem_area(cell_no,1) = sum(mask_mem(:));  % area of the membrane

if sum(mask_nuc(:))~=0
%                   mask_cyt = double(bwmorph(ccell.*bwmorph(~mask_nuc,'thicken',1),'erode',1));  
                    mask_cyt = double(ccell-mask_nuc);%.*bwmorph(~mask_nuc,'thicken',1),'erode',1));
                    % figure;imagesc(mask_cyt)
else
 mask_cyt =nan;
end
                    % visualize:  figure;imagesc(p_nuc)
                    % figure;imagesc(mask_nuc.*I_cell); clim([900 1300])
                    nuc_area(cell_no,1) = sum(mask_nuc(:));
                    cyt_area(cell_no,1) = sum(mask_cyt(:));
                    mem_fl = mask_mem.*I_cell; % figure;imagesc(mem_fl);
                    mean_fl_mem(cell_no,1) = median(mem_fl(mem_fl~=0));
                    std_fl_mem(cell_no,1) = std(mem_fl(mem_fl~=0)); % figure;imagesc(mem+ccell)
                    tot_Fl_mem(cell_no,1) = sum(mem_fl(:));
                    tot_Fl_cyt(cell_no,1) = sum(sum(mask_cyt.*I_cell));
                    tot_Fl_nuc(cell_no,1) = sum(sum(mask_nuc.*I_cell));
                    mean_int_per_area_C(cell_no,1) = sum(sum(mask_cyt.*I_cell))./sum(mask_cyt(:));
                    mean_int_per_area_N(cell_no,1) = sum(sum(mask_nuc.*I_cell))./nuc_area(cell_no,1);
                    nuc_Vol(cell_no,1) = OAM_230905_Get_Sphere_Vol_nuc(mask_nuc);
                    cyt_Vol(cell_no,1) = OAM_230905_Get_Sphere_Vol_cyt(mask_cyt);% figure;imagesc(mask_cyt)
                    cyt_Vol_sub(cell_no,1)  = cell_Vol(cell_no,1) - nuc_Vol(cell_no,1);
                    FL_Conc_T(cell_no,1)     = sum(put_I(:))./ cell_Vol(cell_no,1);
                    FL_Conc_C(cell_no,1)     = tot_Fl_cyt(cell_no,1)./cyt_Vol(cell_no,1);
                    FL_Conc_N(cell_no,1)     = tot_Fl_nuc(cell_no,1)./nuc_Vol(cell_no,1);
                   
                    
                    % Obtain nucleus by threshold constructed with "I_mean_modifier"

                    put_mod = (put_I>(I_mean_modifier.*mean(put_I(put_I>0))));
                    put_mod = bwareaopen(put_mod,5,4);
                    put_mod = imfill(put_mod,'holes');  % figure;imagesc(put_mod)
                    FL_mean_int_N_thr(cell_no,1) = sum(sum(put_mod.*put_I))./sum(put_mod(:));
                    no = put_I.*(~put_mod); % visualize: imagesc(no)
                    FL_mean_int_C_thr(cell_no,1) = sum(no(no>0))./sum(no(no>0)>0);
                 
      else


          cell_Vol(cell_no,1)=0;
          max_nuc_int(cell_no,1)=0;
          mean_cell_Fl(cell_no,1) =0;
          Conc_T_cell_Fl(cell_no,1) =0;
          mem_area(cell_no,1) =0;
          nuc_area(cell_no,1) =0;
          cyt_area(cell_no,1) =0;
          mean_fl_mem(cell_no,1) =0;
          std_fl_mem(cell_no,1) =0;
          tot_Fl_mem(cell_no,1) =0;
          tot_Fl_cyt(cell_no,1) =0;
          tot_Fl_nuc(cell_no,1) =0;
          mean_int_per_area_C(cell_no,1) =0;
          mean_int_per_area_N(cell_no,1) =0;
          nuc_Vol(cell_no,1) =0;
          cyt_Vol(cell_no,1) =0;
          cyt_Vol_sub(cell_no,1) =0;
          FL_Conc_T(cell_no,1) =0;
          FL_Conc_C(cell_no,1) =0;
          FL_Conc_N(cell_no,1) =0;
          FL_mean_int_N_thr(cell_no,1) =0;
          FL_mean_int_C_thr(cell_no,1) =0;
          cell_size(cell_no,1) =0;
       end
% 
   end %parfor

cell_Vol1(:,c_time)=cell_Vol; % imagesc(Mask2{1,1})
max_nuc_int1(:,c_time)=max_nuc_int;
mean_cell_Fl1(:,c_time)=mean_cell_Fl; 
Conc_T_cell_Fl1(:,c_time)=Conc_T_cell_Fl; 
mem_area1(:,c_time)=mem_area; 
nuc_area1(:,c_time)=nuc_area; 
cyt_area1(:,c_time)=cyt_area; 
mean_fl_mem1(:,c_time)=mean_fl_mem; 
std_fl_mem1(:,c_time)=std_fl_mem; 
tot_Fl_mem1(:,c_time)=tot_Fl_mem; 
tot_Fl_cyt1(:,c_time)=tot_Fl_cyt; 
tot_Fl_nuc1(:,c_time)=tot_Fl_nuc; 
mean_int_per_area_C1(:,c_time)=mean_int_per_area_C; 
mean_int_per_area_N1(:,c_time)=mean_int_per_area_N; 
nuc_Vol1(:,c_time)=nuc_Vol; 
cyt_Vol1(:,c_time)=cyt_Vol; 
cyt_Vol_sub1(:,c_time)=cyt_Vol_sub; 
FL_Conc_T1(:,c_time)=FL_Conc_T; 
FL_Conc_C1(:,c_time)=FL_Conc_C; 
FL_Conc_N1(:,c_time)=FL_Conc_N; 
FL_mean_int_N_thr1(:,c_time)=FL_mean_int_N_thr; 
FL_mean_int_C_thr1(:,c_time)=FL_mean_int_C_thr; 
Cell_Size1(:,c_time)=cell_size; 

    else
cell_Vol1(:,c_time)=zeros(no_obj,1);
max_nuc_int1(:,c_time)=zeros(no_obj,1);
mean_cell_Fl1(:,c_time)=zeros(no_obj,1); 
Conc_T_cell_Fl1(:,c_time)=zeros(no_obj,1); 
mem_area1(:,c_time)=zeros(no_obj,1); 
nuc_area1(:,c_time)=zeros(no_obj,1); 
cyt_area1(:,c_time)=zeros(no_obj,1); 
mean_fl_mem1(:,c_time)=zeros(no_obj,1); 
std_fl_mem1(:,c_time)=zeros(no_obj,1); 
tot_Fl_mem1(:,c_time)=zeros(no_obj,1); 
tot_Fl_cyt1(:,c_time)=zeros(no_obj,1); 
tot_Fl_nuc1(:,c_time)=zeros(no_obj,1); 
mean_int_per_area_C1(:,c_time)=zeros(no_obj,1); 
mean_int_per_area_N1(:,c_time)=zeros(no_obj,1); 
nuc_Vol1(:,c_time)=zeros(no_obj,1); 
cyt_Vol1(:,c_time)=zeros(no_obj,1); 
cyt_Vol_sub1(:,c_time)=zeros(no_obj,1); 
FL_Conc_T1(:,c_time)=zeros(no_obj,1); 
FL_Conc_C1(:,c_time)=zeros(no_obj,1); 
FL_Conc_N1(:,c_time)=zeros(no_obj,1); 
FL_mean_int_N_thr1(:,c_time)=zeros(no_obj,1); 
FL_mean_int_C_thr1(:,c_time)=zeros(no_obj,1); 
Cell_Size1(:,c_time)=zeros(no_obj,1); 
    end

end     %cell - for

%if its the first one
ALLDATA{1,it1+1}=channels{it1};
ALLDATA{2,it1+1}=Cell_Size1;
ALLDATA{3,it1+1}=cell_Vol1;
ALLDATA{4,it1+1}=max_nuc_int1; 
ALLDATA{5,it1+1}=mean_cell_Fl1;
ALLDATA{6,it1+1}=Conc_T_cell_Fl1;
ALLDATA{7,it1+1}=mem_area1;
ALLDATA{8,it1+1}=nuc_area1;
ALLDATA{9,it1+1}=cyt_area1;
ALLDATA{10,it1+1}=mean_fl_mem1;
ALLDATA{11,it1+1}=std_fl_mem1;
ALLDATA{12,it1+1}=tot_Fl_mem1;
ALLDATA{13,it1+1}=tot_Fl_cyt1;
ALLDATA{14,it1+1}=tot_Fl_nuc1;
ALLDATA{15,it1+1}=mean_int_per_area_C1;
ALLDATA{16,it1+1}=mean_int_per_area_N1;
ALLDATA{17,it1+1}=nuc_Vol1;
ALLDATA{18,it1+1}=cyt_Vol1;
ALLDATA{19,it1+1}=cyt_Vol_sub1;
ALLDATA{20,it1+1}=FL_Conc_T1;
ALLDATA{21,it1+1}=FL_Conc_C1;
ALLDATA{22,it1+1}=FL_Conc_N1;
ALLDATA{23,it1+1}=FL_mean_int_N_thr1;
ALLDATA{24,it1+1}=FL_mean_int_C_thr1;


end % channels
% add saving part

% save in a folder called "Extracs" besides the positions folders
path_save = [path_h0 fold_names0{exp}  '\FL_extracts\' ];
mkdir(path_save);

   % constrution of file name
    name3= [exp_fold_name{it0} '_FLEX' num1 ]; % '.mat'
   % save results
   save(fullfile(path_save,name3), 'ALLDATA','all_back','I_mean_modifier','peak_cutoff', 'cell_margin');

end % pos

end
toc

delete(gcp('nocreate'));



%% consolidation of time series for GFP: time series are concatenated into a single file with matrice
% where one row is a cell and each column a time point

Ipath0 ='F:\Second\Exps\FL_extracts\';

file1=dir(fullfile(Ipath0, ["*240.mat"]));
file2={file1.name};

for itA=[1:8] % itA=1

    data = load(fullfile(Ipath0, file2{itA}));
 
    % Process the loaded data as needed

if itA==1

consolidatedData = data.ALLDATA; % Initialize consolidated data with the first file's data

num1=size(consolidatedData,1)+1;

consolidatedData{num1,1}='Metadata';
consolidatedData{num1,2}=zeros(size(consolidatedData{2,2},1),13);

consolidatedData{num1,2}(:,1)=itA;

if contains(file2{itA},'uninduced') || contains(file2{itA},'Uninduced')
consolidatedData{num1,2}(:,2)=1;
elseif contains(file2{itA},'Induced') || contains(file2{itA},'induced')
consolidatedData{num1,2}(:,2)=2;
end

else

consolidatedData1 = data.ALLDATA; % size(consolidatedData1{il,2}) 

if contains(file2{itA},'uninduced') || contains(file2{itA},'Uninduced')
consolidatedData1{num1,2}(:,2)=1;
elseif contains(file2{itA},'Induced') || contains(file2{itA},'induced')
consolidatedData1{num1,2}(:,2)=2;
end



consolidatedData1{num1,2}=zeros(size(consolidatedData1{2,2},1),13);

consolidatedData1{num1,2}(:,1)=itA;

if contains(file2{itA},'uninduced') || contains(file2{itA},'Uninduced')
consolidatedData1{num1,2}(:,2)=1;
elseif contains(file2{itA},'Induced') || contains(file2{itA},'induced')
consolidatedData1{num1,2}(:,2)=2;
end



for il=2:num1

    if il~=num1
         if ~isempty(consolidatedData1{il,2})
         consolidatedData{il,2}=[consolidatedData{il,2}(:,1:30); consolidatedData1{il,2}(:,1:30)];
         end
    else

 consolidatedData{il,2}=[consolidatedData{il,2}; consolidatedData1{il,2}];
   
% end
   end

end

end
end

%% plotting GFP

ALLD=cell(1,2);

% ALLD{1,1}=ALLDATA; % induced
% ALLD{1,2}=ALLDATA; % non
thr1=8;% number of timepoints cell is present, excludes cell that enter the microfludics device
SPLOTALL=cell(1,2);
chn=2; % number fluorophores

% save the treatments
A1=cell(2,1);

for ik0=1:2

ALLD=cell(25,2);
cells_OK = consolidatedData{25,2}(:,2)==ik0; 

for ik1=1:num1 % ik=2

  

    if ik1==1
    ALLD{ik1,1} = consolidatedData{ik1,1};
    ALLD{ik1,2} = consolidatedData{ik1,2};
    else
     if ~isempty(consolidatedData{ik1,2})
    ALLD{ik1,1} = consolidatedData{ik1,1};
    ALLD{ik1,2} = consolidatedData{ik1,2}(cells_OK,:);
    end
    end
end

A1{ik0,1} = ALLD;

end


% Align time series so that all cells are analyze aligned to the time of
% birth 
for ik1=1:2

MAT1=A1{ik1,1};% cell array from SPLT_A
% find cells present since the begining 
ok_cells=find(sum(MAT1{2,2}~=0,2)>=thr1 & MAT1{2,2}(:,1)==0)';
% % cells with enough detections
% ok_cells=find(sum(MAT1{2,2}~=0,2)>=thr1)';
% % figure;imagesc(MAT1{2,2}(good_cells,:)~=0)

for ch0=1:chn
    if ch0==1

for ik=1:num1 % ik =23
SPLOT{ik,1}=MAT1{ik,1};
end
    else
        for ik=1:num1

        if ik==1 % ik=3
        SPLOT{ik,ch0}=MAT1{ik,ch0};
        % elseif ik==num1
        % SPLOT{ik,ch0}=MAT1{ik,ch0}(ok_cells,:);
        else %if ik==Gnum
        if ~isempty(MAT1{ik,ch0})    
            SPLOT{ik,ch0}=MAT1{ik,ch0}(ok_cells,:);
        end
        end
        end

    end
end

figure;imagesc(SPLOT{2,2})
plot(mean(SPLOT{2,2},"omitnan"))

%% remove artifacts by seleting only cells without nans in the channel for clustering ****** added on 08/06
% for ik2=2:num1
% SPLOT{ik2,2}(SPLOT{ik2,2}==0)=nan;
% end

for io=2:num1
A_Matrix=nan(size(SPLOT{io,2}));
MAT2=SPLOT{io,2}; 
MAT2(MAT2==0)=nan;
for io2=1:size(MAT2,1)
st_tp=find(~isnan(MAT2(io2,:)),1);
vec1=MAT2(io2,st_tp:end);
A_Matrix(io2,1:size(vec1,2))=vec1;% figure;imagesc(A_Matrix)
% plot(mean(A_Matrix,"omitnan"))
end
SPLOT{io,3}=A_Matrix;
end
SPLOT{1,3}='GFP_Aligned';

SPLOTALL{1,ik1}=SPLOT;
end


%% Plotting by biological replicates
% red is uninduced, green is induced
sae=1;
size_cut_off=500;
colrs=['r','g','k','r','b','c','r','b','m','k','b','c','r','b','m','k','b','c','r','b','m','k','b','c',...
    'g','r','k','r','b','c','r','b','m','k','b','c','r','b','m','k','b','c','r','b','m','k','b','c'];
t=1:size(SPLOT{2,2},2);
for ik4=[2:6 10:12]%  ik4==2; 
  f1=figure;
for  il=1:size(SPLOTALL,2)% loop to plot the average curve for each cluster

okcells1=SPLOTALL{1,il}{2,3}(:,1)<=500;

%okcells1=SPLOTALL{1,il}{2,3}(:,1)>=500;

Mat2=SPLOTALL{1,il}{ik4,3}(okcells1,:);
Mat2(Mat2==0)=nan;

if ik4==2

    Mat2=Mat2*0.108*0.108;

end


    sem = std(Mat2, [], 1, 'omitnan')./sqrt(size(SPLOTALL{1,il}{io,3},1));
    hold on
    s1 = shadedErrorBarV2(t, mean(Mat2,'omitnan'), 2*sem, 'lineprops', colrs(il));
    hold on
    xlim([1 size(Mat2,2)])
    %ylim(limY)
    % title(ik4, interpreter==none)
    %pause
end
titel1=SPLOT{ik4,1};
   title([titel1 '_' num2str(num1)],Interpreter="none");
   
   xlabel('Timepoints');
   ylabel([titel1 '[a.u.]'],Interpreter="none"); 


   if sae==1
saveas(f1,titel1,'pdf')
   else
   end

end
    


