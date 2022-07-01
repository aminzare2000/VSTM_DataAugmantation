function []=convertData2TextAugment_MSR3D_depthmap_timewindowshift()

% fpath = 'F:\VDS2\MSRAction3D\Depthmap_MSR3D\Depthmap_MSR3D.mat';
fpath = 'F:\VDS2\MSRAction3D\mapIndex_Depthmap_AS1_MSR3D\mapIndex_Depthmap_AS1_MSR3D.mat';
% fpath = 'H:\BaiduNetdiskDownload\MSRDailyAct3D\mapIndex_Depthmap_AS1_MSR3D.mat';

% expDir = 'F:\VDS2\MSRAction3D\Depthmap_MSR3D_timewindowshift\';
expDir = 'F:\VDS2\MSRAction3D\mapIndex_Depthmap_AS1_MSR3D_timewindowshift\';
%expDir = 'H:\BaiduNetdiskDownload\MSRDailyAct3D\';

load(fpath);

[~,c1]=size(mapIndex_Depthmap_AS1_MSR3D);
id_id = c1-5;    
file_id = c1-4;    
flip_id = c1-3;
subj_id = c1-2;
cond_id = c1-1;
class_id = c1; 
col_indeces = 1:id_id-1;
fname_AS = '';

fixedheight = 50;

windowsize = 40 ;
step = 3;

clnum = length(unique(mapIndex_Depthmap_AS1_MSR3D(:,class_id)));    
file_num = length(unique(mapIndex_Depthmap_AS1_MSR3D(:,id_id)));


[~,fname,ext1] = fileparts(fpath);

%[expDir,~,~] = fileparts(mfilename('fullpath'));

expTemp = fullfile(expDir,['resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname_AS fname '.mat']);


if exist(expTemp, 'file')
   clear 'mapIndex_Depthmap_AS1_MSR3D'
else
    nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
    data = ones(file_num*nAug_per_file,windowsize*length(col_indeces));        
    cl = zeros(file_num*nAug_per_file,7);
    new_fid = 1;
    gg=1;
    for i=1:file_num
        idxf = mapIndex_Depthmap_AS1_MSR3D(:,id_id) == i;
% % %         if( sum(idxf)==0 ) 
% % %             continue;
% % %         end
        temp = mapIndex_Depthmap_AS1_MSR3D(idxf,:);
        [r_img] = myresize(fixedheight , temp(:,col_indeces));            

%         continue
            
        [newimgs,ids] = timewindowshift(r_img , step , windowsize);
        nim = unique(ids);
        for kk=1:length(nim)
            idx = ids == nim(kk);
            newimg = newimgs(idx,:);
            
            
            % % % % % % % % %  normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % % % % % % % % % % % % %             mi = min(newimg(:));
% % % % % % % % % % % % % % %             ma = max(newimg(:));
% % % % % % % % % % % % % % %             newimg = (newimg-mi)./(ma-mi);                     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            data(gg,:) =  reshape(newimg,1,[]);
            cl(gg,:) = [new_fid temp(1,[id_id file_id flip_id subj_id cond_id class_id])];   


            gg = gg+1;
            

            if(new_fid<300)
                imwrite(newimg , [expDir  'c' num2str(temp(1,class_id)) 'f' num2str(temp(1,id_id)) 'nf' num2str(new_fid) '.tif']);                        
            end
            new_fid = new_fid + 1;           
        end
%         continue;
        i
    end
    
%     save([expDir 'resize_133_window_68_normal_' fname '.mat'],'data','cl','-v7.3' );
    save([expDir 'resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname_AS fname '.mat'],'data','cl','-v7.3' );
end
clear 'mapIndex_Depthmap_AS1_MSR3D'

disp('mat file is saved')

load(expTemp,'cl');
[r1,c1]=size(cl);
new_fid = c1-6;    
id_id = c1-5;    
file_id = c1-4;    
flip_id = c1-3;
subj_id = c1-2;
cond_id = c1-1;
class_id = c1; 
     

file_num  = length(unique(cl(:,new_fid)));

train_person= [1,2,8,9,10];
test_person = [3,4,5,6,7];

    
% train_ind = zeros(file_num,1,'logical');
train_ind = false(file_num,1);
for k=1:length(train_person)
    tp = train_person(k);
    tem_tr_ind = cl(:,subj_id)==tp;    
    train_ind = train_ind | tem_tr_ind;
end

% test_ind = zeros(file_num,1,'logical');
test_ind =  false(file_num,1);
for k=1:length(test_person)
    tp = test_person(k);
    tem_test_ind = ( cl(:,subj_id)==tp & cl(:,flip_id)==0) ;    
    test_ind = test_ind | tem_test_ind;
end

trainlbl = cl(train_ind,class_id);
testlbl = cl(test_ind,class_id);
clear 'cl'
load(expTemp,'data');
% fileid = fopen([expDir 'Train_resize_133_window_68_normal_' fname '.txt'],'w');
fileid = fopen([expDir 'Train_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname_AS fname '.txt'],'w');
r = sum(train_ind);
c = size(data,2);


groundTruth = zeros(r , clnum );
for i=1:r
   groundTruth(i,trainlbl(i))=1;
end

s=[];cc=[];
for i=1:c
    s = [s ' %.4f'];
end
for i=1:clnum
    cc = [cc ' %u'];
end

fprintf(fileid,['|labels' cc ' |features' s '\n'],cat(2,groundTruth,data(train_ind,:))');

fclose(fileid)

disp('train write')


% fileid = fopen([expDir 'Test_resize_133_window_68_normal_' fname '.txt'],'w');
fileid = fopen([expDir 'Test_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname_AS fname '.txt'],'w');
% [r,c]=size(rowData);
r = sum(test_ind);
c = size(data,2);


groundTruth = zeros(r , clnum );
for i=1:r
   groundTruth(i,testlbl(i))=1;
end
%AAA = cat(2,groundTruth,rowData);
s=[];cc=[];
for i=1:c
    s = [s ' %.4f'];
end
for i=1:clnum
    cc = [cc ' %u'];
end

fprintf(fileid,['|labels' cc ' |features' s '\n'],cat(2,groundTruth,data(test_ind,:))');

fclose(fileid)


end

function [imgs,ids] = timewindowshift(img , step , windowsize)
    imgs = [];
    ids = [];
    [r,c] = size(img);
    k =1;
    st = 1;
    en = windowsize;
    while k<=r
        row_idx = st:en;
        imgs = [imgs;img(row_idx,:)];
        ids = [ids;repmat([k],length(row_idx),1)];
        st = st+(step-1);
        en = en+(step-1);
        k = en;
    end
    
end


function [newimg] = myresize( fixedheight , img )
    
   
    [h,w] = size(img);
    
    midfixh = floor(fixedheight/2); 
    midh = floor(h/2); 
    
    if(h>fixedheight)   % crop histogram
        offset1 = floor((h - fixedheight)/2);
        if(offset1==0),offset1=1;end
        newimg = img( offset1:offset1+fixedheight-1 , : );
    else    %interpolation
       RGB = gpuArray(img);
       newimg = imresize(RGB, [fixedheight w]);      
       newimg = double(gather(newimg));
    end
end

function [AS]=ActionSet(ASi,data,class_id)
    AS=[];
    for i=1:length(ASi)
        idx = data(:,class_id)==ASi(i);
        AS = [ AS; data(idx,:) ];
    end
end