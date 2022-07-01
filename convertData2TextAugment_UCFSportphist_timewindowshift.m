function []=convertData2TextAugment_UCFSportphist_timewindowshift()

fpath = 'G:\PHD\Thesis\Code\Abrishami\ucfsport_phist\subtract_phist_flip_newTHR\subtract_phist_flip_newTHR.mat';
expDir = 'G:\PHD\Thesis\Code\Abrishami\ucfsport_phist\subtract_phist_flip_newTHR_timewindowshift\';

load(fpath);

[r1,c1]=size(subtract_phist_flip_newTHR);
file_id = c1-3;    
flip_id = c1-2;
id_id = c1-1;
class_id = c1; 
col_indeces = 1:file_id-1;
fixedheight = 126;

windowsize = 63 ;
step = 3;

clnum = length(unique(subtract_phist_flip_newTHR(:,class_id)));    
file_num = length(unique(subtract_phist_flip_newTHR(:,file_id)));


[~,fname,ext1] = fileparts(fpath);

%[expDir,~,~] = fileparts(mfilename('fullpath'));

expTemp = fullfile(expDir,['resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat']);


if exist(expTemp, 'file')
   clear 'subtract_phist_flip_newTHR'
else
    nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
    data = ones(file_num*nAug_per_file,windowsize*length(col_indeces));        
    cl = zeros(file_num*nAug_per_file,5);
    new_fid = 1;
    gg=1;
    for i=1:file_num
        idxf = subtract_phist_flip_newTHR(:,file_id) == i;
        temp = subtract_phist_flip_newTHR(idxf,:);
        [r_img] = myresize(fixedheight , temp(:,col_indeces));            

%         continue
            
        [newimgs,ids] = timewindowshift(r_img , step , windowsize);
        nim = unique(ids);
        for kk=1:length(nim)
            idx = ids == nim(kk);
            newimg = newimgs(idx,:);
            
            
            % % % % % % % % %  normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            mi = min(newimg(:));
            ma = max(newimg(:));
            newimg = (newimg-mi)./(ma-mi);                     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            data(gg,:) =  reshape(newimg,1,[]);
            cl(gg,:) = [new_fid temp(1,[file_id flip_id id_id class_id])];  
            gg = gg+1;
            

% % %             if(new_fid<300)
% % %                 imwrite(newimg , [expDir 'id' num2str(temp(1,id_id)) 'c' num2str(temp(1,class_id)) 'f' num2str(temp(1,file_id)) 'nf' num2str(new_fid) '.tif']);                        
% % %             end
            new_fid = new_fid + 1;           
        end
%         continue;
        i
    end
    
%     save([expDir 'resize_133_window_68_normal_' fname '.mat'],'data','cl','-v7.3' );
    save([expDir 'resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'data','cl','-v7.3' );
end
clear 'subtract_phist_flip_newTHR'

disp('mat file is saved')
return;
load(expTemp,'cl');
[r1,c1]=size(cl);
new_fid = c1-4;    
file_id = c1-3;    
flip_id = c1-2;
id_id = c1-1;
class_id = c1;     

file_num  = length(unique(cl(:,new_fid)));

% Load  database files
train_person = [5 6 7 8 9 10 11 12 13 14 21 22 23 24 25 26 27 28 29 30 31 32 39 40 41 42 43 44 45 46 47 48 49 50 51 52 55 56 57 58 63 64 65 66 67 68 69 70 75 76 77 78 79 80 81 82 83 88 89 90 91 92 93 94 95 102 103 104 105 106 107 108 109 110 111 112 113 114 115 120 121 122 123 124 125 126 127 128 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 ];
%addedtrain_person_2RH_2SB = [5 6 7 8 9 10 11 12 13 14 21 22 23 24 25 26 27 28 29 30 31 32 39 40 41 42 43 44 45 46 47 48 49 50 51 52 55 56 57 58 61 62 63 64 65 66 67 68 69 70 75 76 77 78 79 80 81 82 83 84 85 88 89 90 91 92 93 94 95 102 103 104 105 106 107 108 109 110 111 112 113 114 115 120 121 122 123 124 125 126 127 128 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 ];
%addedtrain_person_2R_2SB = [5 6 7 8 9 10 11 12 13 14 21 22 23 24 25 26 27 28 29 30 31 32 39 40 41 42 43 44 45 46 47 48 49 50 51 52 55 56 57 58 63 64 65 66 67 68 69 70 72 73 75 76 77 78 79 80 81 82 83 84 85 88 89 90 91 92 93 94 95 102 103 104 105 106 107 108 109 110 111 112 113 114 115 120 121 122 123 124 125 126 127 128 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 ];
%train_person = addedtrain_person_2R_2SB;

test_person = [1 2 3 4 15 16 17 18 19 20 33 34 35 36 37 38 53 54 59 60 61 62 71 72 73 74 84 85 86 87 96 97 98 99 100 101 116 117 118 119 129 130 131 132 133 134 135  ];
    
train_ind = zeros(file_num,1,'logical');
for k=1:length(train_person)
    tp = train_person(k);
    tem_tr_ind = cl(:,id_id)==tp;    
    train_ind = train_ind | tem_tr_ind;
end

test_ind = zeros(file_num,1,'logical');
for k=1:length(test_person)
    tp = test_person(k);
    tem_test_ind = ( cl(:,id_id)==tp & cl(:,flip_id)==1) ;    
    test_ind = test_ind | tem_test_ind;
end

trainlbl = cl(train_ind,class_id);
testlbl = cl(test_ind,class_id);
clear 'cl'
load(expTemp,'data');

fileid = fopen([expDir 'Train_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_2R_2SB.txt'],'w');
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



fileid = fopen([expDir 'Test_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');

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
