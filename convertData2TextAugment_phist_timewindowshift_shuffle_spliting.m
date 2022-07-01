% this function is used for tensorflow code
function []=convertData2TextAugment_phist_timewindowshift_shuffle_spliting()

% fpath = 'D:\dataset_phist\kth_phist\subtract_rpca_phist_newTHR\subtract_rpca_phist_newTHR.mat' 
% expDir = 'D:\dataset_phist\kth_phist\subtract_rpca_phist_newTHR_timewindowshift\';

fpath = 'D:\dataset_phist\kth_phist\subtract_phist_newTHR\subtract_phist_newTHR.mat' 
expDir = 'D:\dataset_phist\kth_phist\subtract_phist_newTHR_timewindowshift\';


load(fpath);
subtract_phist_newTHR = shuffle_clip(subtract_phist_newTHR);

[r1,c1]=size(subtract_phist_newTHR);
file_id = c1-3;    
conditions_id = c1-2;
person_id = c1-1;
class_id = c1; 
col_indeces = 1:file_id-1;

fixedheight = 136;
windowsize = 68 ;

step = 5;

clnum = length(unique(subtract_phist_newTHR(:,class_id)));    
file_num = length(unique(subtract_phist_newTHR(:,file_id)));


[~,fname,ext1] = fileparts(fpath);
% fname = [fname 'dddd']
%[expDir,~,~] = fileparts(mfilename('fullpath'));


expTemp_cl = fullfile(expDir,['cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat']);

k800=200;
if exist(expTemp_cl, 'file')
   clear 'subtract_phist_newTHR'
   Batch_Number = floor((file_num)/k800);
   cntBatch = Batch_Number + 1;
   nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
else
    nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
    cl = zeros(file_num*nAug_per_file,5);
    new_fid = 1;
    gg=1;    
    cntBatch=1;
    data = zeros(k800*nAug_per_file,windowsize*length(col_indeces));        
    chunk = 1;
    Batch_Number = floor((file_num)/k800);    
    rand_index = randperm(max(file_num));
    for k=1:length(rand_index)
        i = rand_index(k)
        % if(i~=97), continue, end
        idxf = subtract_phist_newTHR(:,file_id) == i;
        temp = subtract_phist_newTHR(idxf,:);
        [r_img] = myresize(fixedheight , temp(:,col_indeces));
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
            data(chunk,:) =  reshape(newimg,1,[]);
            cl(gg,:) = [new_fid temp(1,[file_id conditions_id person_id class_id])];  
            gg = gg+1;
            chunk = chunk+1;
            if(new_fid<20)
                imwrite(newimg , [expDir 'p' num2str(temp(1,person_id)) 'c' num2str(temp(1,class_id)) 'con' num2str(temp(1,conditions_id)) 'f' num2str(temp(1,file_id)) 'nf' num2str(new_fid) '.tif']);                        
            end
            new_fid = new_fid + 1;           
            if(chunk>k800*nAug_per_file)               
%                 if(~exist([expDir 'data_resize_133_window_68_normal_' fname '_' num2str(cntBatch) '.mat'],'file'))
%                     save([expDir 'data_resize_133_window_68_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
%                 end
                if(~exist([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'file'))
                    save([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
                end

                save([expDir 'cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl','-v7.3' );                
                chunk = 1;
                cntBatch = cntBatch +1;
                if(cntBatch==Batch_Number+1)
                    clear 'data'
                    tem = file_num - (cntBatch-1)*k800;
                    data = zeros(tem*nAug_per_file,windowsize*length(col_indeces));        
                end
            end
            
        end
        if(k*nAug_per_file~=gg-1)
            gfggf=0;
        end
        k
    end
    if(cntBatch==Batch_Number+1)
        save([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
        save([expDir 'cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl','-v7.3' );
    end
end
return
clear 'subtract_phist_newTHR'
clear 'data'

load([expDir 'cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl');
[r1,c1]=size(cl);
new_fid = c1-4;    
file_id = c1-3;    
conditions_id = c1-2;
person_id = c1-1;
class_id = c1;     

file_num  = length(unique(cl(:,new_fid)));

% Load KTH database files
train_person = [11, 12, 13, 14, 15, 16, 17, 18,19, 20, 21, 23, 24, 25, 1, 4];
test_person = [22, 2, 3, 5, 6, 7, 8, 9, 10 ];

test_ind = false(file_num,1);
% test_ind = zeros(file_num,1,'logical');
for k=1:length(test_person)
    tp = test_person(k);
    tem_test_ind = cl(:,person_id)==tp;    
    test_ind = test_ind | tem_test_ind;
end
train_ind = ~test_ind;


% % % % % fileid = fopen([expDir 'Train_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');
% % % % % for i=1:cntBatch
% % % % %     train_cl_ind = false(file_num,1); %zeros(file_num,1,'logical');
% % % % %     clear 'data'
% % % % %     load([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(i) '.mat'],'data');   
% % % % %     [rd,c] = size(data);
% % % % %     disp(['Batch ' num2str(i) ' size = ' num2str(rd)])
% % % % %     ss = k800*nAug_per_file;
% % % % %     if(i<cntBatch)
% % % % %         batch_train_ind = train_ind((i-1)*ss+1 : i*ss);
% % % % %         if(sum(batch_train_ind)==0)
% % % % %             continue;
% % % % %         end
% % % % %         train_cl_ind((i-1)*ss+1 : i*ss) = batch_train_ind;
% % % % %         trainlbl = cl(train_cl_ind,class_id);
% % % % %     else
% % % % %         batch_train_ind = train_ind((i-1)*ss+1 : end);
% % % % %         if(sum(batch_train_ind)==0)
% % % % %             continue;
% % % % %         end        
% % % % %         train_cl_ind((i-1)*ss+1 : end) = batch_train_ind;
% % % % %         trainlbl = cl(train_cl_ind,class_id);        
% % % % %     end
% % % % %         
% % % % %     r = sum(train_cl_ind);
% % % % %     c = size(data,2);
% % % % % 
% % % % %     groundTruth = zeros(r , clnum );
% % % % %     for p=1:r
% % % % %        groundTruth(p,trainlbl(p))=1;
% % % % %     end
% % % % % 
% % % % %     s=[];cc=[];
% % % % %     for jk=1:c
% % % % %         s = [s ' %.4f'];
% % % % %     end
% % % % %     for jk=1:clnum
% % % % %         cc = [cc ' %u'];
% % % % %     end
% % % % %     fprintf(fileid,['|labels' cc ' |features' s '\n'],cat(2,groundTruth,data(batch_train_ind,:))');    
% % % % %     disp(['Batch ' num2str(i) ' is write'])
% % % % % end
% % % % % fclose(fileid);

disp('train end test strat')

fileid = fopen([expDir 'bbbbTest_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');
for i=1:cntBatch
    test_cl_ind = false(file_num,1);
    clear 'data'
    load([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(i) '.mat'],'data');   
    [rd,c] = size(data);
    disp(['Batch ' num2str(i) ' size = ' num2str(rd)])
    ss = k800*nAug_per_file;
    if(i<cntBatch)
        batch_test_ind = test_ind((i-1)*ss+1 : i*ss);
        if(sum(batch_test_ind)==0)
            continue;
        end
        test_cl_ind((i-1)*ss+1 : i*ss) = batch_test_ind;
        testlbl = cl(test_cl_ind,class_id);
    else
        batch_test_ind = test_ind((i-1)*ss+1 : end);
        if(sum(batch_test_ind)==0)
            continue;
        end        
        test_cl_ind((i-1)*ss+1 : end) = batch_test_ind;
        testlbl = cl(test_cl_ind,class_id);        
    end
        
    r = sum(test_cl_ind);
    c = size(data,2);

    groundTruth = zeros(r , clnum );
    for p=1:r
       groundTruth(p,testlbl(p))=1;
    end

    s=[];cc=[];
    for jk=1:c
        s = [s ' %.4f'];
    end
    for jk=1:clnum
        cc = [cc ' %u'];
    end
    fprintf(fileid,['|labels' cc ' |features' s '\n'],cat(2,groundTruth,data(batch_test_ind,:))');    
    disp(['Batch ' num2str(i) ' is write'])
end
fclose(fileid);


end

function [imgs,ids] = timewindowshift_flip(img , step , windowsize)
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
        imgs = [imgs;fliplr(img(row_idx,:))];
        ids = [ids;repmat([k+1],length(row_idx),1)];
        
        st = st+(step-1);
        en = en+(step-1);
        k = en;
    end
    
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
%        RGB = gpuArray(img);
       newimg = imresize(img, [fixedheight w]);      
       newimg = double(gather(newimg));
    end
end


function [subtract_phist_newTHR]=shuffle_clip(subtract_phist)
    [r1,c1]=size(subtract_phist);
    file_id = c1-3;    
    conditions_id = c1-2;
    person_id = c1-1;
    class_id = c1; 
    col_indeces = 1:file_id-1;
    num_fileid = unique(subtract_phist(:,file_id));
    rand_index = randperm(max(num_fileid));
    subtract_phist_newTHR = zeros(size(subtract_phist));
    st=0;
    en = 0;
    for i=1: length(rand_index)
        temp = [];
        idx = subtract_phist(:,file_id) == rand_index(i);
        temp = subtract_phist(idx,:);
        en = st + size(temp,1);
        subtract_phist_newTHR(st+1:en,:) = temp(:,:);
        st = en;        
    end
end
