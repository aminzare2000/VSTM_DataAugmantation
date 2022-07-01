function []=convertData2TextAugment_OlympicSportsphist_timewindowshift_spli()

fpath = 'E:\dataset_phist\olymplicsports\subtract_phist_newTHR\subtract_phist_newTHR.mat';
load(fpath)

expDir = 'E:\dataset_phist\olymplicsports\subtract_phist_newTHR_timewindowshift\';
[r1,c1]=size(subtract_phist_newTHR);

%                             phistIndces(:,end+1)=total_file;
%                             phistIndces(:,end+1)=class_id;
%                             filenames{total_file} = C{1};

file_id = c1-1;
class_id = c1; 
col_indeces = 1:file_id-1;

fixedheight = 230;
windowsize = 100 ;
step = 5;

clnum = length(unique(subtract_phist_newTHR(:,class_id)));    
file_num = length(unique(subtract_phist_newTHR(:,file_id)));


[~,fname,ext1] = fileparts(fpath);

expTemp_cl = fullfile(expDir,['cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat']);

k800=100;

clear 'cid'
clear 'i'
if exist(expTemp_cl, 'file')
   clear 'newFixedTHR_subtract_phist'
   Batch_Number = floor((file_num)/k800);
   cntBatch = Batch_Number + 1;
   nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
else
    nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
    cl = zeros(file_num*nAug_per_file,3);
    new_fid = 1;
    gg=1;    
    cntBatch=1;
    data = zeros(k800*nAug_per_file,windowsize*length(col_indeces));        
    chunk = 1;
    Batch_Number = floor((file_num)/k800);       
    for i=1:file_num
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
            cl(gg,:) = [new_fid temp(1,[file_id class_id])];  
            gg = gg+1;
            chunk = chunk+1;
            if(new_fid<10)
                imwrite(newimg , [expDir 'f' num2str(temp(1,file_id)) 'nf' num2str(new_fid) '.tif']);                        
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
        if(i*nAug_per_file~=gg-1)
            gfggf=0;
        end
        i
    end
    if(cntBatch==Batch_Number+1)
        save([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
        save([expDir 'cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl','-v7.3' );
    end
end
clear 'subtract_phist_newTHR'
clear 'data'

load([expDir 'cl_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl');
[r1,c1]=size(cl);
new_fid = c1-2;    
file_id = c1-1;    
class_id = c1;     

file_num  = length(unique(cl(:,new_fid)));

[train_fnames]= getsplit_OlympicSports_split('F:\VDS2\OlympicSports\train');
train_file_id = [];
for i=1:length(train_fnames)
    jav = strcmp(filenames ,train_fnames{i} );
    if(sum(jav)~=1)
        error('not found train  in data')
    end
    realfile_idx = find(jav==1);            
    if(length(realfile_idx)~=1)
        error('not found train  in data')
    end
    train_file_id = [realfile_idx;train_file_id];
end


[test_fnames]= getsplit_OlympicSports_split('F:\VDS2\OlympicSports\test');
test_file_id = [];
for i=1:length(test_fnames)
    jav = strcmp(filenames ,test_fnames{i} );
    if(sum(jav)~=1)
        error('not found  test in data')
    end
    realfile_idx = find(jav==1);            
    if(length(realfile_idx)~=1)
        error('not found  test in data')
    end
    test_file_id = [realfile_idx;test_file_id];
end

train_ind = zeros(file_num,1,'logical');
for k=1:length(train_file_id)
    tp = train_file_id(k);
    tem_tr_ind = cl(:,file_id)==tp;    
    train_ind = train_ind | tem_tr_ind;
end

test_ind = zeros(file_num,1,'logical');
for k=1:length(test_file_id)
    tp = test_file_id(k);
    tem_test_ind = (cl(:,file_id)==tp) ;    
    test_ind = test_ind | tem_test_ind;
end


fileid = fopen([expDir 'Train_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');
for i=1:cntBatch
    train_cl_ind = false(file_num,1); %zeros(file_num,1,'logical');
    clear 'data'
    load([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(i) '.mat'],'data');   
    [rd,c] = size(data);
    disp(['Batch ' num2str(i) ' size = ' num2str(rd)])
    ss = k800*nAug_per_file;
    if(i<cntBatch)
        batch_train_ind = train_ind((i-1)*ss+1 : i*ss);
        if(sum(batch_train_ind)==0)
            continue;
        end
        train_cl_ind((i-1)*ss+1 : i*ss) = batch_train_ind;
        trainlbl = cl(train_cl_ind,class_id);
    else
        batch_train_ind = train_ind((i-1)*ss+1 : end);
        if(sum(batch_train_ind)==0)
            continue;
        end        
        train_cl_ind((i-1)*ss+1 : end) = batch_train_ind;
        trainlbl = cl(train_cl_ind,class_id);        
    end
        
    r = sum(train_cl_ind);
    c = size(data,2);

    groundTruth = zeros(r , clnum );
    for p=1:r
       groundTruth(p,trainlbl(p))=1;
    end

    s=[];cc=[];
    for jk=1:c
        s = [s ' %.4f'];
    end
    for jk=1:clnum
        cc = [cc ' %u'];
    end
    fprintf(fileid,['|labels' cc ' |features' s '\n'],cat(2,groundTruth,data(batch_train_ind,:))');    
    disp(['Batch ' num2str(i) ' is write'])
end
fclose(fileid);

disp('train end test strat')

fileid = fopen([expDir 'Test_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');
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

function [train_fnames]= getsplit_OlympicSports_split(splitpath)

    NUM_OF_CLASS = 16;
    listing = dir(splitpath);
    k=1;
    for i=1:length(listing)
        if(strcmp(listing(i).name,'.')~=1 && strcmp(listing(i).name,'..')~=1 && listing(i).isdir==0)              
            fl = 0;
            %                          spl = split(listing(d).name,'.');
            %                          if(size(spl,1)==2 && spl(2)=='avi')
            spl = strsplit(listing(i).name,'.');
            if(length(spl)==2 && strcmp(spl(2),'txt')==1)
                fl=1;
                tempppath = [splitpath '\' listing(i).name];                
                disp(tempppath)
                fid = fopen(tempppath);
                 while 1
                   tline = fgetl(fid);
                   if tline==-1
                     break
                   end   
                   train_fnames{k} = tline;
                   k=k+1;
                 end                
                fclose(fid);
            end
            
        end
    end
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