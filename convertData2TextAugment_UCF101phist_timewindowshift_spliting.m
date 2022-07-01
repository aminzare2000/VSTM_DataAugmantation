function []=convertData2TextAugment_UCF101phist_timewindowshift_spliting()

% loadUCF101()
% return

fpath = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_1_20.mat';%'c:\subtract_Diag_STDphist.mat'; %'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist.mat';
%fpath = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist.mat';


expDir = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist_timewindowshift\';%'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist_timewindowshift\';
%expDir = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist_timewindowshift\';
load(fpath)
[r1,cn]=size(subtract_Diag_STDphist.ApplyEyeMakeup.subtract_phist{1});

% % % phistIndces(:,end+1) = total_file;
% % % phistIndces(:,end+1) = gid;
% % % phistIndces(:,end+1)=cell2mat(classid);                            
          

file_id = cn-2;
g_id = cn-1;
class_id = cn; 


col_indeces = 1:file_id-1;
fixedheight = 190;

windowsize = 160 ;
step = 7;


sactionUCF101 = {'ApplyEyeMakeup','ApplyLipstick','Archery','BabyCrawling','BalanceBeam','BandMarching','BaseballPitch','Basketball','BasketballDunk','BenchPress','Biking','Billiards','BlowDryHair','BlowingCandles',...
    'BodyWeightSquats','Bowling','BoxingPunchingBag','BoxingSpeedBag','BreastStroke','BrushingTeeth','CleanAndJerk','CliffDiving','CricketBowling','CricketShot','CuttingInKitchen','Diving','Drumming','Fencing',...
    'FieldHockeyPenalty','FloorGymnastics','FrisbeeCatch','FrontCrawl','GolfSwing','Haircut','Hammering','HammerThrow','HandstandPushups','HandstandWalking','HeadMassage','HighJump','HorseRace','HorseRiding','HulaHoop',...
    'IceDancing','JavelinThrow','JugglingBalls','JumpingJack','JumpRope','Kayaking','Knitting','LongJump','Lunges','MilitaryParade','Mixing','MoppingFloor','Nunchucks','ParallelBars','PizzaTossing','PlayingCello',...
    'PlayingDaf','PlayingDhol','PlayingFlute','PlayingGuitar','PlayingPiano','PlayingSitar','PlayingTabla','PlayingViolin','PoleVault','PommelHorse','PullUps','Punch','PushUps','Rafting','RockClimbingIndoor','RopeClimbing',...
    'Rowing','SalsaSpin','ShavingBeard','Shotput','SkateBoarding','Skiing','Skijet','SkyDiving','SoccerJuggling','SoccerPenalty','StillRings','SumoWrestling','Surfing','Swing','TableTennisShot','TaiChi','TennisSwing',...
    'ThrowDiscus','TrampolineJumping','Typing','UnevenBars','VolleyballSpiking','WalkingWithDog','WallPushups','WritingOnBoard','YoYo'};

clnum = length(sactionUCF101);
file_num = 13320;
% % for cid=1:length(sactionUCF101)
% %     file_num  = file_num + length(subtract_Diag_STDphist.(sactionUCF101{cid}).subtract_phist);
% % end


[~,fname,ext1] = fileparts(fpath);
fname = 'subtract_Diag_STDphist'
%[expDir,~,~] = fileparts(mfilename('fullpath'));


expTemp_cl = fullfile(expDir,['cl_filename_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat']);

k800=200;
clear 'cid'
clear 'i'
if exist(expTemp_cl, 'file')
   clear 'subtract_Diag_STDphist'
   Batch_Number = floor((file_num)/k800);
   cntBatch = Batch_Number + 1;
   nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
else
    
    nAug_per_file = floor((fixedheight-windowsize)/(step-1)) +1;
    cl = zeros(file_num*nAug_per_file,4);
    filename = cell(file_num*nAug_per_file,1);
    new_fid = 1;
    gg=1;    
    cntBatch=1;
    data = zeros(k800*nAug_per_file,windowsize*length(col_indeces));        
    chunk = 1;
    Batch_Number = floor((file_num)/k800);    
    cntshomar = 0;
    NUM_OF_CLASS_splits = [0 20 40 60 80 101];
    stf=2;
    for cid2=1:length(sactionUCF101)
        for fid2=1:length(subtract_Diag_STDphist.(sactionUCF101{cid2}).subtract_phist)
            cntshomar = cntshomar + 1;
            %idxf = subtract_phist_flip_newTHR(:,file_id) == i;
            %temp = subtract_phist_flip_newTHR(idxf,:);
            temp = subtract_Diag_STDphist.(sactionUCF101{cid2}).subtract_phist{fid2};
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
                cl(gg,:) = [new_fid temp(1,[file_id g_id class_id])]; 
                filename{gg} = subtract_Diag_STDphist.(sactionUCF101{cid2}).filename{fid2};
                if(new_fid<100)
                    A=split(filename{gg},'/');
                    imwrite(newimg , [expDir char(A(2)) '_' num2str(kk) '.tif']);                        
                end
                gg = gg+1;
                chunk = chunk+1;                
                new_fid = new_fid + 1;           
                if(chunk>k800*nAug_per_file)               
    %                 if(~exist([expDir 'data_resize_133_window_68_normal_' fname '_' num2str(cntBatch) '.mat'],'file'))
    %                     save([expDir 'data_resize_133_window_68_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
    %                 end
                    if(~exist([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'file'))
                        save([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
                    end

                    save([expDir 'cl_filename_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl','filename','-v7.3' );                
                    chunk = 1;
                    cntBatch = cntBatch +1;
                    if(cntBatch==Batch_Number+1)
                        clear 'data'
                        tem = file_num - (cntBatch-1)*k800;
                        data = zeros(tem*nAug_per_file,windowsize*length(col_indeces));        
                    end
                end

            end
            if(cntshomar*nAug_per_file~=gg-1)
                gfggf=0;
            end
            cntshomar
        end
        if(cid2~=length(sactionUCF101) && cid2>=NUM_OF_CLASS_splits(stf))
            fpath = ['G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_' num2str(NUM_OF_CLASS_splits(stf)+1) '_' num2str(NUM_OF_CLASS_splits(stf+1)) '.mat']
            load(fpath);
            stf = stf + 1;
        end
    end
    if(cntBatch==Batch_Number+1)
        save([expDir 'data_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '_' num2str(cntBatch) '.mat'],'data','-v7.3' );
        save([expDir 'cl_filename_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl','filename','-v7.3' );
    end
end
clear 'subtract_Diag_STDphist'
clear 'data'
%return
load([expDir 'cl_filename_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.mat'],'cl','filename');
[r1,cn]=size(cl);
new_fid = cn-3;
file_id = cn-2;
g_id = cn-1;
class_id = cn; 

file_num  = length(unique(cl(:,new_fid)));

% splitdir = 'E:\VDS\HMDB51\testTrainMulti_7030_splits';
% [train_fnames,test_fnames]= get_HMDB_split(1,splitdir);  

load 'E:\VDS\UCF101\UCF101TrainTestSplits-RecognitionTask\ucfTrainTestlist\ucfTrainTestlist.mat'
%load 'I:\VDS\UCF101\UCF101TrainTestSplits-RecognitionTask\ucfTrainTestlist\ucfTrainTestlist.mat'
train_fnames=trainlist01;
test_fnames=testlist01;

train_ind = zeros(file_num,1,'logical');
for i=1:length(train_fnames)
    filename1 = char(train_fnames{i,1});
    tem_tr_ind = strcmpi(filename1,filename)';                    
    if(sum(tem_tr_ind)==nAug_per_file)        
        train_ind = train_ind | tem_tr_ind';
    else
        error 'can not find'
    end    
end

test_ind = zeros(file_num,1,'logical');
for i=1:length(test_fnames)
    filename1 = char(test_fnames{i,1});
    tem_ts_ind = strcmpi(filename1,filename)';            

    if(sum(tem_ts_ind)==nAug_per_file)  
        test_ind = test_ind | tem_ts_ind';
    else
        error 'can not find'
    end    
end


fileid = fopen(['E:\ucf101_phist\subtract_Diag_STDphist_timewindowshift\' 'Train_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');
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

fileid = fopen(['E:\ucf101_phist\subtract_Diag_STDphist_timewindowshift\' 'Test_resize' num2str(fixedheight) '_window' num2str(windowsize) '_step' num2str(step-1) '_normal_' fname '.txt'],'w');
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


function [] = loadHMDB()
    saction =      {'brush_hair','cartwheel','catch','chew','clap','climb','climb_stairs',...
          'dive','draw_sword','dribble','drink','eat','fall_floor','fencing',...
          'flic_flac','golf','handstand','hit','hug','jump','kick_ball',...
          'kick','kiss','laugh','pick','pour','pullup','punch',...
          'push','pushup','ride_bike','ride_horse','run','shake_hands','shoot_ball',...
          'shoot_bow','shoot_gun','sit','situp','smile','smoke','somersault',...
          'stand','swing_baseball','sword_exercise','sword','talk','throw','turn',...
          'walk','wave'};
      
  struct_data = struct('brush_hair',[],'cartwheel',[],'catch',[],'chew',[],'clap',[],'climb',[],'climb_stairs',[],...
      'dive',[],'draw_sword',[],'dribble',[],'drink',[],'eat',[],'fall_floor',[],'fencing',[],...
      'flic_flac',[],'golf',[],'handstand',[],'hit',[],'hug',[],'jump',[],'kick_ball',[],...
      'kick',[],'kiss',[],'laugh',[],'pick',[],'pour',[],'pullup',[],'punch',[],...
      'push',[],'pushup',[],'ride_bike',[],'ride_horse',[],'run',[],'shake_hands',[],'shoot_ball',[],...
      'shoot_bow',[],'shoot_gun',[],'sit',[],'situp',[],'smile',[],'smoke',[],'somersault',[],...
      'stand',[],'swing_baseball',[],'sword_exercise',[],'sword',[],'talk',[],'throw',[],'turn',[],...
      'walk',[],'wave',[]);
      
fpath = 'G:\PHD\Thesis\Code\Abrishami\HMDB_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_1_10.mat';
load(fpath);
for cid=1:10
    struct_data.(saction{cid}).subtract_phist = subtract_Diag_STDphist.(saction{cid}).subtract_phist;
    struct_data.(saction{cid}).filename = subtract_Diag_STDphist.(saction{cid}).filename;
end

fpath = 'G:\PHD\Thesis\Code\Abrishami\HMDB_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_11_20.mat';
load(fpath);
for cid=11:20
    struct_data.(saction{cid}).subtract_phist = subtract_Diag_STDphist.(saction{cid}).subtract_phist;
    struct_data.(saction{cid}).filename = subtract_Diag_STDphist.(saction{cid}).filename;
end

fpath = 'G:\PHD\Thesis\Code\Abrishami\HMDB_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_21_30.mat';
load(fpath);
for cid=21:30
    struct_data.(saction{cid}).subtract_phist = subtract_Diag_STDphist.(saction{cid}).subtract_phist;
    struct_data.(saction{cid}).filename = subtract_Diag_STDphist.(saction{cid}).filename;
end


fpath = 'G:\PHD\Thesis\Code\Abrishami\HMDB_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_31_40.mat';
load(fpath);
for cid=31:40
    struct_data.(saction{cid}).subtract_phist = subtract_Diag_STDphist.(saction{cid}).subtract_phist;
    struct_data.(saction{cid}).filename = subtract_Diag_STDphist.(saction{cid}).filename;
end


fpath = 'G:\PHD\Thesis\Code\Abrishami\HMDB_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_41_51.mat';
load(fpath);
for cid=41:51
    struct_data.(saction{cid}).subtract_phist = subtract_Diag_STDphist.(saction{cid}).subtract_phist;
    struct_data.(saction{cid}).filename = subtract_Diag_STDphist.(saction{cid}).filename;
end

subtract_Diag_STDphist=struct_data;
save('G:\PHD\Thesis\Code\Abrishami\HMDB_phist\subtract_Diag_STDphist\subtract_Diag_STDphist.mat','subtract_Diag_STDphist','-v7.3' );
end

function [train_fnames,test_fnames]= get_HMDB_split(isplit,splitdir)
    saction =      {'brush_hair','cartwheel','catch','chew','clap','climb','climb_stairs',...
          'dive','draw_sword','dribble','drink','eat','fall_floor','fencing',...
          'flic_flac','golf','handstand','hit','hug','jump','kick_ball',...
          'kick','kiss','laugh','pick','pour','pullup','punch',...
          'push','pushup','ride_bike','ride_horse','run','shake_hands','shoot_ball',...
          'shoot_bow','shoot_gun','sit','situp','smile','smoke','somersault',...
          'stand','swing_baseball','sword_exercise','sword','talk','throw','turn',...
          'walk','wave'};
      


    for iaction = 1:length(saction)
         itr = 1;
         ite = 1;
         fname = sprintf('%s/%s_test_split%d.txt',splitdir,saction{iaction},isplit);

         fid = fopen(fname);

         while 1
           tline = fgetl(fid);
           if tline==-1
             break
           end
           [tline, u] = strtok(tline,' ');   
           u = str2num(u);

           video = sprintf('%s.avi',tline(1:end-4));

           if u==1 % ignore testing
               train_fnames{iaction}{itr} = tline;
               itr = itr + 1;
           elseif u==2
               test_fnames{iaction}{ite} = tline;
               ite = ite + 1;
           end
         end
         fclose(fid);
    end
end


function [] = loadUCF101()
        load 'E:\VDS\UCF101\UCF101TrainTestSplits-RecognitionTask\ucfTrainTestlist\ucfTrainTestlist.mat'
        %load 'I:\VDS\UCF101\UCF101TrainTestSplits-RecognitionTask\ucfTrainTestlist\ucfTrainTestlist.mat'
        sactionUCF101 = {'ApplyEyeMakeup','ApplyLipstick','Archery','BabyCrawling','BalanceBeam','BandMarching','BaseballPitch','Basketball','BasketballDunk','BenchPress','Biking','Billiards','BlowDryHair','BlowingCandles',...
            'BodyWeightSquats','Bowling','BoxingPunchingBag','BoxingSpeedBag','BreastStroke','BrushingTeeth','CleanAndJerk','CliffDiving','CricketBowling','CricketShot','CuttingInKitchen','Diving','Drumming','Fencing',...
            'FieldHockeyPenalty','FloorGymnastics','FrisbeeCatch','FrontCrawl','GolfSwing','Haircut','Hammering','HammerThrow','HandstandPushups','HandstandWalking','HeadMassage','HighJump','HorseRace','HorseRiding','HulaHoop',...
            'IceDancing','JavelinThrow','JugglingBalls','JumpingJack','JumpRope','Kayaking','Knitting','LongJump','Lunges','MilitaryParade','Mixing','MoppingFloor','Nunchucks','ParallelBars','PizzaTossing','PlayingCello',...
            'PlayingDaf','PlayingDhol','PlayingFlute','PlayingGuitar','PlayingPiano','PlayingSitar','PlayingTabla','PlayingViolin','PoleVault','PommelHorse','PullUps','Punch','PushUps','Rafting','RockClimbingIndoor','RopeClimbing',...
            'Rowing','SalsaSpin','ShavingBeard','Shotput','SkateBoarding','Skiing','Skijet','SkyDiving','SoccerJuggling','SoccerPenalty','StillRings','SumoWrestling','Surfing','Swing','TableTennisShot','TaiChi','TennisSwing',...
            'ThrowDiscus','TrampolineJumping','Typing','UnevenBars','VolleyballSpiking','WalkingWithDog','WallPushups','WritingOnBoard','YoYo'};
        
          struct_data = struct('ApplyEyeMakeup',[],'ApplyLipstick',[],'Archery',[],'BabyCrawling',[],'BalanceBeam',[],'BandMarching',[],'BaseballPitch',[],'Basketball',[],'BasketballDunk',[],'BenchPress',[],'Biking',[],'Billiards',[],'BlowDryHair',[],'BlowingCandles',[],'BodyWeightSquats',[],'Bowling',[],'BoxingPunchingBag',[],'BoxingSpeedBag',[],'BreastStroke',[],'BrushingTeeth',[],'CleanAndJerk',[],'CliffDiving',[],'CricketBowling',[],'CricketShot',[],'CuttingInKitchen',[],'Diving',[],'Drumming',[],'Fencing',[],'FieldHockeyPenalty',[],'FloorGymnastics',[],'FrisbeeCatch',[],'FrontCrawl',[],'GolfSwing',[],'Haircut',[],'Hammering',[],'HammerThrow',[],'HandstandPushups',[],'HandstandWalking',[],'HeadMassage',[],'HighJump',[],'HorseRace',[],'HorseRiding',[],'HulaHoop',[],'IceDancing',[],'JavelinThrow',[],'JugglingBalls',[],'JumpingJack',[],'JumpRope',[],'Kayaking',[],'Knitting',[],'LongJump',[],'Lunges',[],'MilitaryParade',[],'Mixing',[],'MoppingFloor',[],'Nunchucks',[],'ParallelBars',[],'PizzaTossing',[],'PlayingCello',[],'PlayingDaf',[],'PlayingDhol',[],'PlayingFlute',[],'PlayingGuitar',[],'PlayingPiano',[],'PlayingSitar',[],'PlayingTabla',[],'PlayingViolin',[],'PoleVault',[],'PommelHorse',[],'PullUps',[],'Punch',[],'PushUps',[],'Rafting',[],'RockClimbingIndoor',[],'RopeClimbing',[],'Rowing',[],'SalsaSpin',[],'ShavingBeard',[],'Shotput',[],'SkateBoarding',[],'Skiing',[],'Skijet',[],'SkyDiving',[],'SoccerJuggling',[],'SoccerPenalty',[],'StillRings',[],'SumoWrestling',[],'Surfing',[],'Swing',[],'TableTennisShot',[],'TaiChi',[],'TennisSwing',[],'ThrowDiscus',[],'TrampolineJumping',[],'Typing',[],'UnevenBars',[],'VolleyballSpiking',[],'WalkingWithDog',[],'WallPushups',[],'WritingOnBoard',[],'YoYo',[]);

      
fpath = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_1_20.mat';
%fpath = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_1_20.mat';
load(fpath);
for cid=1:20
    struct_data.(sactionUCF101{cid}).subtract_phist = subtract_Diag_STDphist.(sactionUCF101{cid}).subtract_phist;
    struct_data.(sactionUCF101{cid}).filename = subtract_Diag_STDphist.(sactionUCF101{cid}).filename;
end

      
fpath = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_21_40.mat';
%fpath = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_21_40.mat';
load(fpath);
for cid=21:40
    struct_data.(sactionUCF101{cid}).subtract_phist = subtract_Diag_STDphist.(sactionUCF101{cid}).subtract_phist;
    struct_data.(sactionUCF101{cid}).filename = subtract_Diag_STDphist.(sactionUCF101{cid}).filename;
end

fpath = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_41_60.mat';
%fpath = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_41_60.mat';
load(fpath);
for cid=41:60
    struct_data.(sactionUCF101{cid}).subtract_phist = subtract_Diag_STDphist.(sactionUCF101{cid}).subtract_phist;
    struct_data.(sactionUCF101{cid}).filename = subtract_Diag_STDphist.(sactionUCF101{cid}).filename;
end

fpath = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_61_80.mat';
%fpath = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_61_80.mat';
load(fpath);
for cid=61:80
    struct_data.(sactionUCF101{cid}).subtract_phist = subtract_Diag_STDphist.(sactionUCF101{cid}).subtract_phist;
    struct_data.(sactionUCF101{cid}).filename = subtract_Diag_STDphist.(sactionUCF101{cid}).filename;
end

fpath = 'G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_81_101.mat';
%fpath = 'I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist_81_101.mat';
load(fpath);
for cid=81:101
    struct_data.(sactionUCF101{cid}).subtract_phist = subtract_Diag_STDphist.(sactionUCF101{cid}).subtract_phist;
    struct_data.(sactionUCF101{cid}).filename = subtract_Diag_STDphist.(sactionUCF101{cid}).filename;
end
subtract_Diag_STDphist=struct_data;
%save('G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_phist_flip_newTHR\subtract_Diag_STDphist.mat','subtract_Diag_STDphist','-v7.3' );
save('G:\PHD\Thesis\Code\Abrishami\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist.mat','subtract_Diag_STDphist','-v7.3' );
%save('I:\VDS\UCF101\ucf101_phist\subtract_Diag_STDphist\subtract_Diag_STDphist.mat','subtract_Diag_STDphist','-v7.3' );
end