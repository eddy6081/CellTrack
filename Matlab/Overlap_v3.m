% 
min_iou=0.5; %based on dataset. Seems not great.
% this parameter makes it so the program asks user LESS input the lower
% this number is. The closer to 1, the more the two overlapping images must
% exactly match in order to not request user's input.
% also check line 226 for minimum pixel overlap to ask user input.
min_overlap_pixel=50;


%load dataset
% have user point to dataset using uigetdir
% next use dir(*tif) to pull out all the tif files. 
disp('Select directory of binary files for first stack')
d1_path = uigetdir(pwd,'Select directory of binary files for first stack');
d1_ims = dir([d1_path,filesep,'*.tif']);
d1_ims=d1_ims(~ismember({d1_ims.name},{'.','..'}));

%load subsequent dataset
disp('Select directory of binary files for second subsequent stack')
d2_path = uigetdir(pwd,'Select directory of binary files for second subsequent stack');
d2_ims = dir([d2_path,filesep,'*.tif']);
d2_ims=d2_ims(~ismember({d2_ims.name},{'.','..'}));

%load previous layer's overlap. 
%we want to do this because sometimes cells can extend through several
%layers. This would allow us to attach them all. 
%this is new in v3.
disp('If applicable, please select previous layers overlap files')
d3_path = uigetdir(pwd, 'If applicable, please select previous layers overlap files');
if ischar(d3_path)
    d3_ims = dir([d3_path,filesep,'*.tif']);
    d3_ims=d3_ims(~ismember({d3_ims.name},{'.','..'}));
end

idcs2=strfind(d2_path,'/');
%create a directory there for the redone images.
idcs1=strfind(d1_path,'/');
if ~exist([d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_nonoverlap'], 'dir')
   mkdir([d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_nonoverlap'])
   mkdir([d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_nonoverlap',filesep,'BW'])
end

% % % disp('Please select corresponding gray directory')
% % % d3_path = uigetdir(d1_path(1:idcs1(end)-1),'Please select corresponding gray directory');
% % % d3_ims = dir([d3_path,filesep,'*.tif']);
% % % d3_ims=d3_ims(~ismember({d3_ims.name},{'.','..'}));

if ~exist([d2_path(1:idcs2(end)-1),'_',d1_path(idcs1(end-1)+1:idcs1(end)-1),'_nonoverlap'], 'dir')
   mkdir([d2_path(1:idcs2(end)-1),'_',d1_path(idcs1(end-1)+1:idcs1(end)-1),'_nonoverlap'])
   mkdir([d2_path(1:idcs2(end)-1),'_',d1_path(idcs1(end-1)+1:idcs1(end)-1),'_nonoverlap',filesep,'BW'])
end

if ischar(d3_path)
    idcs3=strfind(d3_path,'/');
    if ~exist([d3_path(1:idcs3(end)-1),filesep,'edited'])
        mkdir([d3_path(1:idcs3(end)-1),filesep,'edited'])
    end
end

% % % disp('Please select corresponding gray directory')
% % % d4_path = uigetdir(d2_path(1:idcs2(end)-1),'Please select corresponding gray directory');
% % % d4_ims = dir([d4_path,filesep,'*.tif']);
% % % d4_ims=d4_ims(~ismember({d4_ims.name},{'.','..'}));
% % % %select gray path to look at larger picture

%generate overlap folder
if ~exist([d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_overlap'], 'dir')
   mkdir([d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_overlap'])
   mkdir([d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_overlap',filesep,'BW'])
end

save1_path=[d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_nonoverlap',filesep,'BW'];%[d1_path(1:idcs1(end)-1),'_overlap'];
save2_path=[d2_path(1:idcs2(end)-1),'_',d1_path(idcs1(end-1)+1:idcs1(end)-1),'_nonoverlap',filesep,'BW'];%[d2_path(1:idcs2(end)-1),'_overlap'];
if ischar(d3_path)
    save3_path=[d3_path(1:idcs3(end)-1),filesep,'edited'];
end
overlap_path = [d1_path(1:idcs1(end)-1),'_',d2_path(idcs2(end-1)+1:idcs2(end)-1),'_overlap',filesep,'BW'];

%for each dataset, we want to do bwareaopen and then number the objects.
%Also though do we want to make sure the numbering is consistent? What
%about divisions? I don't think we want to consider this yet. Maybe we do
%take this a frame by frame approach. The problem is, we likely want to
%remove it entirely and put it into its own frame if that is true. 

%maybe for each image we need to run through and do centroid + labeling.
%Ultimately, this doesn't seem practical. I think run the intersection
%thing twice would for sure catch it, so that if it misses early on it will
%be put there later. Or better, record what 

input_frame = str2double(input('What frame should we begin at? (Enter an integer)    ', 's'));
while isnan(input_frame) || fix(input_frame) ~= input_frame || input_frame<1 || input_frame>length(d1_ims)
  input_frame = str2double(input('What frame should we begin at? (Enter an integer)    ', 's'));
end

chs = {'Red','Green','Blue'};
%lets assume a ten micron window is the max (even that isn't right) that a
%centroid can move. at a pixel2dist ratio of .656 microns/pixel, this means
%
rmax=20; %pixels in distance is max distance for prob

for im_num=input_frame:length(d1_ims)
    disp(sprintf('Reading %s',d1_ims(im_num).name))
    I1=imread([d1_path,filesep,d1_ims(im_num).name]);
    %I1=rgb2gray(I1);
    I1=im2bw(I1,0.5);
    bpixels=[I1(:,1)',I1(:,end)',I1(1,:),I1(end,:)];
    if mean(bpixels)>0.5
        %means pixels were flipped.
        I1=imread([d1_path,filesep,d1_ims(im_num).name]); %2016 Matlab doesn't have folder in Dir command.
        I1=im2bw(I1,0.5);
        I1=imcomplement(I1);
    end
    %I1=imfill(I1,'holes');
    I1=bwareaopen(I1,250,4);
    I1=imclearborder(I1);

    I2=imread([d2_path,filesep,d2_ims(im_num).name]);
    I2=im2bw(I2,0.5);
    bpixels=[I2(:,1)',I2(:,end)',I2(1,:),I2(end,:)];
    if mean(bpixels)>0.5
        %means pixels were flipped.
        I2=imread([d2_path,filesep,d2_ims(im_num).name]); %2016 Matlab doesn't have folder in Dir command.
        I2=im2bw(I2,0.5);
        I2=imcomplement(I2);
    end
    %I2=imfill(I2,'holes');
    I2=bwareaopen(I2,250,4);
    I2=imclearborder(I2);

    I1_labeled=bwlabel(I1,4);
    I2_labeled=bwlabel(I2,4);

    overlapped_pixels = cell(1,2);
    overlapped_cellnums = [];
    
    overlappedcells=zeros([size(I1_labeled),3]); %can save as RGB then parse in PreviousCell2D_v2
    
    %find object centroids in previous image (singlecellold)
    %S_prev=regionprops(singlecellold>0,'centroid');
    %centroids_prev = cat(1,S_prev.Centroid);
    
    for obj1=1:max(unique(I1_labeled))
        IOUs = zeros(max(unique(I2_labeled)),1);
        inds1=find(I1_labeled==obj1);
        for obj2=1:max(unique(I2_labeled))
            inds2=find(I2_labeled==obj2);
            %compute intersection of each object: find overlap
            intersection=length(intersect(inds1,inds2));
            union=length(inds1)+length(inds2)-intersection;
            metric=intersection/union;
            if metric>min_iou
                %they are intersecting. 
                if isempty(overlapped_cellnums)
                    overlapped_pixels{1,1} = inds1;
                    overlapped_pixels{1,2} = inds2;
                    overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                else
                    [N,~]=size(overlapped_cellnums);
                    overlapped_pixels{N+1,1} = inds1;
                    overlapped_pixels{N+1,2} = inds2;
                    overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                end
                %I1_labeled_temp=I1_labeled==obj1;
                %I1_labeled_temp(I2_labeled==obj2)=true;
                %I1_labeled_temp = Find_Connected(I1_labeled_temp, overlappedcells==1, min_iou);
                %overlappedcells(I1_labeled_temp==1)=1;
                
                %I2_labeled(I2_labeled==obj2)=0;
                %remove from both
                %obj_del=1;
            
                %we are sometimes missing this.
                %WARNING: If cells may lie entirely on top of each other, this
                %is bad to do. 
            else
                if intersection>0.95*length(inds1) || intersection>0.95*length(inds2)
                    %they are intersecting. 
                    if isempty(overlapped_cellnums)
                        overlapped_pixels{1,1} = inds1;
                        overlapped_pixels{1,2} = inds2;
                        overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                    else
                        [N,~]=size(overlapped_cellnums);
                        overlapped_pixels{N+1,1} = inds1;
                        overlapped_pixels{N+1,2} = inds2;
                        overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                    end
                    %instead, try combining to single image.
                    %if length(inds1)>length(inds2)
                    %    I1_labeled(I2_labeled==obj2)=obj1;
                    %    I2_labeled(I2_labeled==obj2)=0;
                    %else
                    %    I2_labeled(I1_labeled==obj1)=obj2;
                    %    %I1_labeled(I1_labeled==obj1)=0;
                    %    obj_del=1;
                    %end
                    %%overlappedcells(inds1)=1;
                    %I1_labeled_temp=I1_labeled==obj1;
                    %I1_labeled_temp(I2_labeled==obj2)=true;
                    %I1_labeled_temp = Find_Connected(I1_labeled_temp, overlappedcells==1, min_iou);
                    %overlappedcells(I1_labeled_temp==1)=1;
                    %I1_labeled(I1_labeled==obj1)=0;
                    %I2_labeled(I2_labeled==obj2)=0;
                    %just setting these to one may be bad. We likely want to
                    %make sure there is no 
                elseif intersection>min_overlap_pixel %It is possible a protrusion formed away from the cell. 
                    figure(100)
                    Ifuse = imfuse(I1_labeled==obj1,I2_labeled==obj2,'ColorChannels',[1 2 0]);
                    imshow(Ifuse)
                    title(sprintf('%s : first z-stack in RED, second z-stack in GREEN',d1_ims(im_num).name))
                    %ask user if it is one cell or two cells.
                    input_user = str2double(input('Is this showing one cell or two cells? (Enter an integer) 1: one cell, 2: two cells    ', 's'));
                    while isnan(input_user) || fix(input_user) ~= input_user || input_user<1 || input_user>2
                      input_user = str2double(input('Is this showing one cell or two cells? (Enter an integer) 1: one cell, 2: two cells    ', 's'));
                    end
                    if input_user==1
                        %they are intersecting. 
                        if isempty(overlapped_cellnums)
                            overlapped_pixels{1,1} = inds1;
                            overlapped_pixels{1,2} = inds2;
                            overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                        else
                            [N,~]=size(overlapped_cellnums);
                            overlapped_pixels{N+1,1} = inds1;
                            overlapped_pixels{N+1,2} = inds2;
                            overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                        end
                        
                        %one cell. Put them into overlap.
                        %I1_labeled_temp=I1_labeled==obj1;
                        %I1_labeled_temp(I2_labeled==obj2)=true;
                        %I1_labeled_temp = Find_Connected(I1_labeled_temp, overlappedcells==1, min_iou);
                        %overlappedcells(I1_labeled_temp==1)=1;
                        %%overlappedcells(inds2)=1;
                        %%I2_labeled_temp = Find_Connected(I2_labeled==obj2, overlappedcells==1, 0.01);
                        %%overlappedcells(I2_labeled_temp==1)=1;
                        %%I1_labeled(I1_labeled==obj1)=0;
                        %I2_labeled(I2_labeled==obj2)=0;
                        %obj_del=1;
                    %else
                    %    %do nothing, leave them.
                        
                    end
                    
                end
            end
        end
        %if obj_del==1
        %    I1_labeled(I1_labeled==obj1)=0;
        %end
    end
    
    
    %all objects have been identified. At this point, we want to sort them
    %into the 3 channels of overlappedcells image. If they can't fit there,
    %shove them into the I1 or I2 image. Need to delete these images from
    %I1 and I2 if put into another image. as a reminder, we now have 
    %overlapped_cellnums [cellnum in I1, cellnum in I2]
    %overlapped_pixels {[pixels for cellnum in I1}, {pixels for cellnum in I2]}
    
    %I think, first we want to combine objects that have multiple overlaps
    %do a sort rows thing.
    %those that have the same second target we should combine all those
    %pixels. 
    %No, we want this to be smooth. 
    %No, I'm right here. We really want to place the large objects first. 
    
    %find objects where there is more than one target 
    
    [N,~]=size(overlapped_cellnums);
    obj_done = zeros(N,1);
    %Begin with the first object. 
    for obj = 1:N
        added=0;
        dim=1;
        if obj_done(obj)~=1
            %first, see if there are any other indeces that match.
            inds_overlap = find(overlapped_cellnums(:,1)==overlapped_cellnums(obj,1));
            %now for each second element in each of those inds, make those
            %line up.
            i = 1;
            while i<=length(inds_overlap)
                %find indeces 
                t=find(overlapped_cellnums(:,1)==overlapped_cellnums(inds_overlap(i),1));
                t(t==inds_overlap(i))=[];
                inds_overlap = [inds_overlap; t];
                t=find(overlapped_cellnums(:,2)==overlapped_cellnums(inds_overlap(i),2));
                t(t==inds_overlap(i))=[];
                inds_overlap = [inds_overlap; t];
                inds_overlap=unique(inds_overlap);
                i=i+1;
            end
            while added==0
                %check if introducing this object to the image will add a new
                %object or not. If not, add to different channel.
                if dim<4
                    overlap_temp = overlappedcells(:,:,dim);
                    overlappedcells_labeled=bwlabel(overlap_temp,4);
                    for i=1:length(inds_overlap)
                        %obj->inds_overlap(i)
                        overlap_temp(overlapped_pixels{inds_overlap(i),1})=1;
                        overlap_temp(overlapped_pixels{inds_overlap(i),2})=1;
                    end
                    %label, 
                    overlap_temp_labeled = bwlabel(overlap_temp,4);
                    if length(unique(overlappedcells_labeled))<length(unique(overlap_temp_labeled))
                        %make overlappedcells(:,:,dim)=overlap_temp.
                        overlappedcells(:,:,dim)=overlap_temp;
                        %delete from I1 and I2.
                        for i=1:length(inds_overlap)
                            I1(overlapped_pixels{inds_overlap(i),1})=0;
                            I2(overlapped_pixels{inds_overlap(i),2})=0;
                        end
                        added=1;
                    %elseif %statement for if multiple intersections of same object.
                    else
                        dim=dim+1;
                    end
                elseif dim==4
                    %try adding to I1;
                    I1_temp = I1;
                    for i=1:length(inds_overlap)
                        I1_temp(overlapped_pixels{inds_overlap(i),2})=1;
                    end
                    I1_temp_labeled = bwlabel(I1_temp,4);
                    if length(unique(I1_labeled))==length(unique(I1_temp_labeled))
                        %this is good, that means the addition of the object
                        %didn't intersect two different objects, so it can be
                        %added.
                        for i=1:length(inds_overlap)
                            I1(overlapped_pixels{inds_overlap(i),2})=1; %add object from image 2 to image 1
                            I2(overlapped_pixels{inds_overlap(i),2})=0; %delete from image 2
                        end
                        added=1;
                    else
                        dim=dim+1;
                    end
                elseif dim==5
                    %try adding to I2;
                    I2_temp = I2;
                    for i=1:length(inds_overlap)
                        I2_temp(overlapped_pixels{inds_overlap(i),1})=1;
                    end
                    %I2_temp(overlapped_pixels{inds_overlap(i),1})=1;
                    I2_temp_labeled = bwlabel(I2_temp,4);
                    if length(unique(I2_labeled))==length(unique(I2_temp_labeled))
                        %this is good, that means the addition of the object
                        %didn't intersect two different objects, so it can be
                        %added.
                        for i=1:length(inds_overlap)
                            I2(overlapped_pixels{inds_overlap(i),1})=1; %add object from image 1 to image 2
                            I1(overlapped_pixels{inds_overlap(i),1})=0; %delete from image 1
                        end
                        added=1;
                    else
                        dim=dim+1;
                    end
                else
                    %in this case, we don't have a contingency planned. May
                    %need to come up with something. 
                    %best would be to tell user to plot the two overlapping
                    %cells, and then look at each channel and I1 and I2 to show
                    %you cannot put the object anywhere. In this case, the user
                    %needs to figure something out.
                    error('The two overlapping cells cannot fit into any of the 3 channels in overlap image nor the two images. It is recommended the user plot the two overlapping cells using "overlapped_pixels", and then manually check "overlappedcells" images and "I1" and "I2" to make sure there is not a coding error.') 
                end
            end
            obj_done(inds_overlap)=1;
        end
    end
    
    I1_labeled = bwlabel(I1,4);
    I2_labeled = bwlabel(I2,4);
    
    %the next step is to find out if any cells in I2_labeled or
    %overlappedcells overlap with the images from the previous overlap
    %stack
    if ischar(d3_path)
        %than prev_overlap was given.
        %need to load the old image. 
        I3_in = imread([d3_path,filesep,d3_ims(im_num).name]);
        %check if multichannel.
        [~,~,CH]=size(I3_in);
        dims = [];
        overlapped_pixels = cell(1,2);
        overlapped_cellnums = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%CHECKING FOR OVERLAP BETWEEN THE OLD AND NEW OVERLAP IMAGES%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ch=1:CH
            %first, see if there is any overlap 
            I3 = I3_in(:,:,ch);
            I3=im2bw(I3,0.5);
            bpixels=[I3(:,1)',I3(:,end)',I3(1,:),I3(end,:)];
            if mean(bpixels)>0.5
                %means pixels were flipped.
                I3=imcomplement(I3);
            end
            %I2=imfill(I2,'holes');
            I3=bwareaopen(I3,250,4);
            I3=imclearborder(I3);
            
            %now, for each channel in overlappedcells, see if there is
            %overlap, mostly this is repeating the same code as above.
            
            %another for loop, to work through the channels of
            %overlappecells.
            I3_labeled = bwlabel(I3,4);
            for dim=1:3 %channels of overlappecells
                %look for overlap between objects in I3 and
                %overlappedcells(:,:,dim).
                I4 = overlappedcells(:,:,dim);
                I4_labeled = bwlabel(I4,4);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for obj1=1:max(unique(I3_labeled))
                    %IOUs = zeros(max(unique(I2_labeled)),1);
                    inds1=find(I3_labeled==obj1);
                    for obj2=1:max(unique(I4_labeled))
                        inds2=find(I4_labeled==obj2);
                        %compute intersection of each object: find overlap
                        intersection=length(intersect(inds1,inds2));
                        if intersection>min_overlap_pixel %It is possible a protrusion formed away from the cell. 
                            figure(100)
                            Ifuse = imfuse(I3_labeled==obj1,I4_labeled==obj2,'ColorChannels',[1 2 0]);
                            imshow(Ifuse)
                            title(sprintf('%s : Old overlap CH %s in RED + Current Overlap CH %s in GREEN ', d1_ims(im_num).name,chs{ch},chs{dim}))
                            %ask user if it is one cell or two cells.
                            input_user = str2double(input('Is this showing one cell or two cells? (Enter an integer) 1: one cell, 2: two cells    ', 's'));
                            while isnan(input_user) || fix(input_user) ~= input_user || input_user<1 || input_user>2
                              input_user = str2double(input('Is this showing one cell or two cells? (Enter an integer) 1: one cell, 2: two cells    ', 's'));
                            end
                            if input_user==1
                                %they are intersecting. 
                                if isempty(overlapped_cellnums)
                                    overlapped_pixels{1,1} = inds1;
                                    overlapped_pixels{1,2} = inds2;
                                    overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                                    dims=[dims;ch,dim];
                                else
                                    [N,~]=size(overlapped_cellnums);
                                    overlapped_pixels{N+1,1} = inds1;
                                    overlapped_pixels{N+1,2} = inds2;
                                    overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                                    dims=[dims;ch,dim];
                                end
                            end

                        end
                    end

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
            
        %record the dimension of each image the object is in, as well
        %as the object label number.

        %check if there is overlap with these cells and I2. I2 is
        %single channeled. 

        %now, for each channel in overlappedcells, see if there is
        %overlap, mostly this is repeating the same code as above.

        %another for loop, to work through the channels of
        %overlappecells.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Put the images into ovlerlap. The test is now that if we introduce
        %a new object, the number of objects should remain the same. If
        %not, we need to remove both objects and put them into a different
        %channel. Again, if that doesn't work, we need to tell the user
        %there is an error. The fix would be to put it into one of the
        %other images.
        %I think it would be easier to create a new image.
        [N,~]=size(overlapped_cellnums);
        obj_done = zeros(N,1);
        %Begin with the first object. 
        for obj = 1:N
            added=0;
            dim=1;
            if obj_done(obj)~=1
                %first, see if there are any other indeces that match.
                inds_overlap = find(overlapped_cellnums(:,1)==overlapped_cellnums(obj,1) & dims(:,1)==dims(obj,1));
                %now for each second element in each of those inds, make those
                %line up.
                i = 1;
                while i<=length(inds_overlap)
                    %find indeces 
                    t=find(overlapped_cellnums(:,1)==overlapped_cellnums(inds_overlap(i),1) & dims(:,1)==dims(inds_overlap(i),1));
                    t(t==inds_overlap(i))=[];
                    inds_overlap = [inds_overlap; t];
                    t=find(overlapped_cellnums(:,2)==overlapped_cellnums(inds_overlap(i),2) & dims(:,2)==dims(inds_overlap(i),2));
                    t(t==inds_overlap(i))=[];
                    inds_overlap = [inds_overlap; t];
                    inds_overlap=unique(inds_overlap);
                    i=i+1;
                end
                while added==0
                    %this is where we want to see the dimensions where the
                    %original object in dims(obj,2) are 
                    %I suppose it doesn't matter where we add it, but once
                    %added, we need to delete it. %simple, go to
                    %overlapcells(:,:,dims(obj,2))=0;
                    %we are either looking for number of objects + 1 or the
                    %same number... Not always, sometimes you could add it
                    %to the wrong object.
                    
                    %check if introducing this object to the image will add a new
                    %object or not. If not, add to different channel.
                    if dim<4
                        overlap_temp = overlappedcells(:,:,dim);
                        overlappedcells_labeled=bwlabel(overlap_temp,4);
                        for i=1:length(inds_overlap)
                            %obj->inds_overlap(i)
                            overlap_temp(overlapped_pixels{inds_overlap(i),1})=1;
                            overlap_temp(overlapped_pixels{inds_overlap(i),2})=1;
                        end
                        %label, 
                        overlap_temp_labeled = bwlabel(overlap_temp,4);
                        
                        %if we can introduce to a new channel, we should
                        %have N+1 objects in that channel.
                        if length(unique(overlappedcells_labeled))<length(unique(overlap_temp_labeled))
                            %make overlappedcells(:,:,dim)=overlap_temp.
                            overlappedcells(:,:,dim)=overlap_temp;
                            %delete from old overlap image and new overlap
                            %image (check though that
                            %dim~=dims(inds_overlap(i),2)
                            for i=1:length(inds_overlap)
                                I3 = I3_in(:,:,dims(inds_overlap(i),1));
                                I3(overlapped_pixels{inds_overlap(i),1})=0;
                                I3_in(:,:,dims(inds_overlap(i),1))=I3;
                                %above is not memory efficient. 
                                if dim~=dims(inds_overlap(i),2)
                                    %delete from overlappedcells.
                                    I4 = overlappedcells(:,:,dims(inds_overlap(i),2));
                                    I4(overlapped_pixels{inds_overlap(i),2})=0;
                                    overlappedcells(:,:,dims(inds_overlap(i),2)) = I4;
                                end
                            end
                            added=1;
                            
                        %if we can add it to existing channel with same
                        %object.
                        elseif sum(dims(inds_overlap,2)==dim)>0 && length(unique(overlappedcells_labeled))==length(unique(overlap_temp_labeled))
                            %two conditions: if sum(dims(inds_overlap(i),2)==dim)>0 && length(unique(overlappedcells_labeled))==length(unique(overlap_temp_labeled))
                            %add to 
                            overlappedcells(:,:,dim)=overlap_temp;
                            %delete from old overlap image and new overlap
                            %image (check though that
                            %dim~=dims(inds_overlap(i),2)
                            for i=1:length(inds_overlap)
                                I3 = I3_in(:,:,dims(inds_overlap(i),1));
                                I3(overlapped_pixels{inds_overlap(i),1})=0;
                                I3_in(:,:,dims(inds_overlap(i),1))=I3;
                                %above is not memory efficient. 
                                if dim~=dims(inds_overlap(i),2)
                                    %delete from overlappedcells.
                                    I4 = overlappedcells(:,:,dims(inds_overlap(i),2));
                                    I4(overlapped_pixels{inds_overlap(i),2})=0;
                                    overlappedcells(:,:,dims(inds_overlap(i),2)) = I4;
                                end
                            end
                            added=1;
                        else
                            dim=dim+1;
                        end
                        
                    
                    else
                        %in this case, we don't have a contingency planned.
                        %We would want user to add the cell to somewhere,
                        %the problem is we particularly want it in overlap.
                        %
                        error('The two overlapping cells cannot fit into any of the 3 channels in overlap image. It is recommended the user plot the two overlapping cells using "overlapped_pixels", and then manually check "overlappedcells" images and "I1" and "I2" to make sure there is not a coding error.') 
                    end
                end
                obj_done(inds_overlap)=1;
                %delete from overlap_image? only if it added to a different
                %dim. 
                %so the check is if dim~=dims( , 2)
            end
        end
        
        %as of right now, we are returning the following:
        %I3_in is the OLD overlap image that has been changed if there is
        %overlap.
        %overlappecells is the NEW overlap image.
        %I2_labeled = bwlabel(I2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHECKING FOR OVERLAP BETWEEN THE OLD OVERLAP IMAGES AND SECOND Z%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dims = [];
        overlapped_pixels = cell(1,2);
        overlapped_cellnums = [];
        for ch=1:CH
            %first, see if there is any overlap 
            I3 = I3_in(:,:,ch);
            I3=im2bw(I3,0.5);
            bpixels=[I3(:,1)',I3(:,end)',I3(1,:),I3(end,:)];
            if mean(bpixels)>0.5
                %means pixels were flipped.
                I3=imcomplement(I3);
            end
            %I2=imfill(I2,'holes');
            I3=bwareaopen(I3,250,4);
            I3=imclearborder(I3);
            
            %now, for each channel in overlappedcells, see if there is
            %overlap, mostly this is repeating the same code as above.
            
            %another for loop, to work through the channels of
            %overlappecells.
            I3_labeled = bwlabel(I3,4);
            for obj1=1:max(unique(I3_labeled))
                %IOUs = zeros(max(unique(I2_labeled)),1);
                inds1=find(I3_labeled==obj1);
                for obj2=1:max(unique(I2_labeled))
                    inds2=find(I2_labeled==obj2);
                    %compute intersection of each object: find overlap
                    intersection=length(intersect(inds1,inds2));
                    if intersection>min_overlap_pixel %It is possible a protrusion formed away from the cell. 
                        figure(100)
                        Ifuse = imfuse(I3_labeled==obj1,I2_labeled==obj2,'ColorChannels',[1 2 0]);
                        imshow(Ifuse)
                        title(sprintf('%s : Old overlap CH %s in RED + 2nd z stack in GREEN: ', d1_ims(im_num).name,chs{ch}))
                        %ask user if it is one cell or two cells.
                        input_user = str2double(input('Is this showing one cell or two cells? (Enter an integer) 1: one cell, 2: two cells    ', 's'));
                        while isnan(input_user) || fix(input_user) ~= input_user || input_user<1 || input_user>2
                          input_user = str2double(input('Is this showing one cell or two cells? (Enter an integer) 1: one cell, 2: two cells    ', 's'));
                        end
                        if input_user==1
                            %they are intersecting. 
                            if isempty(overlapped_cellnums)
                                overlapped_pixels{1,1} = inds1;
                                overlapped_pixels{1,2} = inds2;
                                overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                                dims=[dims;ch,1];
                            else
                                [N,~]=size(overlapped_cellnums);
                                overlapped_pixels{N+1,1} = inds1;
                                overlapped_pixels{N+1,2} = inds2;
                                overlapped_cellnums = [overlapped_cellnums; obj1,obj2];
                                dims=[dims;ch,1];
                            end
                        end

                    end
                end

            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Put the images into ovlerlap. 
        [N,~]=size(overlapped_cellnums);
        obj_done = zeros(N,1);
        %Begin with the first object. 
        for obj = 1:N
            added=0;
            dim=1;
            if obj_done(obj)~=1
                %first, see if there are any other indeces that match.
                inds_overlap = find(overlapped_cellnums(:,1)==overlapped_cellnums(obj,1) & dims(:,1)==dims(obj,1));
                %now for each second element in each of those inds, make those
                %line up.
                i = 1;
                while i<=length(inds_overlap)
                    %find indeces 
                    t=find(overlapped_cellnums(:,1)==overlapped_cellnums(inds_overlap(i),1) & dims(:,1)==dims(inds_overlap(i),1));
                    t(t==inds_overlap(i))=[];
                    inds_overlap = [inds_overlap; t];
                    t=find(overlapped_cellnums(:,2)==overlapped_cellnums(inds_overlap(i),2) & dims(:,2)==dims(inds_overlap(i),2));
                    t(t==inds_overlap(i))=[];
                    inds_overlap = [inds_overlap; t];
                    inds_overlap=unique(inds_overlap);
                    i=i+1;
                end
                while added==0
                    %this is where we want to see the dimensions where the
                    %original object in dims(obj,2) are 
                    %I suppose it doesn't matter where we add it, but once
                    %added, we need to delete it. %simple, go to
                    %overlapcells(:,:,dims(obj,2))=0;
                    %we are either looking for number of objects + 1 or the
                    %same number... Not always, sometimes you could add it
                    %to the wrong object.
                    
                    %check if introducing this object to the image will add a new
                    %object or not. If not, add to different channel.
                    if dim<4
                        overlap_temp = overlappedcells(:,:,dim);
                        overlappedcells_labeled=bwlabel(overlap_temp,4);
                        for i=1:length(inds_overlap)
                            %obj->inds_overlap(i)
                            overlap_temp(overlapped_pixels{inds_overlap(i),1})=1;
                            overlap_temp(overlapped_pixels{inds_overlap(i),2})=1;
                        end
                        %label, 
                        overlap_temp_labeled = bwlabel(overlap_temp,4);
                        
                        %if we can introduce to a new channel, we should
                        %have N+1 objects in that channel.
                        if length(unique(overlappedcells_labeled))<length(unique(overlap_temp_labeled))
                            %make overlappedcells(:,:,dim)=overlap_temp.
                            overlappedcells(:,:,dim)=overlap_temp;
                            %delete from old overlap image and new overlap
                            %image (check though that
                            %dim~=dims(inds_overlap(i),2)
                            for i=1:length(inds_overlap)
                                I3 = I3_in(:,:,dims(inds_overlap(i),1));
                                I3(overlapped_pixels{inds_overlap(i),1})=0;
                                I3_in(:,:,dims(inds_overlap(i),1))=I3;
                                %above is not memory efficient. 
                                %need to delete from I2.
                                I2(overlapped_pixels{inds_overlap(i),2})=0; %delete from image 2
                            end
                            added=1;
                        else
                            dim=dim+1;
                        end

                    else
                        %in this case, we don't have a contingency planned.
                        %We would want user to add the cell to somewhere,
                        %the problem is we particularly want it in overlap.
                        %
                        error('The two overlapping cells cannot fit into any of the 3 channels in overlap image. It is recommended the user plot the two overlapping cells using "overlapped_pixels", and then manually check "overlappedcells" images and "I1" and "I2" to make sure there is not a coding error.') 
                    end
                end
                obj_done(inds_overlap)=1;
                %delete from overlap_image? only if it added to a different
                %dim. 
                %so the check is if dim~=dims( , 2)
            end
        end

        I2_labeled = bwlabel(I2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            
    end
        
    
    %now write I1_labeled, I2_labeled, I3_in, and overlappecells to file.
    I1_new=zeros(size(I1_labeled));
    I1_new(I1_labeled>0)=1;
    imwrite(I1_new,[save1_path,filesep,d1_ims(im_num).name])
    I2_new=zeros(size(I2_labeled));
    I2_new(I2_labeled>0)=1;
    imwrite(I2_new,[save2_path,filesep,d2_ims(im_num).name])
    %write to file. Keep file name.
    imwrite(overlappedcells,[overlap_path,filesep,d1_ims(im_num).name]) 
    if ischar(d3_path)
        imwrite(I3_in,[save3_path,filesep,d3_ims(im_num).name])
    end
    %issue with compression? Taken care of.
    %overlap is not continuous...is that a problem? Seems like it is. Need
    %to figure that out. 
    
end
disp('Complete :)')
disp('You next need to run PreviousCell_2D_v2 on the created files')
disp('Good luck!')
close all

