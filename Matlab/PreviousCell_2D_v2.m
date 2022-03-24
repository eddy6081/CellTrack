%%%
%Output:
%Input:
%Functions: Max_Fit_Circles and PreivousCell3D defined later in code.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_scale=1; %lets just fix this in the SVM script.
%Check that external functions exist within the directory.
DefaultDir=pwd;
if ~exist([DefaultDir, '/max_inscribed_circle.m'])
    error('The Matlab script "max_inscribed_circle.m" is not in the same directory')
end
if ~exist([DefaultDir, '/inpoly.m'])
    error('The Matlab script "inpoly.m" is not in the same directory')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We want to ask the user to identify the max projected stack of images
%which we will use to identify the cells. 
disp('Select max-projected stack of images  ')
gray_path=uigetdir(pwd,'select max-projected stack of images');
gray_ims=dir([gray_path,filesep,'*png']);
gray_ims=gray_ims(~ismember({gray_ims.name},{'.','..'})); %delete hidden files

if length(gray_ims)==0
    gray_ims=dir([gray_path,filesep,'*tif']);
    gray_ims=gray_ims(~ismember({gray_ims.name},{'.','..'})); %delete hidden files
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select where you want to save the printed images to.
disp('Select directory you wish to save printed images to ')
save_path=uigetdir(pwd,'Select directory you wish to save printed images to');
if strcmp(save_path,gray_path)
    error('You do not want gray path and save path to be the same')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%recursively ask user to identify z stack in consecutive order, until they
%hit cancel. 
last_check=0;
while last_check==0
    p=0;
    Fpaths={};
    while p==0
        disp("Select z-stacks in order. If complete or made an error, press cancel")
        test=uigetdir(pwd,'select z-stacks in order');
        if isstr(test)~=0
            Fpaths{end+1}=test;
        else
            p=1;
        end
    end
    Fpaths=Fpaths';
    disp(Fpaths)
    input_hf = str2double(input('Is the printed order of directories correct? (Enter an integer) 1: yes, 2: no    ', 's'));
    while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 || input_hf>2
      input_hf = str2double(input('Is the printed order of directories correct? (Enter an integer) 1: yes, 2: no    ', 's'));
    end
    if input_hf==1
        last_check=1;
    else
        disp("Having user reselect the z-stacks...")
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Code to process images and correctly identify cell numbers. 

%for each path, we want to analyze the ith frame, measuring each object and
%assigning a z along with it. These objects have cell numbers (indices)
%assigned to them depending only on where they appear on the image (l->r) by default.
%We want to correct these numbers so that we can accurately track the
%trajectory of any individual cell. 

Fnames=cell(1,length(Fpaths));
for fpath_num=1:length(Fpaths)
    d_ims=dir([Fpaths{fpath_num},filesep,'*tif']);
    d_ims=d_ims(~ismember({d_ims.name},{'.','..'})); %delete hidden files
    Fnames{fpath_num}=d_ims;
end

cellmax=0;
maxframe=length(Fnames{1}); %should double check all these filepaths have the same number of files in them


blackoutlist=[]; %%%%%%%%%%%VERY IMPORTANT
allcellnums=[];

%ask user if they want to begin at a certain frame. 
input_hf = str2double(input('Do you wish to begin at a certain frame? (Enter an integer) 1: yes, 2: no    ', 's'));
while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 || input_hf>2
  input_hf = str2double(input('Do you wish to begin at a certain frame? (Enter an integer) 1: yes, 2: no    ', 's'));
end
if input_hf==1
    %uigetfile
    disp('Please guide to the previously saved .mat file')
    [file,path]=uigetfile('*.mat');
    load(fullfile(path,file))
    %all filenames should follow workspace_checkpoint_#.mat so use idcs to
    %look for _ and .mat
    idcPeriod=strfind(file,'.');
    idcUnder=strfind(file,'_');
    
    startframe=max(AllFrames)+1;
    %str2double(file(idcUnder(end)+1:idcPeriod(1)-1));
    prev_stats=AllFrameStats(AllFrames==startframe-1);
    prev_pos=AllFramePos(AllFrames==startframe-1,:);
    prev_z=AllFrameZ(AllFrames==startframe-1);
    prev_cellnums=AllFrameNumbers(AllFrames==startframe-1);
    prev_centers=AllFrameCenters(AllFrames==startframe-1,:);
    allcellnums=AllFrameNumbers(AllFrames==startframe-1);
else
    startframe=1;
end

%ask user if they want to begin at a certain frame. 
check_bright_cells = str2double(input('Do you wish to screen cell manually - new cells will be asked to be confirmed by user? (Enter an integer) 1: yes, 2: no    ', 's'));
while isnan(check_bright_cells) || fix(check_bright_cells) ~= check_bright_cells || check_bright_cells<1 || check_bright_cells>2
  check_bright_cells = str2double(input('Invalid entry. Do you wish to screen cell manually - new cells will be asked to be confirmed by user?    ', 's'));
end

for im_num=startframe:maxframe
    for fpath_num=1:length(Fpaths)
        %load BW, make logical, run stats, positions, circles, frame number,
        %get cell numbers (add cellmax), assign z
        try
            fullpathname=[Fnames{fpath_num}(im_num).folder, filesep, Fnames{fpath_num}(im_num).name];
            fprintf('Reading %s\n',fullpathname)
            Im_in=imread(fullpathname);
        catch
            fullpathname=[Fpaths{fpath_num}, filesep, Fnames{fpath_num}(im_num).name];
            fprintf('Reading %s\n',fullpathname)
            Im_in=imread([Fpaths{fpath_num}, filesep, Fnames{fpath_num}(im_num).name]); %2016 Matlab doesn't have folder in Dir command.
        end
        %for each channel in Image
        [~,~,channels]=size(Im_in);
        for ch=1:channels
            BW=im2bw(Im_in(:,:,ch),0.5);
            bpixels=[BW(:,1)',BW(:,end)',BW(1,:),BW(end,:)];
            if mean(bpixels)>0.5
                %means pixels were flipped.
                BW=im2bw(Im_in(:,:,ch),0.5);
                BW=imcomplement(BW);
            end
            %BW=imfill(BW,'holes');
            %BW=bwareaopen(BW,300,4); %sometimes in images two objects will be connected by 8 connectivity.
            BW=bwareaopen(BW,50,4);
            %BW=imclearborder(BW);


            [h,w]=size(BW);
            %maybe remove objects that are too big?
            %BW=BW-bwareaopen(BW,h*w/80);

            if fpath_num==1
                CC = bwconncomp(BW>0,4);
                stats = regionprops(CC, {'Area','Centroid',...
                    'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
                    'PerimeterOld','FilledArea',...
                    'MajorAxisLength','MinorAxisLength','Eccentricity',...
                    'Solidity','Extent','ConvexArea','PixelIdxList'});
                pos = cat(1,stats.Centroid);
                z = repmat(fpath_num,size(pos,1),1);
                if im_num>1
                    cellnums=[1:length(z)]'+max(prev_cellnums);
                else
                    cellnums=[1:length(z)]';
                end

                %calculate convex hull for convex perimeter
                convperim = Convex_Perimeter(BW);

                [centers,radii] = Max_Fit_Circles(BW); %see function

                FL = Fiber_Length(BW);

                BM = Bleb_Length(BW,centers,radii);

                OneFrameStats=stats;
                OneFramePos=pos;
                OneFrameZ=z;
                OneFrameCellnums=cellnums;
                OneFrameCenters=centers;
                OneFrameRadii = radii;
                OneFrameFL = FL;
                OneFrameBM = BM;
                OneFrameConvPerims=convperim;
            else

                CC = bwconncomp(BW>0,4);
                stats = regionprops(CC, {'Area','Centroid',...
                    'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
                    'PerimeterOld','FilledArea',...
                    'MajorAxisLength','MinorAxisLength','Eccentricity',...
                    'Solidity','Extent','ConvexArea','PixelIdxList'});

                if ~isempty(stats)
                    pos = cat(1,stats.Centroid);
                    z = repmat(fpath_num,size(pos,1),1);
                    if ~isempty(OneFrameCellnums)
                        cellnums=[1:length(z)]'+max(OneFrameCellnums);%[cellnums;[1:length(z)]'+cellnums(end)];
                    else
                        cellnums=[1:length(z)]';
                    end
                    convperim = Convex_Perimeter(BW);
                    FL = Fiber_Length(BW);
                    [centers,radii] = Max_Fit_Circles(BW); %see function
                    BM = Bleb_Length(BW,centers,radii);

                    OneFrameStats=[OneFrameStats;stats];
                    OneFramePos=[OneFramePos;pos];
                    OneFrameZ=[OneFrameZ;z];
                    OneFrameCellnums=[OneFrameCellnums;cellnums];
                    OneFrameCenters = [OneFrameCenters; centers];
                    OneFrameFL = [OneFrameFL; FL];
                    OneFrameRadii = [OneFrameRadii; radii];
                    OneFrameBM = [OneFrameBM; BM];
                    OneFrameConvPerims=[OneFrameConvPerims; convperim];
                end
            end
        end
    end
    
    
    Igray=im2uint8(imread([gray_path,filesep,gray_ims(im_num).name]));
    if im_num>1
        maxcellnum=max(AllFrameNumbers);
        if ~isempty(prev_cellnums)
            [OneFrameCellnums,blackoutlist] = FindPreviousCell3D(BW, Igray,prev_pos,prev_cellnums,prev_stats,prev_z,OneFramePos,OneFrameCellnums,OneFrameStats,OneFrameZ,35,maxcellnum,blackoutlist);
            l=unique(OneFrameCellnums);
            for cnum=1:length(l)
                if sum(OneFrameCellnums==l(cnum))>1
                    disp('two or more cells display the same corrected cell number')
                end
            end
        else
            OneFrameCellnums=[1:length(OneFrameZ)]+maxcellnum;
        end
    end
    
    if OneFrameCellnums==0
        %this is the default return behavior of FindPreviousCell3D. In this
        %case, the user asked to exit and save the data.
        %prompt to user to ask what frame they want to save the data from,
        %then exit.
        input_hf = str2double(input('What frame should we save up to? (Enter an integer)   ', 's'));
        while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 
          input_hf = str2double(input('Invalid entry. What frame should we save up to? (Enter an integer)    ', 's'));
        end
        AllFrameNumbers=AllFrameNumbers(AllFrames<=input_hf);
        AllFrameStats=AllFrameStats(AllFrames<=input_hf);
        AllFramePos=AllFramePos(AllFrames<=input_hf,:);
        AllFrameZ=AllFrameZ(AllFrames<=input_hf);
        AllFrameCenters=AllFrameCenters(AllFrames<=input_hf,:);
        AllFrameConvPerims = AllFrameConvPerims(AllFrames<=input_hf);
        AllFrameFL = AllFrameFL(AllFrames<=input_hf);
        AllFrameRadii = AllFrameRadii(AllFrames<=input_hf);
        AllFrameBM = AllFrameBM(AllFrames<=input_hf);
        AllFrames=AllFrames(AllFrames<=input_hf);
        save(sprintf([save_path,filesep,'workspace_checkpoint_%d.mat'],input_hf),'AllFrameStats','AllFramePos','AllFrameZ','AllFrameCenters','AllFrameNumbers','AllFrames','AllFrameConvPerims','AllFrameFL','AllFrameRadii','AllFrameBM','blackoutlist')
        error('Exiting program, please check binary images if necessary')
    end
    
    %if there are new cells, we want to make sure the user identifies those
    %they do not want to follow.
    %how do you test if there are new cell numbers
    allcellnumspost=unique([OneFrameCellnums;allcellnums]);
    
    if check_bright_cells==1
        if length(allcellnumspost)~=length(allcellnums)
            %we have variables Igray, OneFrameCellnums, OneFrameStats, etc and
            %prev_cellnums, allcellnums,

            BW=zeros(size(Igray));
            T=zeros(size(Igray));
            blackoutpixels=zeros(size(Igray));
            bT=zeros(size(Igray));

            for obj = 1:length(OneFrameCellnums)
                if im_num>1
                    if sum(prev_cellnums==OneFrameCellnums(obj))==0
                        %that is if it is a new cell.
                        BW(OneFrameStats(obj).PixelIdxList)=1;
                        T(OneFrameStats(obj).PixelIdxList)=0.3;
                    end
                else
                    BW(OneFrameStats(obj).PixelIdxList)=1;
                    T(OneFrameStats(obj).PixelIdxList)=0.3;
                end
            end

            x=0;
            y=0;

            while(x>-1 && x<size(BW,2)+1 && y>-1 && y<size(BW,1)+1)
                disp('Please select green colored cells you wish to not be followed. If you may click the red cell to change it back to green. Green cells will be tracked.')
                figure(3000)
                imshow(Igray)
                hold on
                gdisp=imshow(cat(3,zeros(size(BW)),BW,zeros(size(BW))));
                set(gdisp,'AlphaData',T)
                hold on
                rdisp=imshow(cat(3,blackoutpixels,zeros(size(blackoutpixels)),zeros(size(blackoutpixels))));
                set(rdisp,'AlphaData',bT)
                title(sprintf('%s',gray_ims(im_num).name))

                [x,y]=ginput(1);
                x=round(x);
                y=round(y);
                if(x>-1 && x<size(BW,2)+1 && y>-1 && y<size(BW,1)+1)
                    lindex=sub2ind(size(BW),y,x); %form linear index

                    if BW(lindex)==1 || blackoutpixels(lindex)==1 %that is, if the selected point actually matters
                        %determine which cell number it is.
                        if im_num>1
                            possible_cellnums=setdiff(OneFrameCellnums,prev_cellnums); %this returns the values 
                        else
                            possible_cellnums=OneFrameCellnums;
                        end
                        %in OneFrameCellnums that are not in prev_cellnums
                        %determine which one of these cells were selected.
                        %inds is the cell number.
                        for i=1:length(possible_cellnums)
                            in=find(OneFrameCellnums==possible_cellnums(i)); %should only be length 1
                            if ~isempty(find(OneFrameStats(in).PixelIdxList==lindex))
                                %if the selected point is inside that cell
                                %add to blackoutlist, unless already on it. 
                                if sum(blackoutlist==possible_cellnums(i))<1
                                    %it is not in the list yet.
                                    blackoutlist=[blackoutlist,possible_cellnums(i)];
                                    BW(OneFrameStats(in).PixelIdxList)=0;
                                    %Turn display to red.
                                    T(OneFrameStats(in).PixelIdxList)=0;
                                    blackoutpixels(OneFrameStats(in).PixelIdxList)=1;
                                    bT(OneFrameStats(in).PixelIdxList)=0.3;
                                else
                                    %it is already in the list, delete it.
                                    blackoutlist(blackoutlist==possible_cellnums(i))=[];
                                    %Turn display back to green
                                    BW(OneFrameStats(in).PixelIdxList)=1;
                                    T(OneFrameStats(in).PixelIdxList)=0.3;
                                    blackoutpixels(OneFrameStats(in).PixelIdxList)=0;
                                    bT(OneFrameStats(in).PixelIdxList)=0;
                                end

                                break
                            end
                        end
                    else
                        if(x>-1 && x<size(BW,2)+1 && y>-1 && y<size(BW,1)+1)
                            disp('Invalid point.')
                        end
                    end
                end
            end
        end
    end
    
    allcellnums=OneFrameCellnums;
    
    %store for later
    prev_stats=OneFrameStats;
    prev_pos=OneFramePos;
    prev_z=OneFrameZ;
    prev_cellnums=OneFrameCellnums;
    prev_centers=OneFrameCenters;
    
    figure
    %Igray=imread([gray_path,filesep,gray_ims(im_num).name]);
    imshow(Igray)
    if ~isempty(OneFrameCellnums)
        for k=1:size(OneFramePos,1)
            if sum(blackoutlist==OneFrameCellnums(k))==0
                c=OneFrameCenters(k,:);
                text(c(1),c(2),sprintf('%d',OneFrameCellnums(k)),...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','middle',...
                    'FontWeight','bold',...
                    'Color','r');
            end
        end
    end
    axis off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=getframe;
    [framerows,framecols,~]=size(h.cdata);

    %%%FIXING OUTPUT SIZE ISSUE.
    if framerows==721
        im=h.cdata(2:end,:,:);
    else
        im=h.cdata;
    end
    %%%

    close all
    imwrite(im,[save_path,filesep,gray_ims(im_num).name])
    
    frames = repmat(im_num,length(OneFrameZ),1);
    %store all information for the current frame into larger complex
    if im_num==1
        AllFrameStats=OneFrameStats;
        AllFramePos=OneFramePos;
        AllFrameZ=OneFrameZ;
        AllFrameCenters=OneFrameCenters;
        AllFrameNumbers=OneFrameCellnums;
        AllFrames=frames;
        AllFrameConvPerims = OneFrameConvPerims;
        AllFrameFL = OneFrameFL;
        AllFrameRadii = OneFrameRadii;
        AllFrameBM = OneFrameBM;
    else
        AllFrameStats=[AllFrameStats; OneFrameStats];
        AllFramePos=[AllFramePos; OneFramePos];
        AllFrameZ=[AllFrameZ; OneFrameZ];
        AllFrameCenters=[AllFrameCenters; OneFrameCenters];
        AllFrameNumbers=[AllFrameNumbers; OneFrameCellnums];
        AllFrames=[AllFrames; frames];
        AllFrameConvPerims = [AllFrameConvPerims; OneFrameConvPerims];
        AllFrameFL = [AllFrameFL; OneFrameFL];
        AllFrameBM = [AllFrameBM; OneFrameBM];
        AllFrameRadii = [AllFrameRadii; OneFrameRadii];
    end
    
    if mod(im_num,10)==0 || im_num==1
        save(sprintf([save_path,filesep,'workspace_checkpoint_%d.mat'],im_num),'AllFrameStats','AllFramePos','AllFrameZ','AllFrameCenters','AllFrames','AllFrameNumbers','AllFrameConvPerims','AllFrameFL','AllFrameRadii','AllFrameBM','blackoutlist')
    end
end

%Delete the blackout cells we didn't want to track any more.
for obj=1:length(blackoutlist)
    AllFrames(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameZ(AllFrameNumbers==blackoutlist(obj))=[];
    AllFramePos(AllFrameNumbers==blackoutlist(obj),:)=[];
    AllFrameCenters(AllFrameNumbers==blackoutlist(obj),:)=[];
    AllFrameStats(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameConvPerims(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameFL(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameRadii(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameBM(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameNumbers(AllFrameNumbers==blackoutlist(obj))=[];
end

%save information in CSV File format. 
%frame number
%correctedNum
%area
%centroidx
%centroidy
%major
%minor
%eccentricity
%orientation
%convex area
%filledarea
%equivDiameter
%solidity
%extent
%perimeter
%perimeter old
%convex perimeter
%fiber length
%max in radii
%bleb length
%centersx
%centersy
A=cat(1,AllFrameStats.Area) * size_scale * size_scale;
posx=AllFramePos(:,1);
posy=AllFramePos(:,2);
mjr=cat(1,AllFrameStats.MajorAxisLength) * size_scale;
mnr=cat(1,AllFrameStats.MinorAxisLength) * size_scale;
ecc=cat(1,AllFrameStats.Eccentricity);
ori=cat(1,AllFrameStats.Orientation);
cA=cat(1,AllFrameStats.ConvexArea) * size_scale * size_scale;
fA=cat(1,AllFrameStats.FilledArea) * size_scale * size_scale;
eD=cat(1,AllFrameStats.EquivDiameter) * size_scale;
s=cat(1,AllFrameStats.Solidity);
ex=cat(1,AllFrameStats.Extent);
P=cat(1,AllFrameStats.Perimeter) * size_scale;
Pold=cat(1,AllFrameStats.PerimeterOld) * size_scale;
AllFrameConvPerims = AllFrameConvPerims * size_scale;
AllFrameFL = AllFrameFL * size_scale;
AllFrameRadii = AllFrameRadii * size_scale;
AllFrameBM = AllFrameBM * size_scale;
cx = AllFrameCenters(:,1);
cy = AllFrameCenters(:,2);



T=table(AllFrames,AllFrameNumbers, A, posx, posy, mjr, mnr, ecc, ori, ...
    cA, fA, eD, s, ex, P, Pold, AllFrameConvPerims, AllFrameFL, ...
    AllFrameRadii, AllFrameBM, cx, cy);

T.Properties.VariableNames={'Frames','CellNum','Area','CentroidX',...
    'CentroidY','MjrAxisLength','MnrAxisLength','Eccentricity',...
    'Orientation','ConvexArea','FilledArea','EquivDiameter',...
    'Solidity','Extent','Perimeter','PerimOld', 'ConvexPerim',...
    'FibLen', 'InscribeR','BlebLen','CenterX','CenterY'};

T=sortrows(T,2,'ascend');
idcs=strfind(save_path,'/');
writetable(T,[save_path(1:idcs(end)),save_path(idcs(end-1)+1:idcs(end)-1),'.csv'],'Delimiter',',','WriteRowNames',true)


%save the data just in case
save(sprintf([save_path,filesep,'workspace_checkpoint_%d.mat'],im_num),'AllFrameStats','AllFramePos','AllFrameZ','AllFrameCenters','AllFrames','AllFrameNumbers','AllFrameConvPerims','AllFrameFL','AllFrameRadii','AllFrameBM','blackoutlist')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define functions necessary for this algorithm

disp('Complete :)')
disp('Please run Detect_Errors.m and ChangeCellNum.m if necessary')

function [Centers,Rs] = Max_Fit_Circles(BW)
    B=bwboundaries(BW,4,'noholes'); %cell array
    [objs,~]=size(B);
    Centers = zeros(objs,2);
    Rs = zeros(objs,1);
    for bb=1:objs
        grain=false(size(BW));
        for uu=1:length(B{bb})
            grain(B{bb}(uu,1),B{bb}(uu,2))=true;
        end
        grain=uint8(grain)*255;
        [R,cx,cy] = max_inscribed_circle(grain,[]);
        Centers(bb,1)=cx;
        Centers(bb,2)=cy;
        Rs(bb) = R;
    end
end


function ConvP = Convex_Perimeter(BW)
    %the bw image should already be filtered by size.
    %from black and white, create concave objects.
    bw_labeled=bwlabel(BW,4);
    ConvP=zeros(max(max(bw_labeled)),1);
    %we might think that doing bwconvhull(bw,'objects') is indeed faster,
    %but it fails when two objects are nearby and when turned convex they
    %overlap.
    for i=1:max(max(bw_labeled))
        objs = bwconvhull(bw_labeled==i);
        %calculate the perimeter.
        CC = bwconncomp(objs>0,4);
        stats = regionprops(CC, {'Perimeter'});
        ConvP(i)=stats.Perimeter;
    end
end

function FL = Fiber_Length(BW)
    %the bw image should already be filtered by size and in correct binary
    %form.
    bw2=bwlabel(BW,4);
    FL = zeros(max(max(bw2)),1);
    for obj=1:max(max(bw2))
        bw1 = bw2==obj;
        %calculate longest skeleton path, aka fiber length for each object
        bw_skel=bwmorph(bw1,'skel','Inf');
        %how many objects are there?
        ep = find(bwmorph(bw_skel,'endpoints'));
        fl=0;
        for d=1:length(ep)
            D=bwdistgeodesic(bw_skel,ep(d),'quasi');
            if max(max(D))>fl
                fl=max(max(D));
            end
        end
        FL(obj) = fl;
    end
end

%function newcellnums=FindPreviousCell3D(oldpos,oldcellnums,oldstats,oldz,pos,cellnums,stats,z)
%    %use change in shape 

function max_bleb_len = Bleb_Length(BW,Centers,Rs)
    bwl = bwlabel(BW,4);
    max_bleb_len = zeros(max(max(bwl)),1);
    for obj=1:max(max(bwl))
        %find skeleton of obj
        bw = bwl==obj;
        bwskel = bwmorph(bw,'skel','Inf');
        Itemp = zeros(size(BW));
        Itemp(Centers(obj,2),Centers(obj,1))=1;
        Id = bwdist(Itemp);
        bw2 = bwskel;
        bw2(Id<ceil(Rs(obj)))=0;
        %find endpoints
        ep = find(bwmorph(bw2,'endpoints'));
        max_d = 0;
        for endpoint = 1:length(ep)
            bw3 = bwdistgeodesic(bw2,ep(endpoint),'quasi');
            m = nanmax(bw3(~isinf(bw3)));
            if m > max_d
                max_d = m; 
            end
        end
        max_bleb_len(obj) = max_d;
    end
    %     if max(max(bwl))==1
    %         C = imfuse(BW,bwskel,'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);
    %         imshow(C)
    %         hold on
    %         viscircles(Centers,Rs);
    %     else
    %         bwskel = bwmorph(BW,'skel','Inf');
    %         C = imfuse(BW,bwskel,'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);
    %         imshow(C)
    %         hold on
    %         viscircles(Centers,Rs);
    %     end
    %could also consider measuring bleb coherence, that is, the
    %directionality of all protrusions - random, in line or one sided.
end


function [newcellnums,blackoutlist_new] = FindPreviousCell3D(BW,Gray,oldpos,oldcellnums,oldstats,oldz,pos,cellnums,stats,z,maxdist,maxcellnum,blackoutlist)
    %use deltaArea, deltaR, and deltaZ positions to determine if cell
    %number needs to change.
    closestind = zeros(1,length(z));
    closest=zeros(1,length(z));
    newArea=cat(1,stats.Area);
    oldArea=cat(1,oldstats.Area);
    
    %determine the closest cell in the previous frame for each object in
    %the current frame. closest will have 0 whereever the closest cell of
    %the previous frame was not close enough.
    for obj=1:length(z)
        %determine closest proximity to cells in previous frame
        deltar = sqrt(((oldpos(:,1)-pos(obj,1)).^2) + ((oldpos(:,2)-pos(obj,2)).^2));
        %additionally, it is possible two cells at different depths are at
        %the same xy. we should check the areas, and make sure the change
        %in area is not too drastic. 
        deltaA = abs(oldArea-newArea(obj))./oldArea;
        [~,indA]=min(deltaA);
        [~,ind]=min(deltar); %tells you where the closest cell in the previous frame appeared. 
        %this index is for the previous frame stuff.
        if min(deltar)<maxdist %pixels, which is very, very large for a single step.
            if abs(oldz(ind)-z(obj))<2 %that is, most between two stacks.
                
                if indA==ind %easy choice, area and position make sense
                    closestind(obj)=ind;
                    closest(obj)=oldcellnums(ind);
                else
                    %in this case, you have to decide if position is enough
                    %to decide, or if Area is. My gut here says area, but
                    %could be wrong..
                    try
                        closestind(obj)=ind;
                        closest(obj)=oldcellnums(ind);
                        %if deltar(indA)<25 & deltaA(indA)<deltaA(ind)
                        %    %if it both close in proximity and area makes sense
                        %    closestind(obj)=indA;
                        %    closest(obj)=oldcellnums(indA);
                        %else
                        %    %if area makes sense but it is far away
                        %    closestind(obj)=ind;
                        %    closest(obj)=oldcellnums(ind); %should this be zero?
                        %end
                    catch
                        disp('pause')
                    end
                end
            end
        end
    end
    
    %so now, if 2 or more cells say their closest neighbor in the prvious frame is the same
    %cell, which one is correct? 
    %outter loop, may not be necessary. run 3 times?
    outterloop=0;
    while length(nonzeros(unique(closest)))<length(nonzeros(closest))
        outterloop=outterloop+1;

        if ~isempty(nonzeros(closest)) %that is, there are not completely new cells in the next frame
            %closest: vector with length equal to number of objects in current
            %image. Values are 0 if there is no closest neighbor (<20) or has
            %integer number corresponding to closest cell in previous image.
            %closestind: Same * but number is the index where in oldpos,oldstats,etc 
            % the cell was last seen.
            u_closest = nonzeros(unique(closest));%form unique list
            for i=1:length(u_closest)
                if sum(closest==u_closest(i))>1 %two or more cells share neighbor
                    current_cells_ind=find(closest==u_closest(i));
                    %calculate their delta r's
                    if length(current_cells_ind)>2
                        %3 or more cells point to the same one. Hopefully not.
                        disp('pause')
                    else
                        min1=zeros(1,2);
                        min1ind=zeros(1,2);
                        min2=zeros(1,2);
                        min2ind=zeros(1,2);
                        for c_num=1:2
                            deltar=sqrt(((oldpos(:,1)-pos(current_cells_ind(c_num),1)).^2)+...
                                ((oldpos(:,2)-pos(current_cells_ind(c_num),2)).^2));
                            [vs,inds]=sort(deltar);
                            min1(c_num)=vs(1); %distance to the very closest cell (which should be the same for both).
                            min2(c_num)=vs(2); %the distance to the second closest object
                            min1ind(c_num)=inds(1); %index in oldpos, etc where min distance was
                            min2ind(c_num)=inds(2);
                        end
                        
                        
                        if sum(min2>30)>1 % both share only the closest cell in the previous frame; there is no second closer neighbor
                            for c_num=1:2
                                closest(current_cells_ind(c_num))=0;
                            end
                        elseif sum(min2>maxdist)==1 % only one has no other close neighbor, meaning it is likely the correct cell and the other is wrong
                            %change the one that has min2<25
                            closest(current_cells_ind(min2<maxdist))=oldcellnums(min2ind(min2<maxdist));
                            %note, this might point a cell to another that has
                            %already been tested, so may want to rerun this
                            %entire loop. 
                        else
                            %This case is where both cells point to the
                            %same closest cell and the same second closest
                            %cell. In this case, we want to check if cell A
                            %is closer to both or if cell B is closer to
                            %both. In these cases, we want to check deltaA
                            %to tell us which one it should take. 
                            [~,c1]=min(min1);
                            [~,c2]=min(min2);
                            if c1==c2 %this case is a weird one
                                %UI?
                                %ask user to identify which cell is which.
                                
                                for c_num=1:2
                                    figure(2000)
                                    BW_temp=zeros(size(BW));
                                    BW_temp(stats(current_cells_ind(c_num)).PixelIdxList)=1;
                                    I_temp=imfuse(Gray,BW_temp>0,'ColorChannels',[2,1,0]);
                                    imshow(I_temp)
                                    title('Which old cell number is it?')
                                    input_hf = str2double(input(sprintf(['Is this cell %d, %d, or a different one? (enter an integer)    '],oldcellnums(min1ind(c_num)),oldcellnums(min2ind(c_num))), 's'));
                                    while isnan(input_hf) || fix(input_hf) ~= input_hf 
                                      input_hf = str2double(input(sprintf(['Is this cell %d, %d, or a different one? (enter an integer)    '], oldcellnums(min1ind(c_num)), oldcellnums(min2ind(c_num))) , 's'));
                                    end
                                    closest(current_cells_ind(c_num))=input_hf;
                                    close(figure(2000))
                                end
                                
                            else
                                if min2(2)<min2(1)
                                    closest(current_cells_ind(2))=oldcellnums(min2ind(2));
                                else
                                    closest(current_cells_ind(1))=oldcellnums(min2ind(1));
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if outterloop>20
            disp('Engaging User')
            %would be good to put a while loop here to make sure all will
            %be different. 
            while length(nonzeros(unique(closest)))<length(nonzeros(closest))
                %record which cells it is confused by? In this case,
                %we might be stuck here. The problem can occur if one cell has
                %moved such that is is closer to two other objects which
                %than its old position. Typically, you will see three cells
                %with deltar<0. 
                %here, we should create a UI to ask the user to distinuish
                %between the confused cells?
                %ask what number this cell is. If new, add to blackoutlist.
                %because we didn't choose to follow it.
                %display the two cells that point to the same one.
                %I would suggest using transparency.

                nonzeroclosest=nonzeros(closest);
                %if length(unique(nonzeroclosest))==length(nonzeroclosest)-1
                %inds of nonzero elements
                inds_nonzero=find(closest~=0);
                [~,ind]=unique(nonzeroclosest);
                duplicate_ind=setdiff(1:length(nonzeroclosest),ind);
                duplicate_val = nonzeroclosest(duplicate_ind);

                for v=1:length(duplicate_val)
                    val=duplicate_val(v);
                    inds=find(closest==val);
                    if length(inds)==2
                        %these should be the right indeces now... 
                        %we want to display these 

                        %lets find where the two indeces repeat.
                        figure(4000)
                        imshow(Gray)
                        hold on
                        gbw=zeros(size(Gray));
                        gT=zeros(size(Gray));
                        gbw(stats(inds(1)).PixelIdxList)=1;
                        gT(stats(inds(1)).PixelIdxList)=0.3;
                        rbw=zeros(size(Gray));
                        rT=zeros(size(Gray));
                        rbw(stats(inds(2)).PixelIdxList)=1;
                        rT(stats(inds(2)).PixelIdxList)=0.3;
                        gch=imshow(cat(3,zeros(size(gbw)),gbw,zeros(size(gbw))));
                        set(gch,'AlphaData',gT)
                        rch=imshow(cat(3,rbw,zeros(size(rbw)),zeros(size(rbw))));
                        set(rch,'AlphaData',rT)

                        %ask user questions now.
                        input_hf = str2double(input('Which cell number is the green cell? (Enter an integer, >9999 for new and ignore)    ', 's'));
                        while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 
                          input_hf = str2double(input('Invalid entry. Which cell number is the green cell? (Enter an integer, >9999 for new and ignore)    ', 's'));
                        end

                        closest(inds(1))=input_hf;

                        input_hf = str2double(input('Which cell number is the red cell? (Enter an integer, >9999 for new and ignore)    ', 's'));
                        while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 
                          input_hf = str2double(input('Invalid entry. Which cell number is the red cell? (Enter an integer, >9999 for new and ignore)    ', 's'));
                        end

                        closest(inds(2))=input_hf;
                    elseif length(inds)==1
                        figure(4000)
                        imshow(Gray)
                        hold on
                        gbw=zeros(size(Gray));
                        gT=zeros(size(Gray));
                        gbw(stats(inds(1)).PixelIdxList)=1;
                        gT(stats(inds(1)).PixelIdxList)=0.3;
                        gch=imshow(cat(3,zeros(size(gbw)),gbw,zeros(size(gbw))));
                        set(gch,'AlphaData',gT)
                        %ask user questions now.
                        input_hf = str2double(input('Which cell number is the green cell? (Enter an integer)    ', 's'));
                        while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 
                          input_hf = str2double(input('Invalid entry. Which cell number is the green cell? (Enter an integer)    ', 's'));
                        end

                        closest(inds(1))=input_hf;
                    elseif length(inds)>2
                        disp('Three or more cells point to the same old point')
                        for el=1:length(inds)
                            figure(4000)
                            imshow(Gray)
                            hold on
                            gbw=zeros(size(Gray));
                            gT=zeros(size(Gray));
                            gbw(stats(inds(el)).PixelIdxList)=1;
                            gT(stats(inds(el)).PixelIdxList)=0.3;
                            gch=imshow(cat(3,zeros(size(gbw)),gbw,zeros(size(gbw))));
                            set(gch,'AlphaData',gT)
                            %ask user questions now.
                            input_hf = str2double(input('Which cell number is the green cell? (Enter an integer)    ', 's'));
                            while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 
                              input_hf = str2double(input('Invalid entry. Which cell number is the green cell? (Enter an integer)    ', 's'));
                            end

                            closest(inds(el))=input_hf;
                        end
                    end
                end
            end

            close(figure(4000))
            %end
            
            disp('There might be an error with the binary images. Please review the saved files.')
            %ask user questions now.
            input_hf = str2double(input('Do you wish to save the current data and exit? (Enter an integer) 1: yes  2: no    ', 's'));
            while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1  || input_hf>2
              input_hf = str2double(input('Invalid entry. Do you wish to save the current data and exit? (Enter an integer) 1: yes  2: no    ', 's'));
            end
            
            if input_hf==1
                newcellnums=0; %easy default return case
                return 
            end
            
        end
    end

    inds_new = find(closest>maxcellnum & closest<=9999);
    inds_new_bad=find(closest>9999);
    closest(closest>9999)=0;
    newcellnums=closest;
    inds=find(closest>0);
    newcellnums(inds)=closest(inds);
    inds=find(closest==0);
    newcellnums(inds)=[1:length(inds)]+maxcellnum;
    maxcellnum2=max(newcellnums);
    newcellnums(inds_new)=[1:length(inds_new)]+maxcellnum2;
    blackoutlist_new=[blackoutlist, newcellnums(inds_new_bad)];
    newcellnums=newcellnums';
end
