%% Script for cropping images %%
% Developed by Daeyeon Kim

clear
close all

%% SETTING

% SETTING FOLDERS -----------------------------------------------------------------
% setting the main datafolder
%cd('/Volumes/Macintosh HD/PhD/Data/Dig_n_Dive/')

dataDir0='/Users/daeyeonkim/Documents/Box Sync/MATLAB/projects/DnD/data/';

% setting the sub-folder; dont forget the slash to go into the folder
dataFolder='/CS_nocue_agar04_1/';

dataDir=strcat(dataDir0,dataFolder);
outDir1=strcat(dataDir0,'sample_upper/');
outDir2=strcat(dataDir0,'sample_lower/');

%filesArray=dir(strcat(dataDir,'*agar04_1*.bmp'));
filesArray=dir(strcat(dataDir,'*.jpg'));

% -------------------------------------------------------------------------

% SETTGING FOR IMAGE FRAMES TO BE ANALYZED --------------------------------

fps = 1/1;  % Time scale; one frame every 1 seconds
tanal = 120;   % time for analylzing in min
maxFile=60*fps*tanal;  % maximum frames to be analyzed 
%maxFile=10;
% restricting the frames to the maxFile
if length(filesArray)>maxFile
    maxFile=maxFile;
else
    maxFile=length(filesArray);
end

cd(dataDir)

filename=filesArray(1).name;
i0=imread(filename);
ibg=double(i0(:,:,1));%;/maxFile;   % gray image
ibg=uint8(ibg);

H=size(ibg,1);
W=size(ibg,2);
figure, imshow(ibg)
%display('click to crop the middle of image')
title('CROP: Draw the area of the upper device for the analysis.')

h = imrect;
roi_upper = round(getPosition(h));  % [xmin ymin width height]

%display('click to crop the middle of image')
title('CROP: Draw the area of the lower device for the analysis.')

h = imrect;
roi_lower = round(getPosition(h));  % [xmin ymin width height]

iniFile=1;
finalFile=maxFile;

for i=1:2
    if (i==1)
        mkdir(outDir1)
        for file=iniFile:finalFile
            file/finalFile
            filename=filesArray(file).name;
            i0=imread(filename);
            ibg=double(i0(:,:,1));  % gray image
            ibg=uint8(ibg);
            ibgcrop=imcrop(ibg,roi_upper);
            
            %imshow(ibgcrop);
            outfile = strcat(outDir1,filename);
            %print(gcf,'-dbmp',outfile);
            %imwrite(ibgcrop,outfile,'bmp')
            imwrite(ibgcrop,outfile,'jpg')
            %saveas(f1,strcat(outDir1,filename),'bmp');
        end
    else
        mkdir(outDir2)
        for file=iniFile:finalFile
            file/finalFile
            filename=filesArray(file).name;
            i0=imread(filename);
            ibg=double(i0(:,:,1));  % gray image
            ibg=uint8(ibg);
            ibgcrop=imcrop(ibg,roi_lower);
            %imshow(ibgcrop);
            outfile = strcat(outDir2,filename);
            %print(gcf,'-dbmp',outfile);
            %imwrite(ibgcrop,outfile,'bmp')
            imwrite(ibgcrop,outfile,'jpg')
        end        
    end
end
