% function for tracking centroid of larva ---------------------------------

function comLarvae = tracking(dataFolder,filesArray,maxFile,larvalMinSize,larvalMaxSize,scale,thLevel,numLarvae,chamTemp,ibg)
 comLarvae = nan(2,maxFile); % initialization of variable
 for file = 1:maxFile
    
    comTemp = [];
    fileName = filesArray(file).name;
    iTemp = imread(strcat(dataFolder,fileName));   % temporary image
    iTemp = iTemp(:,:,1);   % read only one channel of the image
    h = fspecial('gaussian',5, 0.7);  % creating a Gaussian filter
    iSub = ibg - iTemp; % background substraction    
    iFiltTemp = imfilter(iSub,h,'replicate');  % applying the filter
    %iCropTemp = imcrop(iFiltTemp,chamTemp); % crop the image with size of each chamber
    iFilt = imbinarize(iFiltTemp); % convert image to binary image by automatic thresholding 
    iCropTemp2 = imcrop(iFilt,chamTemp);
    %iCropTemp2 = imbinarize(iCropTemp,thLevel); % convert image to binary image by automatic thresholding 
    iCrop = bwareaopen(iCropTemp2,larvalMinSize); % eliminate the smaller objects than mininum size of larva

    L = bwlabel(iCrop,8);
    %infoL = regionprops(L,'Centroid','Area'); % get the information of centroid and area for the objects
    infoL = regionprops(L,'Centroid','Area','MajorAxisLength'); % get the information of centroid, area, and length of the major axis for the objects
    

     if (isempty(infoL) == 1)
       % mark the centroid of larva in case of no object-detection
       figure, imshow(iTemp)
       hold on
       title(sprintf('Mark the centroid of #%i animal manually by clicking (small larva).',numLarvae))
       comTemp = ginput(1)';
       comLarvae(1,file) = comTemp(1);
       comLarvae(2,file) = comTemp(2);
       plot(comTemp(1),comTemp(2),'+r')
       pause(0.5)
       close all
      else       
    % for normal object detection, find the largest object
        areas = cat(1,infoL.Area);
        majorAxisLength = cat(1,infoL.MajorAxisLength);
        coms = cat(1,infoL.Centroid);
        idxMaxArea = find(areas == max(areas),1);
        comTemp = coms(idxMaxArea,:)';   % local x-z coordinates of the centroid
    % global x-z coordinates of the position of the centroid
        comLarvae(1,file) = chamTemp(1) + comTemp(1);  
        comLarvae(2,file) = chamTemp(2) + comTemp(2);
        
    % indicate the centroid manually if the object is too large
        if (areas(idxMaxArea)>larvalMaxSize || majorAxisLength(idxMaxArea)>larvalMaxSize*scale)
    %if (majorAxisLength(idxMaxArea)>larvalMaxSize)
          figure, imshow(iTemp)
          hold on
          title(sprintf('Mark the centroid of #%i animal manually by clicking (big larva).',numLarvae))
          comTemp = ginput(1)';
          comLarvae(1,file) = comTemp(1);
          comLarvae(2,file) = comTemp(2);
          plot(comTemp(1),comTemp(2),'+r')
          pause(0.5)
          close all
        end
     end

 end
end