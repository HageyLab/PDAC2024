%% Fiber Quantification Algorithm v3
% Created by Shamik Mascharak
% November 28th, 2019 - Created
% December 1st, 2019 - First iteration completed
% October 28th, 2020 - Uploaded to Github
% May 7th, 2021 - Code updated (Jason Guo)

%% Read files and generate binarized images - RUN THIS FIRST

clear
clc
close all

% Get list of all .tif files in the directory
imagefiles = dir('*.tif');
% imagefiles = dir('*.png');
nfiles = length(imagefiles);

%% Run color deconvolution

counter = 0;

% Optional scaling factors (to convert pixels to microns)
scale = 1;

% Define target image for image normalization
target_im = imread(imagefiles(234540).name); % 16_1.tif_7254_41106.tif

verbose = 0; % 1 to display results, 0 to not

% Define stain vector for color deconvolution 
he_vect = [0.6443186 0.7166757 0.26688856;...
    0.09283128 0.9545457 0.28324;...
    0.63595444 0.001 0.7717266];

% Define warning counters for images that are blank/empty after processing
warn_count_pink = 0;
warn_count_pink_bw = 0;
warn_count_pink_bw_branch = 0;
warn_count_pink_bw_skel = 0;
warn_count_pink_bw_strel = 0;
warn_count_pink_noise = 0;

for i = 1:nfiles
    filename = imagefiles(i).name;
    stored_names{i} = filename;
    stored_im{i} = imread(filename); % Store the original image
    
    % Normalize image using RGB histogram method
   stored_norm_im{i} = Norm(stored_im{i},target_im,'RGBHist',verbose);
%   stored_norm_im{i} = stored_im{i}; % To bypass normalization

    % Perform color deconvolution (Ruifrok method)
    [Dch M] = Deconvolve(stored_norm_im{i},he_vect,0);
    
    [purple pink green] = PseudoColourStains(Dch,M);
    
    he_pink = (255-rgb2gray(pink))-(255-rgb2gray(green));
    
    % Remove noise with adaptive filtering
    he_pink_noise = wiener2(he_pink,[3 3]);
    
    stored_he_pink_noise{i} = he_pink_noise;
    
    he_pink_bw = im2bw(imadjust(he_pink_noise));
    
    % Define diamond structuring element
    seD = strel('diamond',1);
    he_pink_bw_strel = imerode(he_pink_bw,seD);
    
    stored_he_pink_strel{i} = he_pink_bw_strel; % Binary map
    
    % Skeletonize, separate at branchpoints, and delete single pixel lines
    he_pink_bw_skel = bwmorph(he_pink_bw_strel,'skel',Inf);
    he_pink_bw_branch = bwmorph(he_pink_bw_skel,'branchpoints');
    he_pink_bw_skel(he_pink_bw_branch == 1) = 0;
    he_pink_bw_skel = bwmorph(he_pink_bw_skel,'clean');
    
    stored_he_pink_skel{i} = he_pink_bw_skel; % Skeletonized image
    stored_he_pink_branch{i} = he_pink_bw_branch;
    
    % Store indices of images that are blank/empty after processing. These
    % indices can be found in the variables containing the 'warnings_' prefix.
    if not(any(any(he_pink)))
        warnings_he_pink(warn_count_pink+1) = i;
        warn_count_pink = warn_count_pink + 1;
    end
    
    if not(any(any(he_pink_bw)))
        warnings_he_pink_bw(warn_count_pink_bw+1) = i;
        warn_count_pink_bw = warn_count_pink_bw + 1;
    end

    if not(any(any(he_pink_bw_branch)))
        warnings_he_pink_bw_branch(warn_count_pink_bw_branch+1) = i;
        warn_count_pink_bw_branch = warn_count_pink_bw_branch + 1;
    end
    
    if not(any(any(he_pink_bw_skel)))
        warnings_he_pink_bw_skel(warn_count_pink_bw_skel+1) = i;
        warn_count_pink_bw_skel = warn_count_pink_bw_skel + 1;
    end   
    
    if not(any(any(he_pink_bw_strel)))
        warnings_he_pink_bw_strel(warn_count_pink_bw_strel+1) = i;
        warn_count_pink_bw_strel = warn_count_pink_bw_strel + 1;
    end    
    
    if not(any(any(he_pink_noise)))
        warnings_he_pink_noise(warn_count_pink_noise+1) = i;
        warn_count_pink_noise = warn_count_pink_noise + 1;
    end
    
    % Display progress
    counter = counter + 1;
    progress = 100*counter/(nfiles)
    
end

stored_names = stored_names';

%% Exclude error-triggering images from current analysis
% Identify all unique error-triggering images
total_warnings = double.empty;
if (warn_count_pink > 0)
    total_warnings = horzcat(total_warnings,warnings_he_pink);
end
if (warn_count_pink_bw > 0)
    total_warnings = horzcat(total_warnings,warnings_he_pink_bw);
end
if (warn_count_pink_bw_branch > 0)
    total_warnings = horzcat(total_warnings,warnings_he_pink_bw_branch);
end
if (warn_count_pink_bw_skel > 0)
    total_warnings = horzcat(total_warnings,warnings_he_pink_bw_skel);
end
if (warn_count_pink_bw_strel > 0)
    total_warnings = horzcat(total_warnings,warnings_he_pink_bw_strel);
end
if (warn_count_pink_noise > 0)
    total_warnings = horzcat(total_warnings,warnings_he_pink_noise);
end
total_warnings = unique(total_warnings);

% Remove error-triggering images from analysis
original_imagefiles = imagefiles;
total_warn_count = length(total_warnings);
nfiles = nfiles-total_warn_count;
stored_he_pink_noise = stored_he_pink_noise(~ismember(1:length(stored_he_pink_noise),total_warnings));
stored_he_pink_strel = stored_he_pink_strel(~ismember(1:length(stored_he_pink_strel),total_warnings));
stored_he_pink_skel = stored_he_pink_skel(~ismember(1:length(stored_he_pink_skel),total_warnings));
stored_he_pink_branch = stored_he_pink_branch(~ismember(1:length(stored_he_pink_branch),total_warnings));
stored_names = stored_names(~ismember(1:length(stored_names),total_warnings));
stored_im = stored_im(~ismember(1:length(stored_im),total_warnings));
stored_norm_im = stored_norm_im(~ismember(1:length(stored_norm_im),total_warnings));
imagefiles = imagefiles(~ismember(1:length(imagefiles),total_warnings));

%% Quantify image features

clc
counter = 0;

for i = 1:nfiles
    
    pink = stored_he_pink_noise{i};
    pink_bw = stored_he_pink_strel{i};
    pink_skel = stored_he_pink_skel{i};
    pink_branch = stored_he_pink_branch{i};
    
    % Calculate values from grayscale images
    props_pink = regionprops(pink_bw,pink,'all');
    % Fiber angle-related values
    alpha_pink = [props_pink.Orientation]*pi/180;
    circ_pink = [mean(alpha_pink) median(alpha_pink) std(alpha_pink) skewness(alpha_pink)...
        kurtosis(alpha_pink) circ_kappa(alpha_pink)];

    % Calculate Haralick features from grayscale images (4 different offsets)
    offset = [0 1; -1 1;-1 0;-1 -1];
    glcm_pink = graycomatrix(pink,'Offset',offset,'Symmetric',true);
    
    grayprops_pink = graycoprops(glcm_pink);
    
    % Calculate values from binary images
    props_pink_bw = regionprops(pink_bw,'all');
    % Fiber angle-related values
    alpha_pink_bw = [props_pink_bw.Orientation]*pi/180;
    circ_pink_bw = [mean(alpha_pink_bw) median(alpha_pink_bw) std(alpha_pink_bw) skewness(alpha_pink_bw)...
        kurtosis(alpha_pink_bw) circ_kappa(alpha_pink_bw)];

    % Calculate values from skeletonized images
    props_pink_skel = regionprops(pink_skel,'all');
    % Fiber angle-related values
    alpha_pink_skel = [props_pink_skel.Orientation]*pi/180;
    circ_pink_skel = [mean(alpha_pink_skel) median(alpha_pink_skel) std(alpha_pink_skel) skewness(alpha_pink_skel)...
        kurtosis(alpha_pink_skel) circ_kappa(alpha_pink_skel)];
    
    % Branchpoint-related values
    [L_pink, num_pink] = bwlabel(pink_skel);
    [~, num_pink_branch] = bwlabel(pink_branch);
    
    % Grayscale regionprops, grayscale circ, Haralick features, binary regionprops,
    % binary circ, skeletonized regionprops, skeletonized circ, branchpoints,
    
    quantified{i} = [mean([props_pink.Area]) std([props_pink.Area])...
        mean([props_pink.MajorAxisLength]) std([props_pink.MajorAxisLength])...
        mean([props_pink.MinorAxisLength]) std([props_pink.MinorAxisLength])...
        mean([props_pink.Eccentricity]) std([props_pink.Eccentricity])...
        mean([props_pink.ConvexArea]) std([props_pink.ConvexArea])...
        mean([props_pink.Circularity]~=Inf) std([props_pink.Circularity]~=Inf)...
        mean([props_pink.FilledArea]) std([props_pink.FilledArea])...
        mean([props_pink.EulerNumber]) std([props_pink.EulerNumber])...
        sum([props_pink.EulerNumber])...
        mean([props_pink.EquivDiameter]) std([props_pink.EquivDiameter])...
        mean([props_pink.Solidity]) std([props_pink.Solidity])...
        mean([props_pink.Extent]) std([props_pink.Extent])...
        mean([props_pink.Perimeter]) std([props_pink.Perimeter])...
        mean([props_pink.PerimeterOld]) std([props_pink.PerimeterOld])...
        mean([props_pink.MeanIntensity]) std([props_pink.MeanIntensity])...
        mean([props_pink.MinIntensity]) std(double([props_pink.MinIntensity]))...
        mean([props_pink.MaxIntensity]) std(double([props_pink.MaxIntensity]))...
        mean([props_pink.MaxFeretDiameter]) std([props_pink.MaxFeretDiameter])...
        mean([props_pink.MaxFeretAngle]) std([props_pink.MaxFeretAngle])...
        mean([props_pink.MinFeretDiameter]) std([props_pink.MinFeretDiameter])...
        mean([props_pink.MinFeretAngle]) std([props_pink.MinFeretAngle])...
        circ_pink...
        [grayprops_pink.Contrast] [grayprops_pink.Correlation] [grayprops_pink.Energy] [grayprops_pink.Homogeneity]...
        mean([props_pink_bw.Area]) std([props_pink_bw.Area])...
        mean([props_pink_bw.MajorAxisLength]) std([props_pink_bw.MajorAxisLength])...
        mean([props_pink_bw.MinorAxisLength]) std([props_pink_bw.MinorAxisLength])...
        mean([props_pink_bw.Eccentricity]) std([props_pink_bw.Eccentricity])...
        mean([props_pink_bw.ConvexArea]) std([props_pink_bw.ConvexArea])...
        mean([props_pink_bw.Circularity]~=Inf) std([props_pink_bw.Circularity]~=Inf)...
        mean([props_pink_bw.FilledArea]) std([props_pink_bw.FilledArea])...
        mean([props_pink_bw.EulerNumber]) std([props_pink_bw.EulerNumber])...
        sum([props_pink_bw.EulerNumber])...
        mean([props_pink_bw.EquivDiameter]) std([props_pink_bw.EquivDiameter])...
        mean([props_pink_bw.Solidity]) std([props_pink_bw.Solidity])...
        mean([props_pink_bw.Extent]) std([props_pink_bw.Extent])...
        mean([props_pink_bw.Perimeter]) std([props_pink_bw.Perimeter])...
        mean([props_pink_bw.PerimeterOld]) std([props_pink_bw.PerimeterOld])...
        mean([props_pink_bw.MaxFeretDiameter]) std([props_pink_bw.MaxFeretDiameter])...
        mean([props_pink_bw.MaxFeretAngle]) std([props_pink_bw.MaxFeretAngle])...
        mean([props_pink_bw.MinFeretDiameter]) std([props_pink_bw.MinFeretDiameter])...
        mean([props_pink_bw.MinFeretAngle]) std([props_pink_bw.MinFeretAngle])...
        circ_pink_bw...
        mean([props_pink_skel.Area]) std([props_pink_skel.Area])...
        mean([props_pink_skel.MajorAxisLength]) std([props_pink_skel.MajorAxisLength])...
        mean([props_pink_skel.MinorAxisLength]) std([props_pink_skel.MinorAxisLength])...
        mean([props_pink_skel.Eccentricity]) std([props_pink_skel.Eccentricity])...
        mean([props_pink_skel.ConvexArea]) std([props_pink_skel.ConvexArea])...
        mean([props_pink_skel.Circularity]~=Inf) std([props_pink_skel.Circularity]~=Inf)...
        mean([props_pink_skel.FilledArea]) std([props_pink_skel.FilledArea])...
        mean([props_pink_skel.EulerNumber]) std([props_pink_skel.EulerNumber])...
        sum([props_pink_skel.EulerNumber])...
        mean([props_pink_skel.EquivDiameter]) std([props_pink_skel.EquivDiameter])...
        mean([props_pink_skel.Solidity]) std([props_pink_skel.Solidity])...
        mean([props_pink_skel.Extent]) std([props_pink_skel.Extent])...
        mean([props_pink_skel.Perimeter]) std([props_pink_skel.Perimeter])...
        mean([props_pink_skel.PerimeterOld]) std([props_pink_skel.PerimeterOld])...
        mean([props_pink_skel.MaxFeretDiameter]) std([props_pink_skel.MaxFeretDiameter])...
        mean([props_pink_skel.MaxFeretAngle]) std([props_pink_skel.MaxFeretAngle])...
        mean([props_pink_skel.MinFeretDiameter]) std([props_pink_skel.MinFeretDiameter])...
        mean([props_pink_skel.MinFeretAngle]) std([props_pink_skel.MinFeretAngle])...
        circ_pink_skel...
        num_pink num_pink_branch];

    counter = counter + 1;
    progress = 100*counter/(nfiles)
    
end

quantified = cell2mat(quantified');

%% UMAP
quantified_with_dummy_label = [ones(1,147);quantified];
writematrix(quantified_with_dummy_label,'Quantified.csv');
umap = run_umap('Quantified.csv','save_template_file','umaptemplate.mat');
