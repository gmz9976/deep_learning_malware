function [gist, param] = fl_gist(D, HOMEIMAGES, param, HOMEGIST)
%
% [gist, param] = LMgist(D, HOMEIMAGES, param);
% [gist, param] = LMgist(filename, HOMEIMAGES, param);
% [gist, param] = LMgist(filename, HOMEIMAGES, param, HOMEGIST);
%
% For a set of images:
% gist = LMgist(img, [], param);
%
% When calling LMgist with a fourth argument it will store the gists in a
% new folder structure mirroring the folder structure of the images. Then,
% when called again, if the gist files already exist, it will just read
% them without recomputing them:
%
%   [gist, param] = LMgist(filename, HOMEIMAGES, param, HOMEGIST);
%   [gist, param] = LMgist(D, HOMEIMAGES, param, HOMEGIST);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling the shape of the scene: a holistic representation of the spatial envelope
% Aude Oliva, Antonio Torralba
% International Journal of Computer Vision, Vol. 42(3): 145-175, 2001.

if nargin==4
    precomputed = 1;
    % get list of folders and create non-existing ones
    %listoffolders = {D(:).annotation.folder};

    %for i = 1:length(D);
    %    f{i} = D(i).annotation.folder;
    %end
    %[categories,b,class] = unique(f);
else
    precomputed = 0;
    HOMEGIST = '';
end

% select type of input
if isstruct(D)
    % [gist, param] = LMgist(D, HOMEIMAGES, param);
    Nscenes = length(D);
    typeD = 1;
end
if iscell(D)
    % [gist, param] = LMgist(filename, HOMEIMAGES, param);
    Nscenes = length(D);
    typeD = 2;
end
if isnumeric(D)
    % [gist, param] = LMgist(img, HOMEIMAGES, param);
    Nscenes = size(D,4);
    typeD = 3;
    if ~isfield(param, 'imageSize')
        param.imageSize = [size(D,1) size(D,2)];
    end
end
    
param.boundaryExtension = 32; % number of pixels to pad

if nargin<3
    % Default parameters
    param.imageSize = 128;
    param.orientationsPerScale = [8 8 8 8];
    param.numberBlocks = 4;
    param.fc_prefilt = 4;
    param.G = createGabor(param.orientationsPerScale, param.imageSize+2*param.boundaryExtension);
else
    if ~isfield(param, 'G')
        param.G = createGabor(param.orientationsPerScale, param.imageSize+2*param.boundaryExtension);
    end
end

% Precompute filter transfert functions (only need to do this once, unless
% image size is changes):
Nfeatures = size(param.G,3)*param.numberBlocks^2;


% Loop: Compute gist features for all scenes
gist = zeros([Nscenes Nfeatures], 'single');
for n = 1:Nscenes
    g = [];
    todo = 1;
    
    % if gist has already been computed, just read the file
    if precomputed==1
        filegist = fullfile(HOMEGIST, D(n).annotation.folder, [D(n).annotation.filename(1:end-4) '.mat']);
        if exist(filegist, 'file')
            load(filegist, 'g');
            todo = 0;
        end
    end
    
    % otherwise compute gist
    if todo==1
        if Nscenes>1 disp([n Nscenes]); end

        % load image
        try
            switch typeD
                case 1
                    img = LMimread(D, n, HOMEIMAGES);
                case 2
                    img = imread(fullfile(HOMEIMAGES, D{n}));
                case 3
                    img = D(:,:,:,n);
            end
        catch
            disp(D(n).annotation.folder)
            disp(D(n).annotation.filename)
            rethrow(lasterror)
        end
        
        % convert to gray scale
        img = single(mean(img,3));

        % resize and crop image to make it square
%         img = imresizecrop(img, param.imageSize, 'bilinear');
        img = imresize(img, param.imageSize, 'bilinear'); %jhhays

        % scale intensities to be in the range [0 255]
        img = img-min(img(:));
        img = 255*img/max(img(:));
        
        if Nscenes>1
            imshow(uint8(img))
            title(n)
        end

        % prefiltering: local contrast scaling
        output    = prefilt(img, param.fc_prefilt);

        % get gist:
        g = gistGabor(output, param);
        
        % save gist if a HOMEGIST file is provided
        if precomputed
            mkdir(fullfile(HOMEGIST, D(n).annotation.folder))
            save (filegist, 'g')
        end
    end
    
    gist(n,:) = g;
    drawnow
end


end

