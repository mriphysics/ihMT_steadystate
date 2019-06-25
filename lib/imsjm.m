%***********************************************************
% imsjm - imagesc tool with a few better features. 18-5-05
%   % 11-8-08: beefed up version
%   % 10-9-08:  use proper cartesian axes w/ 'xy'
%   % 22-10-08: flip second axis for corrct orientation wrt scanner 
% extra arguments:
%       phase or p - phase map rather than absolute
%       mask or m - next argument is mask, multiplies the image
%       slice or sl - next argument is the slice number to display, at the
%       moment only works for 3-D
%       xy => use Cartesian axes
%       sc or scanner - flip second dim for scanner
%       gr or grid - use red gridlines
%       inv - inverts color map
%       nozero - computes range ignoring zero
%       nonan - removes NaN and Inf
%       rot - display rotated by next argument (degrees)
%***********************************************************

function h = imsjm(data,varargin)

data = squeeze(data);
cmap = 'jet';
phase_image = false;
slice = 1; %so far only for 3d
axis_definition = 'ij';
grid_on = false;
invert_map = false;
nozero=false;
nonan = false;
rot_disp=false;

for n=1:nargin-1
    % if a mask is sp[ecified
    if (strcmp(varargin{n},'mask')||strcmp(varargin{n},'m'))
        mask = varargin{n+1};
        data = data .* mask;
        continue
    end
    
    % If any 2 element vector is specified - assume this is the window
    if (numel(varargin{n})==2)&&(~ischar(varargin{n}))
        win = varargin{n};
        continue;
    end
        
    % gray colormap
    if (strcmp(varargin{n},'gray'))
        cmap='gray';
    end
    
    % phase image
    if (strcmp(varargin{n},'phase')||strcmp(varargin{n},'p'))
        phase_image = true;
    end
    
    % select slice
    if (strcmp(varargin{n},'slice')||strcmp(varargin{n},'sl'))
        slice = varargin{n+1};
    end
    
    % grid
    if (strcmp(varargin{n},'gr')||strcmp(varargin{n},'grid'))
        grid_on = true;
    end
    
    % invert color map
    if (strcmp(varargin{n},'inv'))
        invert_map = true;
    end
    
    %  set axis orientation
    if (strcmp(varargin{n},'xy'))
        axis_definition = 'xy';
    end
    %  set axis orientation
    if (strcmp(varargin{n},'scan'))||(strcmp(varargin{n},'scanner'))...
            ||(strcmp(varargin{n},'sc'))
       	data=flipdim(data,2); % flip second axis
%         disp('Flipped second axis for scanner orientation');
    end
    
    % ignore zeros
    if (strcmp(varargin{n},'nozero'))
        nozero = true;
    end
    
    % remove Nan and inf
        
    % ignore zeros
    if (strcmp(varargin{n},'nonan'))
        nonan = true;
    end
    
    % rot90
    if (strcmp(varargin{n},'rot'))
        rot_disp = true;
        rot = varargin{n+1};
    end
end


% squash data to 1 slice if too big    
if ndims(data)>2
    dims = ndims(data);
    switch dims
        case 3
            data = data(:,:,slice);
        case 4
            data = data(:,:,slice,1);
        case 5
            data = data(:,:,slice,1,1);
        case 6
            data = data(:,:,slice,1,1,1);
    end
end
data = squeeze(data);

% rotate
if rot_disp
    data = imrotate(data,rot);
end

if ~isreal(data)
    if ~phase_image
        data = abs(data);
    else
        data = angle(data);
    end
end

if nonan
    data(isnan(data)|isinf(data))=0;
end

% invert colormap
cmap = colormap(cmap);
if invert_map
    cmap = flipdim(cmap,1);
end

if exist('win','var')
   h = imagesc(data,win);colormap(cmap)
else
    if nozero
        % exclude zeros from window calculation
        mask = (data ~= 0);
        w1 = min(abs(data(mask)));
        w2 = max(abs(data(mask)));
        if w2==w1
            w2=w2*1.01;
        end
        imagesc(data,[w1 w2]);colormap(cmap)
    else
        h=imagesc(data);colormap(cmap)
    end
end


% axes
switch axis_definition
    case 'ij'
        axis ij
    case 'xy'
        axis xy
end

% grid
if grid_on
    grid on
    set(gca,'xcolor',[1 0 0])
    set(gca,'ycolor',[1 0 0])
end
   
if nargout==0
    clear h
end
end