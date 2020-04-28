function [O_phase,O_amp]=pty_ramprm_f(outputs, exins,dirbase)
% fast rotate back, crop edges and remove phase ramp
% Zhen Chen, Cornell University

removePhaseRamp=1;

%% load inputs/exins/outputs
fout='MLs_backgroundremove_final_';
if ~exist('dirbase','var')
    dirbase='./';
end
fout=fullfile(dirbase,fout);
%%
if isfield(exins,'scanStepSize_x')
    scanStepSize_x = exins.scanStepSize_x;
    scanStepSize_y = exins.scanStepSize_y;
else
    scanStepSize_y = input('Input scan step size along 1-st dimension of matrix in Angstrom:\n');
    scanStepSize_x = input('Input scan step size along 2-nd dimension of matrix in Angstrom:\n');
end

if isfield(exins,'Nscan')
    npy = exins.Nscan(1);
    npx =exins.Nscan(2);
else
    npy = input('Input scan pixel number along 1-st dimension in Angstrom:\n');
    npx = input('Input scan pixel number along 2-nd dimension in Angstrom:\n');
end

if isfield(exins,'dx')
    dx_x = exins.dx;
else
    dx_x = input('Input pixel size of final reconstructed object in Angstrom:\n');
end

dx_y = dx_x ;
if isfield(exins,'rot_ang')
    rot_ang = exins.rot_ang;
else
    rot_ang = input('Input rotate angle used in the reconstruction:\n');
end

%%
if rot_ang ~= 0
    object_rot = outputs.object{1};
    object_rot=imrotate(object_rot,-rot_ang,'bilinear');
else
    object_rot=outputs.object{1};
end
cx=fix(size(object_rot,2)/2)+1;
cy=fix(size(object_rot,1)/2)+1;

lb_x=cx-round(fix(npx/2)*scanStepSize_x/dx_x);
lb_y=cy-round(fix(npy/2)*scanStepSize_y/dx_y);
ub_x=cx+round((ceil(npx/2)-1)*scanStepSize_x/dx_x);
ub_y=cy+round((ceil(npy/2)-1)*scanStepSize_y/dx_y);

object_crop=double(object_rot(lb_y:ub_y,lb_x:ub_x));

O_amp=abs(object_crop);

%%
if removePhaseRamp == 1
    [ O_phase ] = remove_phase_ramp( object_crop,dx_x,dx_y);
else
    O_phase=angle(object_crop);
end
%%
O_phase_C=mat2gray_enhC(O_phase,0.01,0.999);
imwrite(O_phase_C,strcat(fout,'crop_phase.png'),'png');
imwrite(mat2gray(O_amp),strcat(fout,'crop_amp.png'),'png');

%%
save(strcat(fout,'crop_data.mat'),'dx_x','dx_y','O_amp','O_phase','O_phase_enhC','removePhaseRamp');

%% related functions
function [ output ] = remove_phase_ramp( input,dx_x,dx_y,varargin )
% remove a linear phase ramp

output = zeros(size(input));
Ny = size(input,1);
Nx = size(input,2);
N_recon = size(input,3);
y = linspace(-floor(Ny/2),ceil(Ny/2)-1,Ny);
x = linspace(-floor(Nx/2),ceil(Nx/2)-1,Nx);
[X,Y] = meshgrid(x,y);
X = X.*dx_x;
Y = Y.*dx_y;

if isempty(varargin)
    mask_thresholds = min(min(abs(input)));
elseif length(varargin)==1
    mask_thresholds = ones(N_recon,1)*varargin{1};
elseif length(varargin)==2
    mask_thresholds = ones(N_recon,1)*varargin{1};
    if varargin{2}>0
        figure
        mag_image = abs(input(:,:,1));
        imagesc(mag_image>mask_thresholds(1))
        title('masked magnitude')
    end
end
for i = 1:N_recon
    phase_image = angle(input(:,:,i));
    mag_image = abs(input(:,:,i));
    mask = mag_image>mask_thresholds(i);
    
    %fit
    fo = fit([X(mask), Y(mask)],phase_image(mask),'poly11');
    %background = X*fo.p10+Y*fo.p01+fo.p00;
    background = X*fo.p10+Y*fo.p01;
    output(:,:,i) = phase_image - background;
end
end
%%
    function [imgout,Imin,Imax]=mat2gray_enhC(img,plow,phigh)
    % mat 2 gray image with enhanced contrast by assuming [plow,phigh]
    % percentage intensity outside truncated
    if nargin < 3
        phigh=1;
    end
    if nargin < 2
        plow=0;
    end
    if nargin < 1
        error('input an image matrix')
    end

    temp=sort(img(:));
    temp(isnan(temp))=[];
    Imin=temp(round(plow*length(temp))+1);
    Imax=temp(round(phigh*length(temp)));
    imgout=mat2gray(img,[Imin,Imax]);
    end

end