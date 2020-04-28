function [self,exins] = get_inputs_real_new(self,exins)
% convert data format of collected raw binnary to LSQML format
% Zhen Chen @ Cornell University, July, 2018

% self: the inputs for LSQML_codes
% exins: inputs constaining all experimental related parameters.

import core.*
    
if ~isfield(exins,'fname') || isempty(exins.fname)
    [fname,dirname] = uigetfile({'*.mat;*.raw','Dataset'},'Select *.mat or *.raw dataset file','./');
    fname=fullfile(dirname,fname);
    exins.fname=fname;
end

%%
if strcmpi(exins.fname(end-2:end),'raw')
    cbed=rawread(exins.fname,1);
    cbed(cbed<20)=0;
elseif strcmpi(exins.fname(end-2:end),'mat')
    d=load(exins.fname);
    if exins.dataclean == 1
        exins.rbf=d.rbf;
    end
    if isfield(d,'cbed')
        cbed=d.cbed;
    elseif isfield(d,'dp')
        cbed=d.dp;
        clear dp;
    elseif isfield(d,'CBEDsimu_dose')
        cbed=d.CBEDsimu_dose;
        clear CBEDsimu_dose;
    else
        error('check diffraction variable');
    end    
else
    error('dataset format not supported yet.\n');
end

%% check data fullfillment
if ~isfield(exins,'alpha0')
    exins.alpha0=input('Input convergence angle in mrad: \n');
end
if ~isfield(exins,'voltage')
    exins.voltage=input('Input electron energy in keV: \n');
end
if ~isfield(exins,'rot_ang')
    exins.rot_ang=0;
    disp('WARNING: Assume no relative rotation between diffraction and scan');
end

if ~isfield(exins,'scanStepSize_y')
    exins.scanStepSize_y=input('Input scan step size in 1st dimension in Angstrom: \n');
end

if ~isfield(exins,'scanStepSize_x')
    exins.scanStepSize_x=input('Input scan step size in 2nd dimension in Angstrom: \n');
end
%%
%     load(exins.fname); % load 
if isfield(exins,'Np_p')
    Np_p=exins.Np_p;
else
    Np_p=128;
    disp('Use 128 as the default diffraction size');
end

if isfield(exins,'crop_idx') && length(exins.crop_idx) == 4
    cbed=cbed(:,:,exins.crop_idx(1):exins.crop_idx(2),exins.crop_idx(3):exins.crop_idx(4));
end

if isfield(exins,'Nscan_inter') 
    cbed=cbed(:,:,1:exins.Nscan_inter:end,1:exins.Nscan_inter:end);
    exins.scanStepSize_x = exins.scanStepSize_x * exins.Nscan_inter;
    exins.scanStepSize_y = exins.scanStepSize_y * exins.Nscan_inter;
end

if strcmpi(exins.fname(end-3:end),'.raw') || exins.dataclean ~= 1 
    [sy,sx]=shiftmeasure(cbed,0.8);
    cbed=shift_pix(cbed,[sx,sy]);
    pacbed=mean(mean(cbed,3),4);
    [rdat]=rscan(pacbed,'dispflag',0);
    Icenter=mean(rdat(1:5))*0.5;
    rbf=find(rdat<=Icenter);
    rbf=rbf(1);
    exins.rbf=rbf;  
end

if exins.dataclean ~= 1
    fnameout=exins.fname(1:end-4);
    save(strcat(fnameout,'_mat4pty_clean.mat'),'cbed','pacbed','rbf','sy','sx','-v7.3');
end

[ndpy,ndpx,npy,npx]=size(cbed);
if ndpy < Np_p(1) % pad zeros
    diffraction=padarray(cbed,[(Np_p(1)-ndpy)/2,(Np_p(2)-ndpx)/2,0,0],0,'both');
else
    diffraction=crop_pad(cbed,Np_p);
end
clear cbed;    
% change to electrons
if ~isfield(exins,'ADU')
    exins.ADU=1;
end   

%% find empty diffraction and delete it
if isfield(exins,'rmall0') && exins.rmall0 
    img_s=squeeze(sum(sum(diffraction,1),2));
    exins.rm0_idx=find(img_s == 0);
end

%%
diffraction = diffraction / exins.ADU; 

diffraction=reshape(diffraction,Np_p(1),Np_p(2),[]);
Itot=mean(squeeze(sum(sum(diffraction,1),2)));
    
% low pass filter
if isfield(exins,'cutoff')
    exins.cutoff = min([exins.cutoff,fix(Np_p(1)/2)]); % modified 9/12/2019
    cc=circle3(exins.cutoff,Np_p(1),Np_p(2));
    for ii=1:size(diffraction,3)
       diffraction(:,:,ii)= diffraction(:,:,ii).*cc;
    end
end
scanStepSize_x=exins.scanStepSize_x;
scanStepSize_y=exins.scanStepSize_y;

[~,lambda]=electronwavelength(exins.voltage);
% if ~isfield(exins,'dk')
%     exins.dk=exins.alpha0/1e3/exins.rbf/lambda;
% end
 
if isfield(exins,'alpha0') && isfield(exins,'rbf')
    exins.dk=exins.alpha0/1e3/exins.rbf/lambda;
end

[xx,yy]=scan_position_rot(npx,npy,scanStepSize_x,scanStepSize_y,exins.rot_ang);% in Angstrom
    
% new df
if isfield(exins,'df')
    df=exins.df;
else
    df=0;
    exins.df=df;
end 

% new Cs
if isfield(exins,'cs')
    cs=exins.cs;
else
    cs=0;
    exins.cs=cs;
end
% change pixel
% if isfield(exins,'dx')
%     dx=exins.dx; % pixel size in final object, in A
%     disp('use dx and ignore dk');
%     exins.dk=1/dx/Np_p(1);
% else
    dx=1/Np_p(1)/exins.dk; %% assume Np_p(1)=Np_p(2)
    exins.dx=dx;
%end
xx=xx/dx;
yy=yy/dx;

%%
probe_positions=[xx(:),yy(:)];
probe_positions=single(probe_positions);
probe=generateProbeFunction(dx,Np_p(1),0,0,df,cs,1,exins.voltage,exins.alpha0,0);
probe=probe/sqrt(sum(sum(abs(probe.^2))))*sqrt(Itot)/sqrt(Np_p(1)*Np_p(2));
probe=single(probe);
if ~isfield(exins,'extra')
    exins.extra=0.2;
    disp('Use 20% extra points for object')
end
Np_o = get_object_extent(Np_p, probe_positions, exins.extra);

%%
if isfield(exins,'rm0_idx') && ~ isempty(exins.rm0_idx)
    probe_positions(exins.rm0_idx,:)=[];
end
Npos=length(probe_positions);
% put to self
self.Np_p = Np_p; %  resolution of the difr pattern 
self.Np_o = Np_o; %  resolution of the difr pattern 

self.Npos = Npos; 
self.probe{1} = probe;
self.probe_positions_0 = probe_positions; 
self.probe_positions = [];    

self.reconstruct_ind{1} = 1:self.Npos;

self.diffraction = diffraction;
clear diffraction;
self.object{1} =   single( (rand(self.Np_o))  ...
    .* exp(2i*pi*rand(self.Np_o) * 1e-1 ) ) ;  % not so much random 
self.object_orig = self.object{1};
self.probe_orig = probe;

%% save useful params
exins.Nscan=[npy,npx];
end
%%
function [c]=circle3(r1A,Nx,Ny)%,cx,cy
% generate circle disk, radius r1A, within a matrix, default size 2*r1A+1 

if nargin < 3
    Ny=Nx;
end

if nargin <2
    Nx=2*r1A+1;
end

kx = (-fix(Nx/2):Nx-fix(Nx/2)-1);
ky = (-fix(Ny/2):Ny-fix(Ny/2)-1);
[kyy,kxx] = meshgrid(ky,kx);
k2 = (kxx).^2+(kyy).^2;

c=k2<=r1A^2;
end
%%
function [cbed]=rawread(fname,irot180)
% read raw dataset from EMPAD and crop useful part
% irot180 =1 
%       do rotation 180 degrees 

if nargin < 2
    irot180 = 1;
end
if nargin < 1
    [fname,fdir]=uigetfile('*.raw','Open a *.raw file');
    fname=fullfile(fdir,fname);
end

x_cbed = 128;
y_cbed = 130;  % default size of PAD diffraction

fid=fopen(fname,'r');
fin = fread(fid,'*float32');
fclose(fid);
ntot=length(fin)/x_cbed/y_cbed;

%if ~exist(npx) 
	disp('assume # of scan pixels along x/y the same. ');
	npx=round(sqrt(ntot));
	if mod(ntot,npx) ~= 0
		disp('scan along x/y not the same.');
		npx=input('Input scan pixel along horizontal:\n');
	end
%end	

npy=ntot/npx;
if mod(npy,1)~= 0

%if length(fin) ~= npy * npx *x_cbed * y_cbed
    error('Data size inconsistent, either dataset damage or input wrong scanning pixels, please check!\n');
end

fun= reshape(fin,x_cbed,y_cbed,npy,npx); % reshape for data
clear fin;
cbed = single(fun(3:126,3:126,:,:)); % crop edges
clear fun;
    
if irot180 == 1
    cbedrot=rot90(cbed,2); % rotate CBED 180d
    clear cbed;
    cbedrot=permute(cbedrot,[3,4,1,2]); 
    cbedrot=rot90(cbedrot,2);   % rotate scanning 180d
      
    cbed=permute(cbedrot,[3,4,2,1]); % transpose scan x/y, as there is a storage difference. scanning along  2nd dimension of matrix, horizontal if plotting by imagesc
    clear cbedrot;
else
    cbed=permute(cbed,[3,4,2,1]); 
end
end

%%
function [cx,cy]=shiftmeasure(cbed,thresh)
% detect and shift the center of diffraction via average diffraction
if nargin < 2
    thresh=0.8;
end

pacbed=mean(mean(cbed,3),4);

temp=sort(pacbed(:),'descend');
Imax=mean(temp(1:100))*thresh;
pacbed=pacbed/Imax;
pacbed(pacbed>=1)=1;
pacbed(pacbed<0.1)=0;

[cy,cx]=centroid3(pacbed,1);
cy=round(cy);
cx=round(cx);
end
%%
function [cbed]=shift_pix(cbed,sft_yx)
% shift diffraction by intergers
sft_yx=round(sft_yx);
if sft_yx(1) < 0
    cbed=cbed(1:end+2*sft_yx(1),:,:,:);
elseif sft_yx(1) > 0
    cbed=cbed(2*sft_yx(1)+1:end,:,:,:);
end

if sft_yx(2) < 0
    cbed=cbed(:,1:end+2*sft_yx(2),:,:);
elseif sft_yx(2) > 0
    cbed=cbed(:,2*sft_yx(2)+1:end,:,:);
end

% pad 0 to 124x124 for each diffraction

ny=124-size(cbed,1);
nx=124-size(cbed,2);

if ny > 0
    cbed=padarray(cbed,[ny/2,0,0,0],0,'both');
end
if nx > 0
    cbed=padarray(cbed,[0,nx/2,0,0],0,'both');
end

end
%%
function [rdat,xcoord,ycoord] = rscan(M0,varargin)

% RDAT = RSCAN(M0,VARARGIN)
% Get radial scan of a matrix using the following procedure:
% [1] Get coordinates of a circle around an origin.
% [2] Average values of points where the circle passes through.
% [3] Change radius of the circle and repeat [1] until rprofile is obtained.
%
% For DEMO, run
% >> rscan_qavg();
% or
% >> rscan_qavg('demo','dispflag',0);
% >> plot(ans);
% or
% >> rdat = rscan_qavg();
% >> plot(rdat);
% or
% >> a = peaks(300);
% >> rscan_qavg(a);
% >> rscan_qavg(a,'rlim',50,'xavg',100);
% >> rscan_qavg(a,'rlim',25,'xavg',100,'dispflag',1,'dispflagc',1);
% >> rscan_qavg(a,'rlim',25,'xavg',100,'dispflag',1,'dispflagc',1, ...
%    'squeezx',0.7,'rot',pi/4);
%
% Draw Circle: 
% [ref] http://www.mathworks.com/matlabcentral/fileexchange/
%       loadFile.do?objectId=2876&objectType=file

if nargin < 1 
    disp('This is a Demo');
    M0 = 'demo';
end

if strmatch('demo',lower(M0),'exact')
    %M0 = peaks(200);
    [xx,yy] = meshgrid(linspace(-3,3,201));
    [phi,rho] = cart2pol(xx,yy);
    M0 = besselj(1,rho);
    clear xx yy phi rho;
end

xavg = fix(size(M0,2))/2+1;
yavg = fix(size(M0,1)/2)+1;
Rlim = floor(min(size(M0))/2)-1;
dispFlag = 1;
dispFlagC = 0;
rot = 0; %% radian
squeezx = 1; %% 0.80;
squeezy = 1;
rstep = 1;

if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); end
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'xavg', xavg = varargin{ni+1};
            case 'yavg', yavg = varargin{ni+1};
            case 'squeezx', squeezx = varargin{ni+1};
            case 'squeezy', squeezy = varargin{ni+1};
            case 'rot', rot = varargin{ni+1};
            case 'rlim', Rlim = varargin{ni+1};
            case 'dispflag', dispFlag = varargin{ni+1};
            case 'dispflagc', dispFlagC = varargin{ni+1};
            case 'rstep', rstep = varargin{ni+1};
        end
    end
end


yxz = size(M0);
Rbnd = floor(min(yxz)/2)-1;
if Rlim > Rbnd, Rlim = Rbnd; end
    
for nRho = 1:rstep:floor(Rlim),
NOP = round(2*pi*nRho);
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*round(nRho);
[X,Y] = pol2cart(THETA,RHO);
X = squeezx*X;
Y = squeezy*Y;
[THETA,RHO] = cart2pol(X,Y);
[X,Y] = pol2cart(THETA+rot,RHO);
X = X + xavg;
Y = Y + yavg;

if dispFlag,
    h1 = figure(100);clf;box on;
    set(h1,'position',[10 500 400 300]);
    set(h1,'units','pixels');
    set(gca,'units','pixels');
    imagesc(M0);axis image;hold on;
    % H = plot(X,Y,'c-');
    line([xavg xavg],[1 yxz(1)],'color',[1 1 0],'linewidth',1);
    line([1 yxz(2)],[yavg yavg],'color',[1 1 0],'linewidth',1);
end

%%%
X = round(X);
Y = round(Y);
%%%

clear dat uxy pxy mxy mx nx my ny;
dat = [X;Y];
uxy = diff(dat,1,2);
uxy = [[1;1],uxy];
pxy = union(find(uxy(1,:)~=0),find(uxy(2,:)~=0));
dat = dat(:,pxy);
rdat(nRho) = sum(M0(yxz(1)*(dat(1,:)-1) + dat(2,:)))/length(dat);

if dispFlag,
    H = plot(dat(1,:),dat(2,:),'y-');
    if dispFlagC,
        for nrn = 1:length(dat)
            H = plot(dat(1,nrn),dat(2,nrn),'m.','MarkerSize',12);
            drawnow;
            delete(H);
        end
    end
    drawnow;
end
end

% rdat = rdat/max(rdat);
xcoord = dat(1,:);
ycoord = dat(2,:);

if dispFlag,
h2 = figure(101); clf; hold on;
set(h2,'position',[10 100 400 300]);
plot(rdat,'rx','linewidth',2);axis tight;
% plot(M0(round(xavg)+1,round(yavg)+1:end),'b','linewidth',2);
end
end
%%
function [probe,mask,A] = generateProbeFunction(dx,N, px,py,df ,cs, N_z, Voltage, alpha_max, showFigure )
%Generate probe function
%   dx: pixel size in angstrom
%   N: number of pixels in each side
%   df: defocus in angstrom
%   Nz:
%   Voltage: beam voltage in keV
%   alpha_max: convergence angle in mrad


%Parameter
C3 = cs; %angstrom
C5 = 0;
C7 = 0;
lambda = 12.398/sqrt((2*511.0+Voltage).*Voltage); %angstrom

amax = alpha_max*1e-3; %in rad
amin = 0;

klimitmax = amax/lambda;
klimitmin = amin/lambda;
dk = 1/(dx*N);

if N_z==1
    dff = df;
else
    dff = linspace(-floor(N_z/2),ceil(N_z/2)-1,N_z)*df; %angstrom
end

kx = linspace(-floor(N/2),ceil(N/2)-1,N);
[kX,kY] = meshgrid(kx,kx);

kX = kX.*dk;
kY = kY.*dk;
kR = sqrt(kX.^2+kY.^2);

chi = zeros(N,N,N_z);
probe = zeros(N,N,N_z);

mask = single(kR<=klimitmax).*single(kR>=klimitmin);
%figure
%imshow(mask)
for i= 1:N_z
    chi(:,:,i) = -pi*lambda*kR.^2*dff(i) + pi/2*C3*lambda^3*kR.^4+pi/3*C5*lambda^5*kR.^6+pi/4*C7*lambda^7*kR.^8;
    
    phase = exp(-1i.*chi(:,:,i)).*exp(-2*pi*1i*px.*kX).*exp(-2*pi*1i*py.*kY);
    A = mask.*phase;
    probe(:,:,i) = mask.*phase;
    
    probe(:,:,i) = fftshift(ifft2(ifftshift(probe(:,:,i))));
    
    probe(:,:,i) = probe(:,:,i)/sum(sum(abs(probe(:,:,i))));
    if showFigure
        figure
        imagesc(abs(probe))
        axis image
    end
end

if showFigure
    figure
    imagesc(abs(mask))
    axis image

end
%probe = probe/(L/N)^2; %normalize the probe. Because in reciprocal space
%the zero frequency amplitute is 1, so, the summation is also one.
end
%%
function [wav,wavA] = electronwavelength(kev)
            %wavelength in a0
            a0 = 0.52917720859;
            wav = 12.3986./sqrt((2*511.0+kev).*kev)/a0;
            wavA = wav*a0;
end 
%%
function [ppX_rot,ppY_rot]=scan_position_rot(N_scan_x,N_scan_y,scanStepSize_x,scanStepSize_y,rot_ang)
% scan positions after rotation
    ppx = linspace(-floor(N_scan_x/2),ceil(N_scan_x/2)-1,N_scan_x)*scanStepSize_x;
    ppy = linspace(-floor(N_scan_y/2),ceil(N_scan_y/2)-1,N_scan_y)*scanStepSize_y;
    [ppX,ppY] = meshgrid(ppx,ppy);

    ppY_rot = ppX*(-sind(rot_ang)) + ppY*cosd(rot_ang);
    ppX_rot = ppX*cosd(rot_ang) + ppY*sind(rot_ang);

end

%%

% end