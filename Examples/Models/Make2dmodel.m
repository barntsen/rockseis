clear all
close all

Nx=201;
Nz=301;

Dx=5;
Dz=5;

%%% Constant velocity model
Vp=ones(Nx,1,Nz)*2000;
Vs=ones(Nx,1,Nz)*1500;
Rho=ones(Nx,1,Nz)*1000;

%%% Add reflectors
Vp(:,:,150:end) = 2500;
Vs(:,:,150:end) = 1800;
Rho(:,:,150:end) = 2000;

Header_models = RSShead(Vp);
Header_models.geometry.D(1)=Dx;
Header_models.geometry.D(3)=Dz;

WriteRSS('Vp2d.rss', Vp, Header_models);
WriteRSS('Vs2d.rss', Vs, Header_models);
WriteRSS('Rho2d.rss', Rho, Header_models);

%%% Make wavelet
Nt=8001;
Dt=5e-4;
wav2d = Ricker(15, 0.1, Nt, Dt).';

Head_wavelet = RSShead(wav2d);
Head_wavelet.geometry.D(1) = Dt;
Head_wavelet.type = 2;
Head_wavelet.Nheader = 4;
Head_wavelet.coordinates.s(1) = 0;
Head_wavelet.coordinates.s(2) = 0;
Head_wavelet.coordinates.g(1) = 0;
Head_wavelet.coordinates.g(2) = 0;

WriteRSS('Wav2d.rss',wav2d, Head_wavelet);




