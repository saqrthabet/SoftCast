%% Saqr Thabet

%% 3D SoftCast on 2D videos 
% sending multiple Frames 

clc
clear all
close all

VideoWidth = 352;%352 512
VideoHeight = 288;%288 512
inpara.skip=1;
VideoFrame = 4;%48;
GopSize = 4;   %8;

bw = VideoWidth/8;%4
bh = VideoHeight/8;%3
blockNum = bw*bh*GopSize;
picsize=VideoWidth *VideoHeight;

dimension=8;
chunk_size=dimension*dimension;

tFrame = VideoFrame/ GopSize;
P1 = 1;
TotalPower = P1 *8*8; 
MAX_TIMES = 10;

% read video samples
VideoFile = '/Users/saqr/Downloads/foreman_cif.yuv';

% loading Data from the 2D video 
%fid = fopen(VideoFile,'r');
%  for frame =1:GopSize
%      fseek(fid, (1+frame-2)*inpara.skip*picsize*1.5, -1);  
%      I0 = fread (fid, VideoWidth * VideoHeight, 'uint8');
%      I0 = reshape(I0, VideoWidth, VideoHeight);
%      I(:,:,frame)=double(I0');
%      imshow(uint8(I(:,:,frame)));
%  end
  
  [I,~,~]=yuvRead(VideoFile,VideoWidth,VideoHeight,GopSize);
  for frame =1:GopSize
       imshow(uint8(I(:,:,frame)));
  end
  

snr_vec=[0:20];
xx=zeros(size(I,1),size(I,2),GopSize,length(snr_vec));
xx1=zeros(size(I,1),size(I,2),GopSize,length(snr_vec));
xx2=zeros(size(I,1),size(I,2),GopSize,length(snr_vec));
xx1_mold=zeros(size(xx1,1),size(xx1,2),GopSize);
PSNR=zeros(1,length(snr_vec));

%% Encoder
%horizental stream of(all)dct2 frequency bins
%CHUNKS are CLUSTERING
%         dctM = dct2(I-128);%mirt_dctn
%         %show_frames_FD(dctM,4,2);
% [ChunkStruct,ChunkDctCoeffStdSum] = divideDctCoeffIntoChunks...
%     (dctM, bw, bh, VideoHeight, VideoWidth, GopSize, TotalPower);

x = [];d=0;
 I_dct=mirt_dctn(I-128);%dct2(I-128);
 
 for kk=1:GopSize
  for ii = 1:bh     %1:32  repetition upon vertical number/pick row after another second
    for jj = 1:bw %1:32  repetition upon horizental number/pick colum by colum first
        d=d+1;
        currentBlock = I_dct((ii-1)*8+1:ii*8,(jj-1)*8+1:jj*8,kk);  %chunk dimension of 8X8 everystep, picking chunks as keyboard calculator 7->8->9->4 due to 
        %mean_(d)=mean(mean(currentBlock));currentBlock_zero_mean=currentBlock-mean(mean(currentBlock));
        x = [x ,reshape(currentBlock,chunk_size,1)];   %DCT of each block then convert each block into 64X1 colum vector  FD eventually, 
    
    end                                                                       %  64X1024 every chunk is COLUMN
  end
 end
 
 %% Power Scaling-Optimization problem <variance was General and was not assigned to each chunk>
%
P1=1;
%P1 = sum(sum(x.*x))/(64*blockNum);   %sum(x.*x)=Autocorrelation(lag=0)or power of distribution(E(x^2))___1/(64*blockNum=total number of pixels)=normalization
P = P1 * chunk_size*4; 
mean(x,2);  % mean of each row, mean(x):mean of each column
lamda = mean((x.*x)');    %variance of input signal  obsolete definition-while mean=0
lamda =  lamda';
g = sqrt(P/sum(sqrt(lamda)))./sqrt(sqrt(lamda));     %optimal Scaling factor per chunk that balance power

%% Hadamard 
H = hadamard(chunk_size)/8;         %why hadamard value is -/+(1/8) not -/+ 1

%% Multiplication
% C1 = diag(g);
% ytx = C1 *x;
% sum(sum(ytx.*ytx))/size(ytx,1)/size(ytx,2)
C = H*diag(g);                                                                         % encoding matrix. diag(g) extends row vector (g)into diangonal matrix
ytx = C*x; 

%% Channel 
 % AWGN 

for snr_ind = 1:length(snr_vec)
    snr=snr_vec(snr_ind);                                                     % snr_vect is an array of values

    noise_pow = P1 * 10^(-snr/10);
    noisy = sqrt(noise_pow)*randn(chunk_size,blockNum);
    ynoisy = ytx + noisy;                                                        % additive noise
    %     sum(sum(noisy.*noisy))/size(noisy,1)/size(noisy,2)
    yrx = ynoisy; 
    
    %% ------------------- Decoder --------------------
    %----------------------SoftCast--------------
    sigma = mean((noisy.*noisy)');         % variance of noise for every ROW
    sigma = sigma';
    noise_pow = diag(noise_pow*ones(chunk_size,1));                % power of noise make row vector 64x1 becomes diagonal matrix of 64x64
    z = diag(lamda)*C'*inv(C*diag(lamda)*C'+ noise_pow)*yrx;  % Linear Least Square Estimator (LLSE)
    %z = diag(lamda)*C'\(C*diag(lamda)*C'+ noise_pow)*yrx;  % Linear Least Square Estimator (LLSE)
    % z = inv(C)*yrx;                                                            %At high SNR when noise_pow =0
     %% iDCT
    xx0 = [];xx1=[];
    tt = [];kk=0;
    for ii = 1:blockNum                         % total number of chunks
        temp = z(:,ii);%+mean_(ii);                         % selecting ONe column of matrix(z) which comprises of 64 elements(rows)
        currentBlock = reshape(temp,dimension,dimension);       % rebuilding chunk 8x8
        tt = [tt currentBlock];     %inverse dct2 of chunk block matrix 8x8
        if mod(ii,bw) == 0                              % when construction of one strip of series of chunks is complete,as number of culumns=width=number of chunks laying horizentally
            xx0 = [xx0;tt];                             % it jumps to next line to finish constructing the frame, xx is accumulating the frame 
            tt = [];                                    % tt is cleared to start over 
        end
        if mod(ii,(blockNum/GopSize)) == 0                              % when construction of one strip of series of chunks is complete,as number of culumns=width=number of chunks laying horizentally
            kk=kk+1;
            xx1(:,:,kk) = xx0;                             % it jumps to next line to finish constructing the frame, xx is accumulating the frame 
            tt = [];xx0 = [];                                  % tt is cleared to start over 
        end
    end
    xx00=mirt_idctn(xx1)+128;%idct2(xx1)+128;%
    %% Before setting the Limits
    xx(:,:,:,snr_ind)=xx00;
    
    %% After setting the Limits
   for kk = 1 : size(xx,3)
    for ii = 1 : size(xx,1)
        for jj = 1 : size(xx,2)
            if xx(ii,jj,kk,snr_ind)>255
                xx2(ii,jj,kk,snr_ind) = 255;          % set a lower limit of 0
            elseif xx(ii,jj,kk,snr_ind)<0             % set a upper limit of +255
                xx2(ii,jj,kk,snr_ind) = 0;
            else
                xx2(ii,jj,kk,snr_ind)=round(xx(ii,jj,kk,snr_ind));
            end
        end
     end
   end
       [ssimval, ssimmap] = ssim(uint8(xx2(:,:,:,snr_ind)),uint8(I));
       fprintf('The SSIM value is %0.4f.\n',ssimval);
%        figure(4);
%        ax3=subplot(3,2,4);
%        imshow(ssimmap,[]);
%        title(sprintf('ssim Index Map - Mean ssim Value is %0.4f',ssimval));
       xx2_mold=zeros(size(xx2,1),size(xx2,2),size(xx2,3));
       xx2_mold= xx2(:,:,:,snr_ind);
       MSE=sqrt(mean((double(I(:))-(xx2_mold(:))).^2));
       psnr_programmed= [20*log10(255/MSE)];             %Peak Signal-to-Noise Ratio (PSNR)
       fprintf('The PSNR value is %0.4f.\n', psnr_programmed);    %The acceptable range of SNR is 30db to 50db.
       PSNR(snr_ind)=psnr_programmed;
       
end
% images test
out_xx=image_test( xx,snr_vec );
figure(3);imshow(out_xx,gray(256)),title('before setting limits SNR=[1,20]');

out_xx2=image_test( xx2,snr_vec );
figure(4); imshow(out_xx2,gray(256)),title('after setting limits SNR=[1,20]');   

figure(7);
plot(snr_vec,PSNR,'-*b');
xlim([0 20]),ylim([5 50])
title('PSNR at Rx')
legend('2D-SoftCast','Direct transmission')

