%% Saqr Thabet

%% The modification:
% 1- block size to 32X32=1064 &total No.of blocks 64
% 2- every chunk in Old code is represented by ONE COLUMN while in theory it 
% has to be represented as a ONE ROW.
% 3- 128 can be subtracted in TD just before DCT and can be returned
% directly after IDCT
% 4- we also subtracted the mean from x
% 5- 2D-DCT two methods: Block based & Frame based

% old functions will be commented(%-)

clc
clear;
close all;
datestr(now)
%% 
RGB = imread('visual_samples/lena512color.tiff'); % supposed to be omited
I=rgb2gray(RGB);
I = double(I);                         % input frame/stream
[height, width] = size(I);     
%function  out=soft_cast(I,snr_vec)        % snr_vect affects yrx because we will have different number of noise =length(snr_vect)
%% DCT Block Dimension
dimension=8;
chunk_size=dimension*dimension;
bw = width/dimension;                         % it tells how many chunks Horizentally    
bh = height/dimension;                        % it tells how many chunks Vertically 
blockNum = bw*bh;               % total number of chunks in single frame GoP/4

snr_vec=[0:20];
xx=zeros(size(I,1),size(I,2),length(snr_vec));
xx1=zeros(size(I,1),size(I,2),length(snr_vec));
xx1_mold=zeros(size(xx1,1),size(xx1,2));
PSNR=zeros(1,length(snr_vec));


%% ----------------- Encoder------------
%horizental stream of(all)dct2 frequency bins
%CHUNKS are CLUSTERED

% Method 1 Chunk based 2D-DCT
% x = [];                                         % build a stream(row vector) of chunks 
%  for ii = 1:bh     %1:8  repetition upon vertical number/pick row after another second
%     for jj = 1:bw %1:8  repetition upon horizental number/pick colum by colum first
%        
%         currentBlock = I((ii-1)*dimension+1:ii*dimension,(jj-1)*dimension+1:jj*dimension);
%         currentBlock=dct2(currentBlock-128);
%         currentBlock_zero_mean=currentBlock-mean(mean(currentBlock));
%         %-x = [x; reshape(dct2(currentBlock-128),1,1024)];
%         x = [x, reshape(currentBlock_zero_mean,chunk_size,1)];% -mean(reshape(dct2(currentBlock-128),1,1024))]; 
%          
%     end                                                                       %  64X1024 every chunk is COLUMN
%   end                                                                                                  % while in theory each chunk has to be a ROW
  
% Method 2 Frame based 2D-DCT

x = [];
 I_dct=dct2(I-128);%
 for ii = 1:bh     %1:32  repetition upon vertical number/pick row after another second
    for jj = 1:bw %1:32  repetition upon horizental number/pick colum by colum first
        currentBlock = I_dct((ii-1)*8+1:ii*8,(jj-1)*8+1:jj*8);  %chunk dimension of 8X8 everystep, picking chunks as keyboard calculator 7->8->9->4 due to 
        %currentBlock_zero_mean=currentBlock-mean(mean(currentBlock));
        x = [x ,reshape(currentBlock,64,1)];   %DCT of each block then convert each block into 64X1 colum vector  FD eventually, 
    end                                                                       %  64X1024 every chunk is COLUMN
  end
%%
% HOW TO REMOVE THE MEAN FROM EACH CHUNK AND MAKE IT 'ZERO'
% by subtracting mean(x)of every row from all row elements, and send that mean(x)of every row in 'METADATA'
% 'METADATA' are just 'mean' and 'variance' of each chunk, sent via Huffman coding with BPSK & half-rate convolution code(FEC)
% its overhead impact is insignifficant 0.014bits/pixel

%% Power Scaling-Optimization problem <variance was General and was not assigned to each chunk>
%
P1=1;
%P1 = sum(sum(x.*x))/(64*blockNum);   %sum(x.*x)=Autocorrelation(lag=0)or power of distribution(E(x^2))___1/(64*blockNum=total number of pixels)=normalization
P = P1 * chunk_size; 
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
ytx = C*x;                                                                                 % encode input stream by Encoding matrix  'Slice' is output
% sum(sum(ytx.*ytx))/size(ytx,1)/size(ytx,2)

%% Channel 
 % AWGN 

for snr_ind = 1:length(snr_vec)
    snr=snr_vec(snr_ind);                                                     % snr_vect is an array of values

    noise_pow = P1 * 10^(-snr/10);
    noisy = sqrt(noise_pow)*randn(chunk_size,blockNum);
    ynoisy = ytx + noisy;                                                        % additive noise
    %     sum(sum(noisy.*noisy))/size(noisy,1)/size(noisy,2)
    yrx = ynoisy;                                                                 % ytx with additive noise
                        %At high SNR when noise_pow =0

    %% Decoder
    sigma = mean((noisy.*noisy)');         % variance of noise for every ROW
    sigma = sigma';
    noise_pow = diag(noise_pow*ones(chunk_size,1));                % power of noise make row vector 64x1 becomes diagonal matrix of 64x64
    z = diag(lamda)*C'*inv(C*diag(lamda)*C'+ noise_pow)*yrx;  % Linear Least Square Estimator (LLSE)

     %% iDCT
    xx0 = [];
    tt = [];
    for ii = 1:blockNum                         % total number of chunks
        temp = z(:,ii);                         % selecting ONe column of matrix(z) which comprises of 64 elements(rows)
        currentBlock = reshape(temp,dimension,dimension);       % rebuilding chunk 8x8
        tt = [tt currentBlock];  %idct2(currentBlock)+128    %inverse dct2 of chunk block matrix 8x8
        if mod(ii,bw) == 0                              % when construction of one strip of series of chunks is complete,as number of culumns=width=number of chunks laying horizentally
            xx0 = [xx0;tt];                             % it jumps to next line to finish constructing the frame, xx is accumulating the frame 
            tt = [];                                    % tt is cleared to start over 
        end
    end
    xx00=idct2(xx0)+128;
    %% Before setting the Limits
    xx(:,:,snr_ind)=xx00;
    
    %% After setting the Limits
    for ii = 1 : size(xx,1)
        for jj = 1 : size(xx,2)
            if xx(ii,jj,snr_ind)>255
                xx1(ii,jj,snr_ind) = 255;          % set a lower limit of 0
            elseif xx(ii,jj,snr_ind)<0             % set a upper limit of +255
                xx1(ii,jj,snr_ind) = 0;
            else
                xx1(ii,jj,snr_ind)=round(xx(ii,jj,snr_ind));
            end
        end
    end

       [ssimval, ssimmap] = ssim(uint8(xx1(:,:,snr_ind)),uint8(I));
       fprintf('The SSIM value is %0.4f.\n',ssimval);
%        figure(4);
%        ax3=subplot(3,2,4);
%        imshow(ssimmap,[]);
%        title(sprintf('ssim Index Map - Mean ssim Value is %0.4f',ssimval));
       xx1_mold=zeros(size(xx1,1),size(xx1,2));
       xx1_mold= xx1(:,:,snr_ind);
       MSE=sqrt(mean((double(I(:))-(xx1_mold(:))).^2));
       psnr_programmed= [20*log10(255/MSE)];             %Peak Signal-to-Noise Ratio (PSNR)
       fprintf('The PSNR value is %0.4f.\n', psnr_programmed);    %The acceptable range of SNR is 30db to 50db.
       PSNR(snr_ind)=psnr_programmed;
end

    figure(5);
    plot(snr_vec,PSNR,'-*b');
    xlim([0 20]),ylim([5 50])
    legend('BDCT 8x8') % it is important to write the Dimension of BDCT
    title('PSNR at Rx')
    
    % images test
    out_xx=image_test( xx ,snr_vec);
    figure(3);imshow(out_xx,gray(256)),title('before setting limits SNR=[1,20]');

    % images test
    out_xx1=image_test( xx1,snr_vec);
   figure(4); imshow(out_xx1,gray(256));title('after setting limits SNR=[1,20]');
