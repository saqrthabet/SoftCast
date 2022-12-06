% Read video data from file, each time for ONE GOP
function VIn = ReadVideoDataFromFile(videoFileName, videoWidth, videoHeight, GopSize, tFrame, videoOrder)

sz = size(videoOrder);     % 1  2
SeqNum = sz(2);      %[ 1  5] No.of columns=2

VIn = zeros(videoHeight, videoWidth, GopSize, SeqNum); 

for s = 1 : SeqNum
   SeqName = sprintf(videoFileName,videoOrder(s));
   fp = fopen(SeqName, 'rb');    % r=read , b=Big-endian ordering
   for i = 1 : GopSize        
       fseek(fp, ((tFrame-1)*GopSize+(i-1))*videoWidth*videoHeight*3/2, 'bof');
       tmp = fread(fp, [videoWidth, videoHeight]);
       VIn(:,:,i,s)= uint8(tmp');
       %imshow(uint8(VIn(:,:,i,s))); % Just for Debug
   end
   fclose(fp);
end
end