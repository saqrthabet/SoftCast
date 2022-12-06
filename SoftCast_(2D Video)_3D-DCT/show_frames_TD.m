function show_frames_TD(frame,gop,o)
%SAQR
figure(o);
for i = 1 : gop
        subplot(gop/2,2,i);
        %imshow(mat2gray(frame(:, :, i)));
        imshow(uint8(frame(:, :, i)));
        title(['Frame #', num2str(i)])
        xlabel(size(frame(:, :, i),2)); ylabel(size(frame(:, :, i),1));
end

end

