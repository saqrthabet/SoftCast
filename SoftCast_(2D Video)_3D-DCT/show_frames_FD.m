function show_frames_FD(frame,gop,o)
%SAQR
figure(o);
for i = 1 : gop
        subplot(gop/2,2,i);
        %imshow(mat2gray(frame(:, :, i)));
        %imshow(uint8(frame(:, :, i)));
        imshow(log(abs(frame(:, :, i))),[]);colormap(gca,jet(64)); h = colorbar; caxis([-4 10]) 
        title(['Frame #', num2str(i)])
        xlabel(size(frame(:, :, i),2)); ylabel(size(frame(:, :, i),1));
end

end

