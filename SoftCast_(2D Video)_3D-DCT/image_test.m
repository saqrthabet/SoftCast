function out_xx=image_test( xx,snr_vec )
out_xx=[];
im=[];

  for i=1:length(snr_vec)
    im=[im,xx(:,:,i)];
      if mod(i,4) == 0                % when construction of one strip of series of chunks is complete,as number of culumns=width=number of chunks laying horizentally
            out_xx = [out_xx;im];     % it jumps to next line to finish constructing the frame, xx is accumulating the frame 
            im = [];                  % tt is cleared to start over 
      end
  end

end

