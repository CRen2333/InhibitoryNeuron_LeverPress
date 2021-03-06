function WarpedImage = WarpImage(ImageMatrix, N, tform)

    for kk = 1:size(ImageMatrix,2)
        temp_image = ImageMatrix(:,kk);
        temp_image = reshape(temp_image,[N N]);
        Rfixed = imref2d([N N]);
        temp_image = imwarp(temp_image,tform,'OutputView',Rfixed);
        WarpedImage(:,kk) = temp_image(:);
    end
    
end