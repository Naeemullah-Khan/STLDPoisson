function im_tile=image_tiling2(im)
[l,w,h]=size(im);
count=ceil((0.5*h)^0.5);
for i=1:count
    U(:,:,i)=zeros(l,2*count*w);
    for j=1:2*count
        
        if 2*count*(i-1)+j>h
            break
        end
        U(:,1+w*(j-1):w*j,i)=im(:,:,2*count*(i-1)+j)/(max(max(max(im)))+eps);
        umat(i,j)=2*count*(i-1)+j;
    end
end

im_tile=[];
for i=1:count
    im_tile=cat(1,im_tile,U(:,:,i));
end
% count
% i
% j
% umat