clear
close all
%% load images
%Ni_canny_image = ground truth (GT)
%Nd_fuzzy_image = extracted edge (EE)

Nd_fuzzy_image = double(imread('ACO_sphere.bmp'))/255;
Ni_canny_image = double(imread('canny_sfphere.bmp'))/255;
[m,n] = size(Nd_fuzzy_image);
maxD = sqrt((m-1)^2+(n-1)^2);

%% show images
figure
imshow(Nd_fuzzy_image)

figure
imshow(Ni_canny_image)

%% build edge pixels index matrix


%build non-edge pixels index matrix
%find pozitions in which GT != EE
for x = 1:m
    for y = 1:n
        my_diff(x,y) = sqrt((Nd_fuzzy_image(x,y)-Ni_canny_image(x,y))^2);
        
        if (my_diff(x,y) ~=0 )
            distances = ones(1,1);
            wi=1;
            indexes = 130*ones(wi,2);
            for i=1:m
                for j=1:n
                    if Nd_fuzzy_image(i,j) == Ni_canny_image(x,y)
                        indexes(wi,1) = i;
                        indexes(wi,2) = j;
                        wi = wi+1;
                    end
                end
            end
            for kk = 1:wi-1
                distances(kk) = sqrt((x-indexes(kk,1))^2+(y-indexes(kk,2))^2);
            end
            sort_distances =sort(distances);
            minimum_distance = min(distances);
            my_diff(x,y) = minimum_distance;
        end
    end
end

%% show distance matrix
figure
imshow(my_diff)

FD = double(my_diff)/maxD; 

card_Ni = sum(sum(Ni_canny_image));
card_Nd = sum(sum(Nd_fuzzy_image));




%% matlab utilities
FTP = min(Ni_canny_image,Nd_fuzzy_image); %fuzzy intersection
card_FTP = sum(sum(FTP));
% bounded difference between Ni and Nd
FFP = max(0,Ni_canny_image-Nd_fuzzy_image); 
card_FFP = sum(sum(FFP));
% bounded difference between Nd and Ni
FFN = max(0,Nd_fuzzy_image-Ni_canny_image); 
card_FFN = sum(sum(FFN));
for i = 1:m
    for j = 1:n
        matrix_SUM(i,j) = 1/(1+FD(i,j));
    end
end
SUM = sum(sum(matrix_SUM));

% [m,n] = size(Nd_fuzzy_image);

FI = 1/(m*n*max(card_Ni,card_Nd)) * (card_FTP*SUM-card_FFP-card_FFN);%fuzzy_index


