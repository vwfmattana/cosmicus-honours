% Pierre Roux 22258116                             Vincent Mattana 21128707
% =========================================================================
% =========================================================================
% START of Image Enhancement
% =========================================================================
% =========================================================================
%--------------------Input Image-------------------------------------------
clear;
[filename, pathname] = ...
     uigetfile({'*.bmp;*.jpg;',...
 'Picture Files (*.bmp,*.jpg)';},'Select source image for processing');
originalImage=imread(strcat(strcat(pathname,'\'),filename) );
f=rgb2hsv(double((imread(strcat(strcat(pathname,'\'),filename) ))));

%--------------------Paths-------------------------------------------------
path='F:\Ground Truth\Generate\'; 

%--------------------Baseline Perimeters-----------------------------------
d=3; %----gap reach used in (SCALELINE finder) L.114
hm_index=0.04; %hourmarker visibility parameter

dimension=size(f);
ndim=dimension(1);
mdim=dimension(2);
% % figure, imshow(f(:,:,3)/255);figure(gcf);

%--------------------Image Features----------------------------------------
H=fspecial('prewitt');
%create horizontal prewitt mask

g_features=imfilter(f,H,'replicate');
g_features_scaled(:,:,3)=abs(g_features(:,:,3));
%binarize image using otsu's method(greythresh)
threshold=graythresh(g_features_scaled(:,:,3)/255);

%===+++   BINARIZE   +++===%
g_features_scaled_bw = im2bw(g_features_scaled(:,:,3)/255, threshold+log(threshold+1)*0.1);  

%applies vertical prewitt filter to f, leaving sharp-mask of vertical lines
g_features_sideways=imfilter(f,transpose(H),'replicate');
g_features_scaled_sideways(:,:,3)=abs(g_features_sideways(:,:,3));
g_features_scaled_bw_sideways = im2bw(g_features_scaled_sideways(:,:,3)/255, graythresh(g_features_scaled_sideways(:,:,3)/255));
imwrite(g_features_sideways(:,:,3),strcat(path,'0deg prewitt.bmp'),'bmp');
imwrite(g_features_scaled_sideways(:,:,3),strcat(path,'0deg abs.bmp'),'bmp');
imwrite(g_features_scaled_bw_sideways,strcat(path,'0deg bw.bmp'),'bmp');

g_bw=g_features_scaled_bw;
g_fill=zeros([ndim mdim]);

%--------------------Gap Filling-------------------------------------------
%--------------(Horizontal & Vertical)-------------------------------------
g_fill=imdilate(g_bw,strel('line',2,90));
g_fill=imerode(g_fill,strel('line',2,90));
g_fill=imerode(g_fill,strel('line',2,90));
g_fill=imdilate(g_fill,strel('disk',1,0));
g_fill=imerode(g_fill,strel('line',2,90));
imwrite(g_bw,strcat(path,'90deg 01 binarized.bmp'),'bmp');
imwrite(g_fill,strcat(path,'90deg 02 filled.bmp'),'bmp');

%--------------------Data line Removal-------------------------------------
g_dataclear=imerode(g_fill,strel('line',7,0));
g_dataclear=imdilate(g_dataclear,strel('line',2,0));

imwrite(g_dataclear,strcat(path,'90deg 03 dataclear.bmp'),'bmp');
% =========================================================================
% END of image enchancement
% =========================================================================
% #########################################################################

% =========================================================================
% START of property(SCALELINE) Approximation
% =========================================================================
%--------------------Initialization-------------------------------------
g_counter=zeros([ndim mdim]);
g_counter=padarray(g_counter, [d d], 0);
g_counter=double(g_counter);

%--------------------Approximation-------------------------------------
for m=1+d:mdim-1-d,
prev=g_fill(1+d,m);
counter=0;
    for n=1+d:ndim-1-d,
        counter=counter+1;
        record(n,m)=counter;
        if (g_fill(n,m)~=prev)
            if (prev==1)
                g_counter(n-1,m)=counter;
            else 
                g_counter(n,m)=counter;
            end
            prev=g_fill(n,m);    
            counter=0;
        end
    end
end
% figure, imshow(g_counter/64);figure(gcf);
% image shows scaleline positions + weight,
% where weight=the distance from upper(previous) scale line
% *weight=1,2,3...  ;probably the bottom of some scale line
%       =    +-12   ;probably the top of a scale line
%       =   16+     ;probably scaleline,below some missing scaleline(s)

for m=1+d:mdim+2*d,
     scaleline_count_per_column(m-d)=nnz(g_counter(:,m));
end     %nnz=number-of-nonzero elements,
% count the number of possible scalelines represented by each column
% low count = column at position of hourmarker ?
% medium count = ...at pos of garbage/unreadable piece
% high count = should contain information on every scale line present

candidate_scaleline_colnr = find(scaleline_count_per_column==max(scaleline_count_per_column),1, 'first');
% find("where some array satisfies logic"  ,  "how many items?"  , "which element to return if multiples found" )
% returns the column number where most scalelines dwell
scaleline_imageaprox(:,1)=g_counter(:,candidate_scaleline_colnr);
% numel(scaleline_imageaprox)
% nnz(scaleline_imageaprox)
scaleline_imageaprox = imresize(scaleline_imageaprox, [ndim mdim], 'nearest');

scaleline_imageaprox=im2bw(scaleline_imageaprox, 0.1);
% create full size image of scalelines by using the candidate as template
% binarize scale line image, since counter left elements > [0:1]
imwrite(scaleline_imageaprox,strcat(path,'scaleline_imageaprox.bmp'),'bmp');

sl_fillaprox=imdilate(scaleline_imageaprox,strel('line',4,90));
sl_fillaprox=imdilate(sl_fillaprox,strel('line',4,90));
sl_fillaprox=imerode(sl_fillaprox,strel('line',4,90));
% fill gaps between scaleline markers(from g_counter)



sl_skelaprox=imerode(sl_fillaprox,strel('line',2,90));
sl_skelaprox=imerode(sl_skelaprox,strel('line',2,90));
sl_skelaprox=bwmorph(sl_skelaprox,'thin',5);
imwrite(sl_skelaprox,strcat(path,'sl_skelaprox.bmp'),'bmp');


num_of_SL=nnz(sl_skelaprox(:,candidate_scaleline_colnr));
c_last_SL=0;
p_between_sl=round(ndim/num_of_SL);

for n=1:ndim,       %-----add missing SL based on average SL seperation
    c_last_SL=c_last_SL+1;
    if (sl_skelaprox(n,candidate_scaleline_colnr)>0)
        c_last_SL=0;
    end
    
    if (c_last_SL>(p_between_sl)*1.5)
        sl_skelaprox(n-round(c_last_SL*0.333),:)=1;
        c_last_SL=0;
    end
end
imwrite(sl_skelaprox,strcat(path,'90deg 05_skelaprox.bmp'),'bmp');
%==========================================================================
%--------------------End of Scaleline approximation------------------------
%==========================================================================
% #########################################################################
%==========================================================================
% START of scaleline REGRESSION
%==========================================================================


%--------------------Scaleline sampler-------------------------------------
sl_pop=g_dataclear; %(population/sample from statistics) sl_pop is the full
                    %dataset from where samples will be taken
                    c_SL=candidate_scaleline_colnr;
                    c_SL=15;
                    sl_skel=sl_skelaprox;
                    imwrite(sl_skel,strcat(path,'SL_skeleton.bmp'),'bmp');
                    
                   
sl_pop=sl_skel;
seek_range=round(p_between_sl/2);
sl_count=0;
sl_regressed=double(zeros([ndim mdim]));
    for n=1:ndim,
        if (sl_skelaprox(n,c_SL)>0)    
            sl_count=sl_count+1;
            
            obs_count=0;
            for m=1:mdim,
                for i=-seek_range:seek_range,   %refer 5 lines upward
                    if ((n+i>0) & (n+i<ndim) & (sl_pop(n+i,m)>0)) % (n+i>0) border control (1+(-6)=-5 outofbounds matrix)
                        obs_count=obs_count+1;                    % (n+i<ndim) same thing
                        resp(obs_count,sl_count)=n+i;             %
                        pred(obs_count,sl_count)=m;
                    end
                end
            end
            if (size(resp,2)<sl_count)
                sl_count=sl_count+1;
                break;
            end
            
%--------------------Draw regression lines part 1-------------------------------------
            index=find(pred(:,sl_count)==0,1, 'first')-1;
            [r,slope,b]=regression(pred(1:index,sl_count),resp(1:index,sl_count),'one');
            for m=1:round(mdim/3),
                if (round(slope*m+b)<1)
                        %border exception
                elseif (round(slope*m+b)>ndim)
                         %border exception
                else                                % !! ONLY draw when inside image ranges
                    sl_deci=slope*m+b-floor(slope*m+b);
                    if (sl_deci<0.5),
                        sl_regressed(round(slope*m+b),m)=1-sl_deci; %bilinear interpolation of offcenter pixel values
                        sl_regressed(round(slope*m+b)+1,m)=sl_deci;
                    else
                        sl_regressed(round(slope*m+b),m)=sl_deci;
                        sl_regressed(round(slope*m+b)-1,m)=1-sl_deci;
                    end
                end                                         %ignore the temptation to write outside matrix bounds
            end

%--------------------Draw regression lines part 2-------------------------------------
        index=find(pred(:,sl_count)==0,1, 'first')-1;
        [r,slope,b]=regression(pred(1:index,sl_count),resp(1:index,sl_count),'one');
        for m=round(mdim/3)+1:round(2*mdim/3),
            if (round(slope*m+b)<1)
                    %border exception
            elseif (round(slope*m+b)>ndim)
                     %border exception
            else                                % !! ONLY draw when inside image ranges
                sl_deci=slope*m+b-floor(slope*m+b);
                if (sl_deci<0.5),
                    sl_regressed(round(slope*m+b),m)=1-sl_deci; %bilinear interpolation of offcenter pixel values
                    sl_regressed(round(slope*m+b)+1,m)=sl_deci;
                else
                    sl_regressed(round(slope*m+b),m)=sl_deci;
                    sl_regressed(round(slope*m+b)-1,m)=1-sl_deci;
                end
            end                                         %ignore the temptation to write outside matrix bounds
        end
            
%--------------------Draw regression lines part 3-------------------------------------
            index=find(pred(:,sl_count)==0,1, 'first')-1;
            [r,slope,b]=regression(pred(1:index,sl_count),resp(1:index,sl_count),'one');
            for m=round(2*mdim/3)+1:mdim,
                if (round(slope*m+b)<1)
                        %border exception
                elseif (round(slope*m+b)>ndim)
                         %border exception
                else                                % !! ONLY draw when inside image ranges
                    sl_deci=slope*m+b-floor(slope*m+b);
                    if (sl_deci<0.5),
                        sl_regressed(round(slope*m+b),m)=1-sl_deci; %bilinear interpolation of offcenter pixel values
                        sl_regressed(round(slope*m+b)+1,m)=sl_deci;
                    else
                        sl_regressed(round(slope*m+b),m)=sl_deci;
                        sl_regressed(round(slope*m+b)-1,m)=1-sl_deci;
                    end
                end                                         %ignore the temptation to write outside matrix bounds
            end            
            
        end
    end

sl_regressed_old=sl_regressed;    
c_last_SL=0;    
for n=1:ndim,       %-----add missing SL based on average SL seperation
    c_last_SL=c_last_SL+1;
    if (sl_regressed_old(n,1)>0)
        c_last_SL=0;
    end
    
    if (c_last_SL>(p_between_sl)*1.5)
        index=find(pred(:,sl_count)==0,1, 'first')-1;
        [r,slope,b]=regression(pred(1:index,1),resp(1:index,1),'one');
        for m=1:mdim,
            if (n+round(slope*m)>0 & n+round(slope*m)<ndim)     % !! ONLY draw when inside image ranges
                sl_deci=slope*m+b-floor(slope*m+b);
                if (sl_deci<0.5),
                    sl_regressed(round(slope*m+b)-round(c_last_SL*0.333)+n,m)=1-sl_deci;
                    sl_regressed(round(slope*m+b)+1-round(c_last_SL*0.333)+n,m)=sl_deci;
                else
                    sl_regressed(round(slope*m+b)-round(c_last_SL*0.333)+n,m)=sl_deci;
                    sl_regressed(round(slope*m+b)-1-round(c_last_SL*0.333)+n,m)=1-sl_deci;
                end
            end                                         %ignore the temptation to write outside matrix bounds
        end
        c_last_SL=0;
    end
end
    imwrite(sl_regressed,strcat(path,'0deg 08SL_regressed.bmp'),'bmp');
    imwrite(f(:,:,3)/255-sl_regressed,strcat(path,'reg fitted to source.bmp'),'bmp');    
%==========================================================================
% END of scaleline REGRESSION
%==========================================================================    
% #########################################################################
% =========================================================================
% START of property(HOURLINE) finder
% =========================================================================
%---------------Initialization---------------------------------------------
hourmarker=zeros([ndim mdim]);
hm_sort=sort(scaleline_count_per_column);
hm_quantile=quantile(hm_sort,hm_index);

for m=1:mdim,
    if (scaleline_count_per_column(m)<hm_quantile)
        hourmarker(1:ndim,m)=1;
    end
end
hourmarker=imdilate(hourmarker,strel('diamond',9));
imwrite(hourmarker,strcat(path,'90deg 06hourmarker.bmp'),'bmp');

%---------------Find center of hourmarkers---------------------------------
hm_skel=bwmorph(hourmarker,'thin',500);
hm_skel=imresize(hm_skel(100:200,:), [ndim mdim]) ;
hm_final=imdilate(hm_skel,strel('diamond',5));

imwrite(hm_final,strcat(path,'90deg 06 1hm_final.bmp'),'bmp');
% =========================================================================
% END of property(HOURLINE) finder
% =========================================================================
% #########################################################################
% =========================================================================
% START of property(DATALINE) finder
% =========================================================================

%===============Line Reduction=====================
%---------------Preliminaries----------------------------------
[filename, pathname] = ...
     uigetfile({'*.bmp;*.jpg;',...
 'Picture Files (*.bmp,*.jpg)';},'Select source image for processing');
data=rgb2hsv(double((imread(strcat(strcat(pathname,'\'),filename) ))));


%---------------Preliminaries----------------------------------
data_bw = double(im2bw(data(:,:,3)/255,0.99 ));
data_bw = medfilt2(data_bw, [5 5]);
data_bw = imcomplement(data_bw);
dimension=size(data_bw);
ndimLINE=dimension(1);
mdimLINE=dimension(2);

data_bw = bwmorph(data_bw, 'spur', 15);
data_bw = medfilt2(data_bw, [2 2]);
data_line=zeros([ndimLINE mdimLINE]);

%---------------Determining Mean Line Thicknes----------------------------------
   count=0;
for m=1:mdimLINE-1,
            for n=1:ndimLINE-1,
                if (data_bw(n,m) > 0)
                count = count+1;
                    if (data_bw(n+1,m)==0),
                        data_line(round(n-(count/2)),m) = 1;
                        count=0;
                    end
                end
            end
            count=0;           
end

%---------------Seperating data--------------------------------------------

data_line_heat=data_line;
data_line_heat(250:ndimLINE,1:mdimLINE)=0;

data_line_bar=data_line;
data_line_bar(1:550,1:mdimLINE)=0;

data_line(1:250,1:mdimLINE)=0;
data_line(550:ndimLINE,1:mdimLINE)=0;

% =========================================================================
% END of property(DATALINE) finder


% =========================================================================
% START of property(SPROCKET) finder
% =========================================================================

%===============Finding Sprocket Locations=====================

%---------------Preliminaries----------------------------------
g_bw = im2bw(f(:,:,3)/255,0.95 );
graythresh(f(:,:,3)/255);
g_bw = im2bw(f(:,:,3)/255,0.95 );
g_bw=double(g_bw);
H=fspecial('gaussian', 5);
g_sprocketTEMP = g_bw;
for n=60:ndim-60,
            for m=1:mdim,
                g_sprocketTEMP(n,m) = 0;
            end
end
%---------------Filtering--------------------------------------
    g_sprocket=imfilter(g_sprocketTEMP,H,'replicate');

%---------------Erosion----------------------------------------
for i=1:10,
    g_sprocket=imerode(g_sprocketTEMP,strel('diamond',1));
    g_sprocketTEMP = g_sprocket;
end
 
g_sprocket = im2bw(g_sprocket,0.01);

%---------------Closing & Clustering---------------------------
g_sprocket = imclose(g_sprocket, strel('line',5, 90));
[labels, num] = bwlabel(g_sprocket, 8);
g_sprocket_locations = logical(zeros(size(g_sprocket)));

for i = 1:num
   [r, c] = find(labels==i); 
   x = round(mean(c));
   y = round(mean(r));
   g_sprocket_locations(y,x) = 1;
end

%---------------Ideal Sprocket----------------------------------
%for all pixels
%    if pixel >0 and is no more than 60 pixels away from a top or bottom
%    edge
%        all pixels within 22 distance => white. from 1-25 (if distance
%        from pixel <23 = white.)
%        delete white pixels in rows more than 14, but less than 40 pixels to the left or right of centre pixel

g_ideal_sprockets = logical(zeros(size(g_bw)));

for m=23:mdim-1,
    for n=23:ndim-1,
        if (g_sprocket_locations(n,m)>0) && ((n<60) || (ndim-n <60))
            for i=n-23:n+23
                for j=m-23:m+23
                distance = sqrt(((i)-n)^2+((j)-m)^2);
                    if  ((j>m-14) && (j<m+14)) && (distance <= 22)
                        g_ideal_sprockets(i,j)=1;
                    end
                end
            end
        end
    end
end       
imwrite((g_ideal_sprockets),(strcat(path,'2) Ideal Sprockets.bmp')),'bmp');
%---------------Shifting and un-unsharpening------------------------
g_ideal_sprockets = g_ideal_sprockets(1:ndim,1:mdim);
g_sprocket_padded = padarray(g_ideal_sprockets, [4 4], 0);
g_sprocket_shifted = g_sprocket_padded(4:ndim+3,4:mdim+3);
unshadow=xor(g_ideal_sprockets,g_sprocket_shifted);
shadow = g_ideal_sprockets-unshadow;
shadow_blur=imfilter(double(shadow),fspecial('gaussian',5,2));

g_ideal_sprockets_shadows=g_ideal_sprockets+shadow_blur*1.2;
sprock_test = g_ideal_sprockets_shadows(1:ndim,1:mdim);

imwrite((sprock_test),(strcat(path,'2) Distorted Sprockets.bmp')),'bmp');

% =========================================================================
% END of property(SPROCKET) finder
%==========================================================================


%==========================================================================
% START of background TEXTURE
%==========================================================================
%---------------Scalelines and Hourmarkers combined------------------------
temp1=imdilate(sl_regressed,strel('diamond',2));
temp2=imdilate(hm_skel,strel('diamond',5));
tex_grid=temp1 | temp2 | g_ideal_sprockets(1:ndim,1:mdim);
%---------------Find texture between grid----------------------------------
tex_piece=xor(f(:,:,3) , tex_grid);
tex_piece=double(tex_piece.*f(:,:,3))/255;

imwrite(f(:,:,3)/255-tex_grid,strcat(path,'tex_grid to source.bmp'),'bmp'); 

%---------------Filling of empty spaces between texture blocks-------------
% tex_tile=tex_piece;
% pit=[1,1];  %direction "data pit", where to find samples of colors
% 
% for i=1:ndim,                         %works 100% , but slow
%     if (50*pit(1)+i+11 > ndim) 
%         pit(1)=pit(1)*(-1);
%     end
%     
%     for j=1:mdim,
%         if (50*pit(2)+j+11) > mdim 
%             pit(2)=pit(2)*(-1);
%         end
%             
%         if (tex_tile(i,j)==0),
%             for adj=1:50,   %adjuster
%                 rdm=random('unid',10,2,1)-1;    %random offset of where to search for pixels
%                 
%                 if (i+adj*pit(1) +rdm(1) ~= 1:ndim)
%                     rdm(1)=rdm(1)*(-1);
%                 end
%                 if (j+adj*pit(2) +rdm(2) ~= 1:mdim)
%                     rdm(2)=rdm(2)*(-1);
%                 end
%                 if (tex_piece(i+adj*pit(1)+rdm(1) , j+adj*pit(2)+rdm(2)) ~=0)
%                     tex_tile(i,j)=tex_piece(i+adj*pit(1) +rdm(1) , j+adj*pit(2) +rdm(2));
%                     break;
%                 end
%             end
%         end
%     end
%     pit(2)=1;
% end
% imwrite(tex_tile,strcat(path,'tex_tile.bmp'),'bmp');
%---------------Shortcut---------------------------------------------------
tex_tile=imread(strcat(path,'tex_tile.bmp'));       %!!!->=="SHORTCUT TO SKIP ABOVE ALGO(saves 20 minutes)"==<-!!!
tex_tile=double(tex_tile)/255;                      %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
imwrite(tex_tile,strcat(path,'tex_tile.bmp'),'bmp'); 

%---------------Collect blant texture without detail-----------------------
tex_tile_contrast=histeq(tex_tile);
tex_tile_contrast=-tex_tile_contrast+tex_tile;
% tex_tile_contrast=tex_tile_contrast+min(min(tex_tile_contrast)) 
tex_ghost=log10(tex_tile_contrast+10)*0.5/log10(11);
tex_ghost=tex_ghost*(mean(mean(f(:,:,3)/255))/mean(mean(nonzeros(tex_piece))));
imwrite(tex_ghost,strcat(path,'tex_ghost.bmp'),'bmp');


%---------------Collect contrast details-----------------------------------
tex_viginete_b=im2bw(tex_ghost,mean(quantile(tex_ghost,0.92)));

tex_viginete_w=1-im2bw(tex_ghost,mean(quantile(tex_ghost,0.08)));


%---------------Remove ghostly data line from texture----------------------
for n=1:ndim,   
    tex_viginete_w(n,:)=tex_viginete_w(n,:)*(1-sin(n*pi/ndim)); %sinus gradient fill
end

%---------------Create motion blur of contrast detail----------------------
H = fspecial('motion', 50, 90);
tex_viginete_blur_b=-(imfilter(double(tex_viginete_b),H));
tex_viginete_blur_w=(imfilter(double(tex_viginete_w),H));
tex_final=tex_ghost+tex_viginete_blur_b/16-tex_viginete_b/32+tex_viginete_blur_w/16+tex_viginete_w/32;

%==========================================================================
% END of background TEXTURE
%==========================================================================

%==========================================================================
% START of DISTORTION
%==========================================================================
%---------------Rotate hourmarkers perpendicular to scalelines-------------
slope_array=zeros([size(pred,2) 1]);
for i=1:size(pred,2),
    index=find(pred(:,i)==0,1, 'first')-1;
    [r,slope_array(i),b]=regression(pred(1:index,i),resp(1:index,i),'one');
end

hm_distort=imrotate(double(hm_final),-radtodeg(atan(mean(slope_array))),'bilinear','crop');

%---------------Distort cosmic ray data------------------------------------
H = fspecial('motion', 8, 90);
data_line_final=imfilter(double(data_line),H);
data_line_final=log10(data_line_final+1);
data_line_final=data_line_final*1/max(max(data_line_final));
data_line_final=log10(data_line_final+1);
data_line_final=data_line_final*1/mean(quantile(nonzeros(data_line_final),0.001));
data_line_final=min(data_line_final,1);
data_line_final=imfilter(data_line_final,H);

%---------------Distort barometric data------------------------------------
data_line_bar_final=imfilter(double(data_line_bar),H);
data_line_bar_final=log10(data_line_bar_final+1);
data_line_bar_final=data_line_bar_final*1/max(max(data_line_bar_final));
data_line_bar_final=log10(data_line_bar_final+1);
data_line_bar_final=data_line_bar_final*1/mean(quantile(nonzeros(data_line_bar_final),0.001));
data_line_bar_final=min(data_line_bar_final,1);
data_line_bar_final=imfilter(data_line_bar_final,H);

%---------------Distort temperature data-----------------------------------
H = fspecial('motion', 5, 90);
data_line_heat_final=imfilter(double(data_line_heat),H);
data_line_heat_final=log10(data_line_heat_final+1);
data_line_heat_final=data_line_heat_final*1/max(max(data_line_heat_final));
data_line_heat_final=log10(data_line_heat_final+1);
data_line_heat_final=data_line_heat_final*1/mean(quantile(nonzeros(data_line_heat_final),0.001));
data_line_heat_final=min(data_line_heat_final,1);
data_line_heat_final=imfilter(data_line_heat_final,H);

%---------------Intergration of component images---------------------------
H = fspecial('gaussian', 5, 1);
grid=imfilter(sl_regressed+(hm_distort-(im2bw(sl_regressed,0.01)&hm_distort))/4,H);

ADL=data_line_bar_final+data_line_heat_final+data_line_final+grid*5;
H = fspecial('gaussian', 2, 2);
ADL_blur_1=imfilter(ADL,H);
H = fspecial('gaussian', 3, 2);
ADL_blur_2=imfilter(ADL,H);

ADL_ring=ADL_blur_2-ADL_blur_1-ADL;
ADL_ring=max(ADL_ring,0);

%==========================================================================
% END of DISTORTION
%==========================================================================

%==========================================================================
% START of OPTICAL NOISE
%==========================================================================
%---------------Finalizing result------------------------------------------   
H = fspecial('gaussian', 5, 0.1);
data_distort_d=(data_line_final-0.4*imfilter((im2bw(sl_regressed,0.01)&im2bw(data_line_final,0.01)),H))/3;
data_distort_t=1.2*(data_line_heat_final-0.4*imfilter((im2bw(sl_regressed,0.01)&im2bw(data_line_heat_final,0.01)),H))/3;
data_distort_b=(data_line_bar_final-0.4*imfilter((im2bw(sl_regressed,0.01)&im2bw(data_line_bar_final,0.01)),H))/6;


groundTruth=grid+data_line_final+data_line_heat_final-data_line_bar_final+sprock_test+(ADL_ring*0.3);
synthetic=tex_final+grid+data_line_final+data_line_heat_final-data_line_bar_final+sprock_test+(ADL_ring*0.3);

imwrite(groundTruth,strcat(path,'groundTruth.bmp'),'bmp');
imwrite(synthetic,strcat(path,'synthetic.bmp'),'bmp');
imwrite(imnoise(synthetic,'poisson'),strcat(path,'synthetic_noise.bmp'),'bmp');
%==========================================================================
% END of OPTICAL NOISE
%==========================================================================

%==========================================================================
% END of Cross section
%==========================================================================
display('Completed.');
display('----------');
display('- - - - - ');