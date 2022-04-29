
clc;
clear all;
V = 'flicker.mp4';      %Video Name 
xyloObj = VideoReader(V);   %Using video reader reading video
   %Extracting frames
   T= xyloObj.NumberOfFrames            % Calculating number of frames
   for g=1:T
           p=read( xyloObj,g);          % Retrieve data from video
           if(g~=  xyloObj.NumberOfFrames)
                 J=read( xyloObj,g+1);
                 th=difference(p,J);        %To calculate histogram difference between two frames 
                 X(g)=th;
           end
   end
   %calculating mean and standard deviation and extracting frames
   mean=mean2(X)
   std=std2(X)
   threshold=std+mean*4
    k = 0 ; %is short-term flicker metric
   for g=1: T
       p =  read(xyloObj,g);
       if(g~=xyloObj.NumberOfFrames)
        J=   read(xyloObj,g+1);
        th=difference(p,J);
                if(th>mean)    % Greater than threshold select as a key frame     

                k(g)=light_flickermeter_metric_PstLM(reshape(J,1,[]), 48000);                       

                end 
       end
   end 
   X=1:783;
   Y=k(1:783);
figure
plot(X,Y)
title('Short Term Flicker Values for Key Frames')
xlabel('Frame Number')
ylabel('Short Term flicker Values')
