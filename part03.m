
disp(['Part 03: fitting 2D Gaussian constrcting TabulatedData1 and 2; ' ...
    'Amp, x, y, sigma, mean bkg '])

tolerance = config.tolerance;
PSF=config.PSF; % a possible value to give a code a starting point to work on

convrems1=[];
convrems2=[];
domain_red = config.domain_rad ;

%%%% CHANNEL 1 %%domain_red%%%%%%%%%
% This part of the code picks the highest peak->fit 2D Gaussian-> Deflate 
% (deconvolve)->repeate until noise (set by tolerance parameters) level 
% is reached

t=1;
I1=imcrop(imageData1(:,:,t),rectfilament);
class(I1);
I1 = double(I1);
class(I1);

Iback1=double (imcrop(imageData1(:,:,t),rectnoise));
noise1=mode(I1((I1.*mask2)>0))+std(Iback1(:))*randn(size(mask2));
Iorig=(noise1.*(~mask2)+I1.*mask2);

stdbefore=std(Iorig(:));
stdnow=std(Iback1(:));
j=0;


try
    while stdnow-stdbefore<tolerance
        stdbefore=std(Iorig(:));
        j=j+1;
        Y = fft2(Iorig);

        % remove specific peak from the signal define the delta filter out
        % Find the approximate parameter based on symmetric fit
        [pospeaky, pospeakx]=find(max(Iorig(:))==Iorig);
        x_dom =(pospeakx-domain_red):(pospeakx+domain_red);
        y_dom =(pospeaky-domain_red):(pospeaky+domain_red);

        working_chunk = Iorig(y_dom, x_dom);

        %if working_chunk
        options = optimset('Display', 'off');
        %fit Gaussian around the peak neighbourhood to extract the peak position and
        %amplitude...
        a0(1)=max(Iorig(:));
        a0(2)=pospeakx;
        a0(3)=pospeaky;
        a0(4)=PSF/2;
        a0(5)=mean(Iback(:));

        [xdomain,ydomain]=meshgrid(x_dom,y_dom);
        weightfn=exp(-(xdomain-pospeakx).^2/((size(working_chunk,1))^2) ...
            ).*exp(-(ydomain-pospeaky).^2/((size(working_chunk,2))^2));
        ydata=working_chunk.*weightfn;
        lb=[];
        ub=[];

        [a, ~,~, exitFlag]=lsqcurvefit(@gauss_2DSymmetric,a0, ...
            {x_dom,y_dom},ydata,lb,ub,options);

        convrem=gauss_2DSymmetric(a,{1:size(Iorig,2),1:size(Iorig,1)});

        Y2 = fft2(convrem);

        %subtract the peak from the Fourier signal
        newfsig=Y-Y2;

        % invert the new signal and plot
        newsig = ifft2(newfsig);
        TabulatedData1(j,:)=a;
        Iorig=newsig;
        stdnow=std(Iorig(:));

        Iorigs1(:,:,j)=Iorig;
        newsigs1(:,:,j)=newsig;
        convrems1(:,:,j)=convrem;

    end
catch
end

%%%% CHANNEL 2 %%%%%%%%%%%
% a) Select the signal with maximum intensity (highest peak) 
% b) Fit it with a 2D Gaussian-> 
% c) Deflate (deconvolve) from the remaining signal ->
% d) Repeat a-c until noise (set by tolerance parameters) level is reached

I2=imcrop(imageData2(:,:,t),rectfilament); I2 = double(I2);

Iback2=double(imcrop(imageData2(:,:,t),rectnoise));
noise2=mode(I2((I2.*mask2)>0))+std(Iback2(:))*randn(size(mask2));

Iorig=(noise2.*(~mask2)+I2.*mask2);

stdbefore=std(Iorig(:));
stdnow=std(Iback2(:));
j=0;


try
    while stdnow-stdbefore<tolerance
        stdbefore=std(Iorig(:));
        j=j+1;
        Y = fft2(Iorig);
        % remove specific peak from the signal define the delta filter out
        [pospeaky pospeakx]=find(max(Iorig(:))==Iorig);
        x_dom = (pospeakx-domain_red):(pospeakx+domain_red);
        y_dom = (pospeaky-domain_red):(pospeaky+domain_red);

        working_chunk = Iorig(y_dom, x_dom);
        options = optimset('Display', 'off');
        %fit Gaussian around the peak neighbourhood to extract the peak position and
        %amplitude...
        a0(1)=max(Iorig(:));
        a0(2)=pospeakx;
        a0(3)=pospeaky;
        a0(4)=PSF/2;
        a0(5)=mean(Iback2(:));        
        
        [xdomain,ydomain]=meshgrid(x_dom,y_dom);
        weightfn=exp(-(xdomain-pospeakx).^2/((size(working_chunk,1)/1)^2) ...
            ).*exp(-(ydomain-pospeaky).^2/((size(working_chunk,2)/1)^2));
        ydata=working_chunk.*weightfn;
        lb=[];
        ub=[];
       
        [a, ~,~, exitFlag]=lsqcurvefit(@gauss_2DSymmetric,a0,{x_dom,y_dom},ydata,lb,ub);
        %exitFlag

        convrem=gauss_2DSymmetric(a,{1:size(Iorig,2),1:size(Iorig,1)});

        Y2 = fft2(convrem);
        %subtract the peak from the Fourier signal
        newfsig=Y-Y2;
        % invert the new signal and plot
        newsig = ifft2(newfsig);
        TabulatedData2(j,:)=a;
        Iorig=newsig;

        stdnow=std(Iorig(:));

        Iorigs2(:,:,j)=Iorig;
        newsigs2(:,:,j)=newsig;
        convrems2(:,:,j)=convrem;

    end
catch
end

