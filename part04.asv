
disp('Part 4: Removing peaks outside the filament ROI and return TabulatedData1good and 2..., still 5 columns')

% Initiate found1 and found2
found1=zeros(size(convrems1));
found2=zeros(size(convrems2));

for i=1:size(convrems1,3)
    found1(:,:,i)=sum(convrems1(:,:,1:i),3);
end

for i=1:size(convrems2,3)
    found2(:,:,i)=sum(convrems2(:,:,1:i),3);
end

for i=1:length(xCh1)-1
    seglength(i)=sqrt((xCh1(i+1)-xCh1(i))^2 +(yCh1(i+1)-yCh1(i))^2);
end

totallength=sum(seglength)*92.3/1000; %pixel size with 100x objective is 93.2 nm, ad with 10x objective is 91.3 nm

%get noise mean & correct raw images to have 0 offset off of the filament
avenoise1=mean(Iback1(:));
avenoise2=mean(Iback2(:));

I1prime=I1-avenoise1;
I2prime=I2-avenoise2;

% Calculate the total intensity inside the filament area of each channel
sumCh1=sum(nonzeros(mask2.*I1prime));
sumCh2=sum(nonzeros(mask2.*I2prime));

% ratio of the signal in Channel 1 or2 (green or red signal) to the total signal
Ratr= sumCh1/(sumCh1+sumCh2);
Ratg= sumCh2/(sumCh1+sumCh2);

% Remove peaks outside the filament
innonzeromask=find(mask); % find index of all non-zero elements in the mask area
dumind=0;

    for n=1:size(TabulatedData1,1); % get the number of rows
        if round(TabulatedData1(n,3))<size(mask2,1) && round(TabulatedData1(n,2))<size(mask2,2)... 
            && round(TabulatedData1(n,3))>0 && round(TabulatedData1(n,2))>0
            indpeak=sub2ind(size(mask2),round(TabulatedData1(n,3)), round(TabulatedData1(n,2))); % replave sub to index

            findEro =find(innonzeromask==indpeak); %find(innonzeromask==indpeak,1)

            X=isempty(findEro); %ISEMPTY(X) returns 1 if X is an empty array and 0 otherwise, if peak is inside or out side the mask

            if X==0   % if =0 means: peak is inside the mask include it in new tabulated data called TabulatedData1good
                dumind=dumind+1;
                TabulatedData1good(dumind,:)=TabulatedData1(n,:);
            elseif X==1
                TabulatedData1good(n,:)=TabulatedData1(n,:);
            end
        end
    end

    dumind=0;
    for n=1:size(TabulatedData2,1);
        if round(TabulatedData2(n,3))<size(mask2,1) && round(TabulatedData2(n,2))<size(mask2,2) ...
                && round(TabulatedData2(n,3))>0 && round(TabulatedData2(n,2))>0
            indpeak=sub2ind(size(mask2),round(TabulatedData2(n,3)), round(TabulatedData2(n,2)));

            findEro =find(innonzeromask==indpeak);
            X=isempty(findEro);
            if X==0
                dumind=dumind+1;
                TabulatedData2good(dumind,:)=TabulatedData2(n,:);
            elseif X==1
                TabulatedData2good(n,:)=TabulatedData2(n,:);
            end
        end
    end



%######### Integral of single peak ########################<<<<<<<------
      
SMIntegralPeakCh1Avg=AmpSMch1Avg*(SigmaSMch1Avg.^2)*pi;  %Ch1
SMIntegralPeakCh2Avg=AmpSMch2Avg*(SigmaSMch2Avg.^2)*pi;  %ch2 

%######### Integral of cluster ########################<<<<<<<------
IntegralPeak_ch1= TabulatedData1good(:,1).*((TabulatedData1good(:,4)).^2)*pi; %Ch1
IntegralPeak_ch2= TabulatedData2good(:,1).*((TabulatedData2good(:,4)).^2)*pi; %Ch2

%######### Number of Tpms per peak in row 6 ########################<<<<<<<-------
TabulatedData1good(:,6)=IntegralPeak_ch1/SMIntegralPeakCh1Avg; %Ch1
TabulatedData2good(:,6)=IntegralPeak_ch2/SMIntegralPeakCh2Avg; %Ch2





