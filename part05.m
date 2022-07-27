
%########################################################################
%######## This part of the code 
%             a) Filter out sigma outside 1-3 pixel range 
%             b) Select Amplititide larger than 100 and get 
%             c) calculate the distance between adjacent peaks
% 
%             d) Combine the peaks if the distance is smaller than 1/2 of sigma
%             e) Calculate the new centres
%             f) Calculate the distances between the selected peaks
%##########################################################################



disp('Part6 and add TabulatedData1good to Tabsorted...')


cutoffDist1=SigmaSMch1Avg./2; % change to match PSF for each Fluorophore
cutoffDist2=SigmaSMch2Avg./2;

fileSet=1;%

TabdataSorted1=TabulatedData1good;
TabdataSorted2=TabulatedData2good;


% F I L T E R sigma and Amp

% if sigma is between 1-3 and amp is larger
% than 100 get the index and add the value to Tab2 table


index1=find(TabdataSorted1(:,4)<=3 & (TabdataSorted1(:,4)>=1) & ...
    TabdataSorted1(:,1)>=100);%find the one that sum of RedGreen is equal or smaller than max possible
    Tab_lim_Sig_Amp_Ch1=TabdataSorted1(index1,:);% apply the index to 'TabdataSorted1' make new table

index2=find((TabdataSorted2(:,4)<=3) & (TabdataSorted2(:,4)>=1)& ...
    TabdataSorted2(:,1)>=100);%find the one that sum of RedGreen is equal or smaller than max possible
    Tab_lim_Sig_Amp_Ch2=TabdataSorted2(index2,:);%


% % C A L C U L A T E  D I S T A N C E  between adjasent peaks
for h=2:size(Tab_lim_Sig_Amp_Ch1,1)

    Tab_lim_Sig_Amp_Ch1(h,7)=((Tab_lim_Sig_Amp_Ch1(h,2)-Tab_lim_Sig_Amp_Ch1(h-1,2))^2 ...
        +(Tab_lim_Sig_Amp_Ch1(h,3)-Tab_lim_Sig_Amp_Ch1(h-1,3))^2)^0.5;
end

for h=2:size( Tab_lim_Sig_Amp_Ch2,1)

    Tab_lim_Sig_Amp_Ch2(h,7)=(( Tab_lim_Sig_Amp_Ch2(h,2)- Tab_lim_Sig_Amp_Ch2(h-1,2))^2 ...
        +( Tab_lim_Sig_Amp_Ch2(h,3)- Tab_lim_Sig_Amp_Ch2(h-1,3))^2)^0.5;
end



% Groups the peaks too close to each other
%CHANNEL 1

clearvars dummy* indcombine groupedpeaks Tab3DataCombined1

dummy1=1;
totpeaks=size(Tab_lim_Sig_Amp_Ch1,1);

i=0;
groupedpeaks1=[];
while dummy1<=totpeaks
    dummy2=dummy1+1;
    indcombine=dummy1;
    while  dummy2<=totpeaks && Tab_lim_Sig_Amp_Ch1(dummy2,7)<cutoffDist1
        indcombine=[indcombine;dummy2];
        dummy2=dummy2+1;
    end

    dummy1=dummy2;
    i=i+1;
    groupedpeaks1{i}=indcombine;
    Tab3DataCombined1(i,1)=sum(Tab_lim_Sig_Amp_Ch1(indcombine,1));
    Tab3DataCombined1(i,2)=sum(Tab_lim_Sig_Amp_Ch1(indcombine,1).* ...
        Tab_lim_Sig_Amp_Ch1(indcombine,2))./sum(Tab_lim_Sig_Amp_Ch1(indcombine,1));
    Tab3DataCombined1(i,3)=sum(Tab_lim_Sig_Amp_Ch1(indcombine,1).* ...
        Tab_lim_Sig_Amp_Ch1(indcombine,3))./sum(Tab_lim_Sig_Amp_Ch1(indcombine,1));

    mindist=((Tab3DataCombined1(i,2)-Tab_lim_Sig_Amp_Ch1(min(indcombine),2))^2 + ...
        (Tab3DataCombined1(i,3)-Tab_lim_Sig_Amp_Ch1(min(indcombine),3))^2)^0.5;
    maxdist=((Tab3DataCombined1(i,2)-Tab_lim_Sig_Amp_Ch1(max(indcombine),2))^2 + ...
        (Tab3DataCombined1(i,3)-Tab_lim_Sig_Amp_Ch1(max(indcombine),3))^2)^0.5;

    Tab3DataCombined1(i,4)=(mindist+maxdist)/2+(Tab_lim_Sig_Amp_Ch1(min(indcombine),4)+ ...
        Tab_lim_Sig_Amp_Ch1(max(indcombine),4))/2; %fix the sigma: if mix sig+sig2/2
    Tab3DataCombined1(i,5)=sum(Tab_lim_Sig_Amp_Ch1(indcombine,5)); % leave or add the bkgs
    Tab3DataCombined1(i,6)=sum(Tab_lim_Sig_Amp_Ch1(indcombine,6));
    Tab3DataCombined1(1,7)=0; %intiate the distance col


end

% Re-calculating the distances

for h=2:size(Tab3DataCombined1,1)

    Tab3DataCombined1(h,7)=((Tab3DataCombined1(h,2)-Tab3DataCombined1(h-1,2))^2 ...
     +(Tab3DataCombined1(h,3)-Tab3DataCombined1(h-1,3))^2)^0.5;

end

T3_Tab3DataCombined1 = array2table(Tab3DataCombined1,...
    'VariableNames',{'Amp','x','y','Sigma','bkgInt','TmNoCom','NewDistance'});

% Groups the peaks too close to each other
%CHANNEL 2

clearvars dummy* indcombine groupedpeaks Tab3DataCombined2

dummy1=1;
totpeaks=size(Tab_lim_Sig_Amp_Ch2,1);

i=0;
groupedpeaks2=[];
while dummy1<=totpeaks
    dummy2=dummy1+1;
    indcombine=dummy1;
    while  dummy2<=totpeaks && Tab_lim_Sig_Amp_Ch2(dummy2,7)<cutoffDist2;
        indcombine=[indcombine;dummy2];
        dummy2=dummy2+1;
    end
    dummy1=dummy2;
    i=i+1;
    groupedpeaks2{i}=indcombine;

    Tab3DataCombined2(i,1)=sum(Tab_lim_Sig_Amp_Ch2(indcombine,1));
    Tab3DataCombined2(i,2)=sum(Tab_lim_Sig_Amp_Ch2(indcombine,1).* ...
        Tab_lim_Sig_Amp_Ch2(indcombine,2))./sum(Tab_lim_Sig_Amp_Ch2(indcombine,1));
    Tab3DataCombined2(i,3)=sum(Tab_lim_Sig_Amp_Ch2(indcombine,1).* ...
        Tab_lim_Sig_Amp_Ch2(indcombine,3))./sum(Tab_lim_Sig_Amp_Ch2(indcombine,1));

    mindist=((Tab3DataCombined2(i,2)-Tab_lim_Sig_Amp_Ch2(min(indcombine),2))^2 + ...
        (Tab3DataCombined2(i,3)-Tab_lim_Sig_Amp_Ch2(min(indcombine),3))^2)^0.5;
    maxdist=((Tab3DataCombined2(i,2)-Tab_lim_Sig_Amp_Ch2(max(indcombine),2))^2 + ...
        (Tab3DataCombined2(i,3)-Tab_lim_Sig_Amp_Ch2(max(indcombine),3))^2)^0.5;

    Tab3DataCombined2(i,4)=(mindist+maxdist)/2+(Tab_lim_Sig_Amp_Ch2(min(indcombine),4)+ ...
        Tab_lim_Sig_Amp_Ch2(max(indcombine),4))/2;
    Tab3DataCombined2(i,5)=sum(Tab_lim_Sig_Amp_Ch2(indcombine,5));
    Tab3DataCombined2(i,6)=sum(Tab_lim_Sig_Amp_Ch2(indcombine,6));

    Tab3DataCombined2(i,7)=0;

end


% RE-calculating the distance between Red and green peaks
for h=2:size(Tab3DataCombined2,1)
    Tab3DataCombined2(h,7)=((Tab3DataCombined2(h,2)-Tab3DataCombined2(h-1,2))^2+ ...
        (Tab3DataCombined2(h,3)-Tab3DataCombined2(h-1,3))^2)^0.5;
end

%Put the value in a titled table
T3_Tab3DataCombined2 = array2table(Tab3DataCombined2,...
    'VariableNames',{'Amp','x','y','Sigma','bkgInt','TmNoCom','NewDistance'});

% ratio of Tm A compare to Tm B on one filament
RatioTmPerFilament=mean(T3_Tab3DataCombined1.TmNoCom)/mean(T3_Tab3DataCombined2.TmNoCom);


