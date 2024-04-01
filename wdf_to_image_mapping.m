addpath('./Matlab/')
filename='./map/S6T_1.wdf';
wdf = WdfReader(strcat(filename,'.wdf'), 'rb');
[WhiteLightImage,XCoords,YCoords]=wdf.GetWLImage();
imwrite(WhiteLightImage,strcat(filename,"_",...
    sprintf('%.f',XCoords(1)),"_",...
    sprintf('%.f',XCoords(length(XCoords))),"_",...
    sprintf('%.f',length(XCoords)),"_",...
    sprintf('%.f',YCoords(1)),"_",...
    sprintf('%.f',YCoords(length(YCoords))),"_",...
    sprintf('%.f',length(YCoords)),".png"));
maplength=wdf.Count;
mapping=wdf.GetSpectra(1,maplength);
[x,y]=wdf.GetOriginCoords(1,maplength);
title={'X','Y'};
title=string([title,num2cell(flip(wdf.GetXList(),2))]);
T = array2table([x,y,flip(mapping,2)]);
T.Properties.VariableNames=title;
writetable(T,strcat(filename,'.csv'))
wdf.Close()