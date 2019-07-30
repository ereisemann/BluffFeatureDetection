% ---- From Lauren Dunkin 02/2017
% ---- edited by Eve Eisemann 03/2017
% ---- edited by Michael Hartman 04/2017
% ---- edited by Eve Eisemann 04/2017, 05/2017, 04/2018

% clear,clc
% addpath(genpath('E:\Eisemann\Projects\SWG_Feature_Extraction\data\'))

%%%%%%%%%%%%%%%%%  USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%
% --- Adding folder and subfolders to search path --- %

% ---- Save a file? y/n ----%
savefile = 'y';

% ---- Make transect figures? y/n ---%
figures = 'n';

% --- Change this to skip transects - meters
% pSkip = 10;

% ---- Input Smoothing Parameter ---- %
% S = 20;
S = 1;
% ---- helps to minimize jaggedness (5-20) 

% ---- Cross shore limit (maximum distance from shoreline to look for
% bluff in meters)
Lim = 400;

% --- Bluff Percent Slope Minimum (decimal)
slopemin = 0.35;

% ---- slope maximum at top of bluff (used line ~311)
% ---- added 04/22/2019
slopemax = 0.15;

% ---- Minimum bluff height (elevation)
minht = 190; 

indir = 'F:\GreatLakes\';
outdir = 'F:\GreatLakes\Bluffs\MI_5_1_UPDATE\';
% outdir = 'F:\GreatLakes\Bluffs\Testing\';
SLbase = 'MI_5_1_10m_rerun_';
addpath(genpath(indir))
addpath(genpath(outdir))

% ---- Loading ALL transects in one, will find within each segemnt
transpath = ('F:\GreatLakes\Bluffs\Transects_II\MI_5_1_10m.shp');
trans = shaperead(transpath,'UseGeoCoords',true);
Ntrans = length(trans);
addpath(genpath('F:\GreatLakes\Bluffs\MI_5_1_UPDATE\'))

%--GRIDS--%
GRIDSpath = ([indir 'Post_11_12_13\Michigan_2012\2012_NCMP_WI_Michigan_BareEarth_1mGrid\']);
%GRIDSpath = ([indir 'Post_11_12_13\Huron_2013\2013_NCMP_MI_Huron_BareEarth_1mGrid\']);
GRIDSstruct = dir([GRIDSpath '*.TIF']);
GRIDS = {GRIDSstruct(:).name}; % making cell array of file names
GRIDS = GRIDS';
addpath(genpath(GRIDSpath))

%% CREATING LOOP FOR FILES

nn = length(GRIDS);
Transects = [];

% for jj = 15:22; %MI_1b
% for jj = 22:42  % MI_c2
% for jj = 106:108  %HU_1
% for jj = 72:95 % MI_3
% for jj = 48:55 %MI_3 44085*
for jj = 29:50
% for jj = 1:nn

   tiff = char(GRIDS(jj));
   blocknum = tiff(23:24); % CHANGE for different file names
%    blocknum = tiff(23:30); % For MI_3 44085*
%    blocknum = jj;
   blocknum
   
   tiffpath = strcat(GRIDSpath, tiff);

   %--Loading TIFF--%
   %%%%%%%%%%%%%%%%%%
   [z, r, bbox] = geotiffread(tiffpath);
   info = geotiffinfo(tiffpath);

   % ---- To view the grid 
%    mapshow(z,r,'DisplayType', 'mesh')

%% Clipping Grids 

[nR, nC] = size(z);
noDataValue = min( min( z ) );

% ---- Make location vector
x1 = ( r(3,1) : r(2,1) : r(3,1) + (nC-1)*r(2,1) );
y1 = ( r(3,2) : r(1,2) : r(3,2) + (nR-1)*r(1,2) )';

% ---- Make location matrices
x = repmat(x1,nR,1);
y = repmat(y1,1,nC);

% ---- Find empty space in top, bottom and sides of grids
rInd = [];
dirs = {'fwd','bwd'};
for iDir = 1:2
    
    switch dirs{iDir}
        case 'fwd'
            rind = 1:nR;
            cind = 1:nC;
        case 'bwd'
            rind = nR:-1:1;
            cind = nC:-1:1;
    end
    
    rInd = [];
    for R = rind
        dataInd = find( z(R,:) ~= noDataValue, 1, 'first' );
        if isempty( dataInd )
            rInd = [ rInd, R];
        end
    end
    
    cInd = [];
    for c = cind
        dataInd = find( z(:,c) ~= noDataValue, 1, 'first' );
        if isempty( dataInd )
            cInd = [ cInd, c];
        end
    end
end

% ---- Replace no data values with nan
z( z == noDataValue ) = NaN; 
z = double(z);   

    Zsn = [];
    Xsn = [];
    Ysn = [];
   
% -- looping through transects in segment ---%    
pInd = 1;    
toex = []; 
toey = [];
toez = [];

bluffz = [];
bluffx = [];
bluffy = [];

% --- Finding transects within segment
in_1 = [];
in = [];
for tt  = 1:Ntrans
    in_1 = inpolygon(trans(tt).Lon(1), trans(tt).Lat(1), bbox(:,1),bbox(:,2));
    in = vertcat(in, in_1);
end 

trans_ind = find(in  == 1);

% --- Cycling through only those transects 

for tt = min(trans_ind):max(trans_ind)
%  for tt = [7460, 7461, 7473, 7474, 7475, 7537, 7538, 7555, 7556]
% for tt = [6747, 6748, 6749, 6796, 6809, 6810, 6828, 6829]
% for tt = [5473:5476]
% for tt = 7555
     %----------NEW 20180712----------%
     transLon = trans(tt).Lon(:);
     transLat = trans(tt).Lat(:);
     transLon(isnan(transLon)) = [];
     transLat(isnan(transLat)) = [];

     ptAX = transLon(1);
     ptAY = transLat(1);
     ptBX = transLon(end);
     ptBY = transLat(end);
     %--------------------------------%
     
% -------- Replaced By the Above ------------- %     
%      [ptAX, ind] = max(trans(tt).Lon(:));
%        ptAY = trans(tt).Lat(ind);
%       [ptBX, ind] = min(trans(tt).Lon(:));
%        ptBY = trans(tt).Lat(ind);
% -------------------------------------------- %

        % ---- Looking at shore-perpendicular transects
         disp([num2str(tt),'/',num2str(Ntrans)])
        [ZZ, ~, YY, XX] = mapprofile(z, r, [ptAY, ptBY], [ptAX, ptBX]);
        test = isempty(ZZ);
        
        Xi = XX';
        Yi = YY';
        Zi = ZZ';
        
        if test == 1
            continue
% -------- Removed 20180712 ------------- %                 
%         else
%       % --- convert column vector to row vector
%             Yi = fliplr(YY');
%             Xi = fliplr(XX');
%             Zi = fliplr(ZZ');
% -------- Removed 20180712 ------------- %  
        end
                       
%% Extracting Features for each transect

    iProfile = tt;
    pInd = tt;

pos = find(Zi>176); % finding elevations below Lake Level (176 m)

X = Xi(pos); % excluding bathy data from X Y Z
Y = Yi(pos);
Z = Zi(pos);

% ---- Plots transects on map
% figure
% hold on
% yyaxis left
% title(['MI_1', num2str(iProfile)])
% plot(X,Y,'k-')
% plot(X(1), Y(1), 'r*')
% ---- Save transects to structure array
% Trans = struct('ID',iProfile,'Geometry','Line','X',X,'Y',Y);
% Transects = [Transects,Trans];

% % ---- replacing Nan with empty set
nanind = find(isnan(Z));
Z(nanind) = [];
X(nanind) = [];
Y(nanind) = [];
if length(Z)<20 % --- Remove if little data
    Z = [];
    X = [];
    Y = [];
continue % --- exits loops if empty
end

% ---- Creating meters cross-shore vector
format long
lat1 = Y(1); % -- Most shoreward point
lon1 = X(1);
lat2 = Y(end); % -- Most inland point
lon2 = X(end);
[arclen, ~ ] = distance (lat1, lon1, lat2, lon2); % ---- Find arc distance, degrees
pdist =  distdim (arclen, 'degrees', 'meters'); % ---- Convert to meters
CSM1 = (0 : pdist/(length(X)-1) : pdist); % --- Cross Shore Meters
% ---- zero is seaward edge of data

% ---- Smooth each profile z values 
z2 = smoothn(Z, S);

% ---- Trimming to search window - cross shore limit
xshorelim = find(CSM1<Lim);
CSM = CSM1(xshorelim);
z2 = z2(xshorelim);

%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%
if figures == 'y'
figure(tt)
hold on
yyaxis left
ylabel('Elevation (m)')
plot(CSM,z2,'b-')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Calculating first derivative of profile
D1 = diff(z2);
dCSM = CSM(1:end-1); %derivative is length minus one, new CSM

% ---- finding greatest slope
[dpks, dlocs] = findpeaks(D1,dCSM);

testpk = isempty(dpks);
if testpk == 1
    continue
end

maxpk = max(dpks); %peak value
maxloc = dlocs(find(dpks == max(dpks))); %cross shore location
maxlocind = find(CSM == maxloc); %FOR TOE

% ---- finding next trough
pkind = find(dCSM == maxloc);

test = length(dCSM(pkind:end));
if test < 3
    continue
end
  
[dtrough, dtroughloc] = findpeaks(-D1(pkind:end), dCSM(pkind:end));
choice = find(-dtrough < slopemax, 1, 'first'); %changed to slopemax to make sure the trough is enough of a dip in slope (20190423)

test = isempty(choice);
if test == 1
    continue
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% below section moved to within the c = choice selection loop below (20190424)
% %incase D1 is not long enough to do an average over the next 20 m
% if length(D1)<(find(dCSM == dtroughloc(choice))+20)
%     disp('past maximum inland distance')
%     continue 
% %     dd = length(D1)-troughind;
% %     m = mean(D1(troughind:troughind+dd));
% end

%incase not many troughs found
if length(dtrough(choice:(length(dtrough))))==1
    disp('only one trough selected')
    choicePlus = choice;
elseif length(dtrough(choice:(length(dtrough))))<10
% if length(dtrough(choice):dtrough(end))<10 
    choicePlus = (choice + (length(dtrough(choice):dtrough(end))));
    disp('fewer than 10 additional troughs')
else
    choicePlus = (choice + 10);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = choice:choicePlus %indices of other troughs, up to ten more, or other extent from above loop
    trough = -dtrough(c);
    troughloc = dtroughloc(c);
    troughind = find(dCSM == troughloc);
    %incase D1 is not long enough to do an average over the next 20 m
        if length(D1)<(find(dCSM == dtroughloc(c))+20)
            disp('past maximum inland distance')
            m = NaN;
%             continue
            break
        %     dd = length(D1)-troughind;
        %     m = mean(D1(troughind:troughind+dd));
        end
    m = mean(D1(troughind:troughind+20)); %added average slope across 25 m must be below maxslope (20190423)
    disp(m)
if(m <= slopemax)
    break 
end
end

test = isnan(m); % added to avoid fitting if booted from previous loop (20190424)
if test == 1
    continue
end

% ---- Fitting Curve to D1 between crest and trough
D1polyfit = polyfit(dCSM(pkind:troughind),D1(pkind:troughind),3); %might want to expand (+10) to give alittle more space for curve to cross below min slope 
CSMpoly = dCSM(pkind:troughind);
polyD1 = polyval(D1polyfit,CSMpoly); %for plotting 

% ---- Find spot slope drops off below slopemax
crestind = find(polyD1<slopemax,1,'first');
crestloc = CSMpoly(crestind);

test = isempty(crestind);
if test == 1
    continue
end

CSMind = find(CSM == crestloc); %finding index in entire transect

% ---- Placing Bluff Crest
if maxpk < slopemin
           bluff_z = nan;
           bluff_x = nan;
           bluff_y = nan;
            
elseif Z(CSMind) < minht
           bluff_z = nan;
           bluff_x = nan;
           bluff_y = nan;
            
else 
          bluff_z = Z(CSMind);
          bluff_x = X(CSMind);
          bluff_y = Y(CSMind);
          
          bluffx = horzcat(bluffx,bluff_x); 
          bluffy = horzcat(bluffy,bluff_y);
          bluffz = horzcat(bluffz,bluff_z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIGURE%%%%%%%%%%%%%%%%%
if figures == 'y'
% figure(tt)
% hold on
% yyaxis left
% ylabel('Elevation (m)')
% plot(CSM1,Z,'b-')
plot(CSM(CSMind), bluff_z, 'k*')
% plot(CSM, z2, 'c-')
yyaxis right
ylabel('Slope')
plot(dCSM, D1, 'r-')
xlabel('distance from shoreline (m)')
title(['Lk. Michigan 3, transect ', num2str(iProfile)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIGURE%%%%%%%%%%%%%%%%%

% ---- Finding Bluff Toe
D1toe = D1(1:maxlocind);
dCSMtoe = dCSM(1:maxlocind);
toeloc = find(D1toe < 0.05, 1, 'last'); %last time slope is shallow. 

test = isempty(toeloc);
if test == 1
    continue
end

toeCSMind = find(CSM == dCSMtoe(toeloc));

% ---- Placing Bluff Toe


testz = isempty(bluff_z);
testa = isempty(toeloc);

% if testz == 1
%      toe_x = NaN;
%      toe_y = NaN;
%      toe_z = NaN;
% end
%      
% if testa == 1 % --- if toe_loc empty
%      toe_x = NaN;
%      toe_y = NaN;
%      toe_z = NaN;
% else
%     toe_z = Z(toeCSMind);
%     toe_x = X(toeCSMind);
%     toe_y = Y(toeCSMind); 
%     
%     if testz == 1
%      toe_x = NaN;
%      toe_y = NaN;
%      toe_z = NaN;
%     end

    if maxpk < slopemin
           toe_z = nan;
           toe_x = nan;
           toe_y = nan;
            
    elseif Z(CSMind) < minht
           toe_z = nan;
           toe_x = nan;
           toe_y = nan;
    else
        toe_z = Z(toeCSMind);
        toe_x = X(toeCSMind);
        toe_y = Y(toeCSMind); 
     
        
%%%%%%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures == 'y'
yyaxis left
plot(CSM(toeCSMind), toe_z, 'ko')
hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
%     if testa == 1 % --- if toe_loc empty
%      toe_x = NaN;
%      toe_y = NaN;
%      toe_z = NaN;
%     end
    
        toex = horzcat(toex,toe_x); 
        toey = horzcat(toey,toe_y);
        toez = horzcat(toez,toe_z);
% end

    end
hold off
end

bluffx(bluffx == 0) = NaN;
bluffy(bluffy == 0) = NaN; % For some reason there are zeros
toex(toex == 0) =  NaN;
toey(toey == 0) = NaN;


%% Writing Files
if savefile == 'y'

    mpCrest = mappoint(bluffx',bluffy','Bluff_Crest', bluffz');
    test = isempty(mpCrest);
    if test == 1
        warning('No features detected in block.')
        continue 
    end
    shapewrite(mpCrest, [outdir 'Crest\' SLbase blocknum '_BluffCrest_02'])

    mpToe = mappoint(toex', toey', 'Dune_Toe', toez');
    shapewrite(mpToe, [outdir 'Toe\' SLbase blocknum '_BluffToe_02'])
    
% --- FOR SAVING TEST FILES --- %   

%     shapewrite(mpCrest, [outdir SLbase blocknum '_BluffCrest'])
% 
%     mpToe = mappoint(toex', toey', 'Dune_Toe', toez');
%     shapewrite(mpToe, [outdir SLbase blocknum '_BluffToe'])

% --- end if yes to save file
end
end
