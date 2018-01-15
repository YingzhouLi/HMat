%%
%General stuff starting the script, you can ignore this part

%Loading the current directory and all folders contained
% pathstring = pwd;                                   %Get the address of the current working directory
% parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
% [~,n] = find(~cellfun(@isempty,strfind(parts,'Yingzhou'))); %finding the scripts root directory and its location (n)
% addpath(genpath(strjoin(parts(1,1:n),'\')));                      %adding all folders to the path that are in this dir 
% cd(pathstring)                                       %jumping back to the current dir
addpath('../src');
%getting rid of any arrays/figures hanging around
clear; close all

%%
%My Default equation system

%Loading 4 dense coefficient matricies and two vectors. 
load 'Data_CoeffMatriciesUnConc.mat'  

%Arranging these as dense matrix 'A'
A=  [NdTn,SdTn;
     NdTs,SdTs];

%Now making into a large column vector 'B'
B = [Tn;Ts];  

size=1:numel(B);

%Running to find unknowns 'D' 
D = A\B;


%%
%The same system using the Hmat script

%Converting these to Hmats (from func down base of script)
[NdTn,SdTn,NdTs,SdTs]=CreateHMat(NdTn,SdTn,NdTs,SdTs);

%Appending these 
AA=[NdTn,SdTn;NdTs,SdTs];
% % A1=horzcat(NdTn,SdTn);A2=horzcat( NdTs,SdTs );AA=vertcat(A1,A2);

%Running to find unknowns 'DD'
DD = AA\B; %DD = mldivide(AA,B);

% %If I use the pre appended matrix A it works fine. 
% [A]=CreateHMat(A);
% DD = A\B; %DD = mldivide(AA,B);

%Drawing a comparison figure
figure;
plot(size,D,'k');
hold on
plot(size,DD,'.');









