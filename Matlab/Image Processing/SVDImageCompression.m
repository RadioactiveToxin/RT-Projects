%% SVD Image Compression
clear; close all; clc;

% First we use uigetfile to open a dialog box
[file,path]=uigetfile('C:\Users\Asus N552VW\Desktop\*.jpg; *.JPG; *.PNG; *.*','Select an Image File');
% And then fullfile gives the full path of the chosen file
selectedfile=fullfile(path,file);

% now we assign our inputs to M,A,R,G,B
M=imread(selectedfile); A=im2double(M);
R=A(:,:,1); [UR,SR,VR]=svd(R);
G=A(:,:,2); [UG,SG,VG]=svd(G);
B=A(:,:,3); [UB,SB,VB]=svd(B);
C=im2double(rgb2gray(M));
[m,n]=size(C);

% Then we set up the dialog box properties
prompt={"Enter The K Values Less Than "+string(min(m,n))+":"}; dlgtitle='Input'; dims=[1 60];
% inputdlg Creates dialog box to gather user input {inputdlg(prompt,dlgtitle,dims)}
Input=inputdlg(prompt,dlgtitle,dims);
K=str2num(string(Input));
K=K(K < min(m,n));

% R is the number of rows in the subplot
if (floor((length(K)+2)/4))*4==(length(K)+2)
    R=(length(K)+2)/4;
else
    R=(floor((length(K)+2)/4))+1;
end

figure("Name","SVD Image Compression","NumberTitle","off");
subplot(R,4,1), imshow(M); title("Original");
subplot(R,4,2), imshow(C); title("OriginalGrayscale");

plotnum=3;
for k=K
    X(:,:,1)=UR(:,1:k)*SR(1:k,1:k)* VR(:,1:k)';
    X(:,:,2)=UG(:,1:k)*SG(1:k,1:k)* VG(:,1:k)';
    X(:,:,3)=UB(:,1:k)*SB(1:k,1:k)* VB(:,1:k)';
    subplot(R,4,plotnum), imshow(X); title("Approximation."+string(plotnum)+",k="+string(k));
    plotnum=plotnum+1;
end
