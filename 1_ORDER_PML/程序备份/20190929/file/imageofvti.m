clc,clear;
%% 参数设置
i_nt=75;  %时间采样
i_nz=400; %空间采样Nz
i_nx=400; %空间采样Nx
c_sei=zeros(i_nz,i_nx,i_nt);

%% 文件读取
[a,state0]=dataRead('Pall.dat',i_nt*i_nx*i_nz,1);

%将文件排成一个三维矩阵
i=1;
for t=1:i_nt
    for x=1:i_nx
        for z=1:i_nz
            c_sei(z,x,t)=a(i);
            i=i+1;
        end
    end
end

cmin=min(min(c_sei(:)));
cmax=max(max(c_sei(:)))  ;
% cmin=-1*cmax;

% 
% %% 文件保存为图片并生成GIF
% path_in='./pic/';  %路径设置
% for i=1:50
%     figure(1),
%     %     imagesc(c_sei(:,:,i)),colormap(gray),axis equal,axis([0 400 0 400]),hold on
%     imagesc(c_sei(:,:,i),[cmin,cmax] ),hold on
%     colormap(gray)
%     axis([0 400 0 400])
%     axis equal
% 
%     %     caxis([cmin cmax])
%     xlabel('Nx'),ylabel('Nz'),
%     %     saveas(gca,[path_in,num2str(i),'Vx.jpg']); %把单个图片保存在path_in路径下
%     drawnow;
%     MakeGif('Vx3.gif',i)  %保存为动态图
%     hold on
%     close
% end
% hold off
% 
% %    mesh(c_sei(:,:,i)),axis([0 400 0 400 -4*10^-13  4*10^-13])
% 
% 
% 
% 
%% 子函数
%% （子函数01）二进制文件读取
function [data,state]=dataRead(fileName,M,N)
fid=fopen(fileName,'rb');
if(fid>0)
    [data,count]=fread(fid,[M,N],'float');
    if(count==size(data,1)*size(data,2))
        state=1;
    else
        state=-1;
    end
end
fclose(fid);
end
%% （子函数02）gif动态图保存
function MakeGif(filename,i)
f = getframe(gcf);
imind = frame2im(f);
[imind,cm] = rgb2ind(imind,256);
if i==1
    imwrite(imind,cm,filename,'gif',...
        'Loopcount',inf,'DelayTime',0.01);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01);
end
end
