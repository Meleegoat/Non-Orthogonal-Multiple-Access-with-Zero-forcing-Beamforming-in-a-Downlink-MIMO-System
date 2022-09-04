clear
%number of Tx antenna
for round=1:3
p=0.95;
nTx=2;
N=100;%number of users
R_NOMAtot=[];
R_NConvtot=[];
R_NOMA=[];
R_Weaktot=[];
R_Convtot=[];
R_Conventional=[];
for ite=1:800
    m=1;
for nUser=[4,16,36,64,128]
radius=500;
h=[];
r = abs(2*radius*rand(1,nUser)-radius)/1e3;
pl_dB = 128.1+37.6*log10(r);
pl = 10.^(0.1*pl_dB);
h=1/sqrt(2)*(randn(nTx,nUser)+1i*randn(nTx,nUser));
Bw=4.32*1e6;
N0= Bw*10^(0.1*(-169-30));
P0=10^((43-30)/10);
h=h.*pl.^(-1/2);
pairs=[];
Gain=[];
Gaincopy=[];
CorrelationArray=[];
i=1;
for iUser = 1 : nUser
    Gain(iUser)=norm(h(:,iUser));
    for jUser=iUser+1:nUser
    correlation=abs(dot(h(:,iUser),h(:,jUser)))/(norm(h(:,iUser))*norm(h(:,jUser)));
    CorrelationArray(iUser,jUser)=correlation;
    if correlation>p
       pairs(1,i)=iUser;
       pairs(2,i)=jUser;
       GainDif=abs(norm(h(:,iUser))-norm(h(:,jUser)));
       %gain difference=||h1|-|h2||
       pairs(3,i)=GainDif;
       i=i+1;
     end
     end
end
Gaincopy=Gain;
pairscopy=pairs;
k=1;
GainDif=[];
selected=[];
if isempty(pairs)==0
for j=1:2
    GainDif=pairs(3,:);
    if max(GainDif)~=0
    index=find(GainDif==max(GainDif));
    %create a new array storing selected pairs
    selected(1,k)=pairs(1,index);
    selected(2,k)=pairs(2,index);
    selected(3,k)=pairs(3,index);     %gain difference
    selected(4,k)=Gain(selected(1,k));%chanel gain for selected user 1
    selected(5,k)=Gain(selected(2,k));%chanel gain for selected user 2
    pairs(:,index)=0;
    %delete the used Tx and Rx
    Common1=find(pairs(1,:)==selected(1,k));
    Common2=find(pairs(1,:)==selected(2,k));
    Common3=find(pairs(2,:)==selected(1,k));
    Common4=find(pairs(2,:)==selected(2,k));
    Common=[Common1,Common2,Common3,Common4];
    pairs(1,Common)=0;
    pairs(2,Common)=0;
    pairs(3,Common)=0;
    k=k+1;
    end
end
end
%identify the user use for beamforming vector
stronguser=zeros(1,2);
weakuser=zeros(1,2);
H=[];
H_weak=[];
W=[];
sz2=size(selected);
%Conventional beamform clustering
for i=1:2
    convuser1(i)=find(Gaincopy==max(Gaincopy));
    Gaincopy(convuser1(i))=0; 
    convuser2(i)=find(Gaincopy==max(Gaincopy));
    Gaincopy(convuser2(i))=0;
    H_conv1(:,i)=h(:,convuser1(i));
    H_conv2(:,i)=h(:,convuser2(i));
end
%NOMA beamform clustering
for i=1:2
    if i<=sz2(2)
    if selected(4,i)>selected(5,i)
        stronguser(i)=selected(1,i);
        Gain(stronguser(i))=0;
        weakuser(i)=selected(2,i);
        Gain(weakuser(i))=0;

    else
        stronguser(i)=selected(2,i);
        Gain(stronguser(i))=0;
        weakuser(i)=selected(1,i);
        Gain(weakuser(i))=0;
    end
    else
        stronguser(i)=find(Gain==max(Gain));
        Gain(stronguser(i))=0; 
        weakuser(i)=find(Gain==max(Gain));
        Gain(weakuser(i))=0;
    end 
    H(:,i)=h(:,stronguser(i));
    H_weak(:,i)=h(:,weakuser(i));
end
W=[];
W_1=[];
W=H*(H'*H)^-1;
W_1=H_weak*(H_weak'*H_weak)^-1;
ii=1;
switch round
    case 1 
    a1=1/sqrt(1+P0/N0*norm(H_weak(:,1))^2)-((sqrt(1+P0/N0*norm(H_weak(:,1))^2)-1)*(norm(H_weak(:,1)'*W(:,2))^2)*(P0/N0)+1)/(norm(H_weak(:,1)'*W(:,1))^2*(P0/N0)*sqrt(1+P0/N0*norm(H_weak(:,1))^2));
    if a1>=0
        a1=a1;
    else
        a1=0.10;
    end
        a2=1/sqrt(1+P0/N0*norm(H_weak(:,2))^2)-((sqrt(1+P0/N0*norm(H_weak(:,2))^2)-1)*(norm(H_weak(:,2)'*W(:,1))^2)*(P0/N0)+1)/(norm(H_weak(:,2)'*W(:,2))^2*(P0/N0)*sqrt(1+P0/N0*norm(H_weak(:,2))^2));
     if a2>=0
        a2=a2;
    else
        a2=0.10;
        
     end
     a1copy(ii)=a1;
     ii=ii+1;
    case 2
        a1=0.15;
        a2=0.15;
    case 3
        a1=0.25;
        a2=0.25;
end
a3=1-a1;
a4=1-a2;
SINR_str1=norm(H(:,1))^2*a1*P0/N0;
SINR_str2=norm(H(:,2))^2*a2*P0/N0;
SINR_W1=norm(H_weak(:,1)'*W(:,1))^2*a3*P0/(norm(H_weak(:,1)'*W(:,1))^2*a1*P0+norm(H_weak(:,1)'*W(:,2))^2*P0+N0);
SINR_W2=norm(H_weak(:,2)'*W(:,2))^2*a4*P0/(norm(H_weak(:,2)'*W(:,2))^2*a2*P0+norm(H_weak(:,2)'*W(:,1))^2*P0+N0);
SINR_Nconv1=norm(H_weak(:,1))^2*P0/N0;
SINR_Nconv2=norm(H_weak(:,2))^2*P0/N0;
%Conventional beamforming
SINR_conv1=norm(H_conv1(:,1))^2*P0/N0;
SINR_conv2=norm(H_conv1(:,2))^2*P0/N0;
SINR_conv3=norm(H_conv2(:,1))^2*P0/N0;
SINR_conv4=norm(H_conv2(:,2))^2*P0/N0;
switch sz2(2)
    case 0
        R_NOMA(m)=1/2*1/1e6*Bw*(log2(1+SINR_conv1)+log2(1+SINR_conv2)+log2(1+SINR_conv3)+log2(1+SINR_conv4));
    case 1
        R_NOMA(m)=1/1e6*Bw*(log2(1+SINR_str1)+log2(1+SINR_W1))+1/1e6*1/2*Bw*(log2(1+SINR_str2/a2)+log2(1+SINR_Nconv2));
    case 2
        R_NOMA(m)=1/1e6*Bw*(log2(1+SINR_str1)+log2(1+SINR_str2)+log2(1+SINR_W1)+log2(1+SINR_W2));
end
 if ite==1 
 R_NOMAtot(m)=R_NOMA(m);
 else 
 R_NOMAtot(m)=(R_NOMAtot(m)*(ite-1)+R_NOMA(m))/ite;
 end
 m=m+1;
end
end
if round==1
name1=plot([4,16,36,64,128],R_NOMAtot,'-*');
hold on

end
if round==2
name3=plot([4,16,36,64,128],R_NOMAtot,'-o');
hold on

end
if round==3
name5=plot([4,16,36,64,128],R_NOMAtot,'-^');
hold on

legend([name1,name3,name5],["NOMA-BF(dynamic \alpha)","NOMA-BF(fixed \alpha=0.15)","NOMA-BF(fixed \alpha=0.25)"]);
end
end
xlabel('Number of users');
ylabel('Sum Capacity [Mbps]');

ax = gca; 
ax.FontSize = 16; 