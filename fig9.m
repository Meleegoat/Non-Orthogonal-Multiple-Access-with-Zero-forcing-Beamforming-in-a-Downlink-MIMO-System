clear
%number of Tx antenna
for p= 0.85
nTx=2;
N=100;%number of users
R_NOMAtot=[];
R_NEXH=[];
R_NOMA=[];
R_Weaktot=[];
R_Convtot=[];
R_Conventional=[];
for ite=1:300
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
       pairs(4,i)=log2(1+Gain(iUser)^2*P0/N0)+log2(1+norm(h(:,jUser))^2*P0/N0);
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
    selected(6,k)=log2(1+selected(4,k)^2*P0/N0)+log2(1+selected(5,k)^2*P0/N0);

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
H_1=[];
H_weak1=[];
H_2=[];
H_weak2=[];
W=[];
sz2=size(selected);
sz3=size(pairscopy);
rdnnumber=[];
%%random select
rdnnumber=randperm(nUser,4);
H_2(:,1)=h(:,rdnnumber(1));
H_weak2(:,1)=h(:,rdnnumber(2));
H_2(:,2)=h(:,rdnnumber(3));
H_weak2(:,2)=h(:,rdnnumber(4));

%exhaustive search
for i=1:2 
    if nUser~=4
  if i<=sz2(2) && max(pairscopy(4,:))~=0

  number=find(pairscopy(4,:)==max(pairscopy(4,:)));
  if Gain(pairscopy(1,number))> Gain(pairscopy(2,number))
    H_1(:,i)=h(:,pairscopy(1,number));
    H_weak1(:,i)=h(:,pairscopy(2,number));
    pairscopy(4,number)=0;


  else
    H_1(:,i)=h(:,pairscopy(2,number));
    H_weak1(:,i)=h(:,pairscopy(1,number));
    pairscopy(4,number)=0;


  end
    Common5=find(pairscopy(1,:)==pairscopy(1,number));
    Common6=find(pairscopy(1,:)==pairscopy(2,number));
    Common7=find(pairscopy(2,:)==pairscopy(1,number));
    Common8=find(pairscopy(2,:)==pairscopy(2,number));
    Common=[Common5,Common6,Common7,Common8];
    pairscopy(1,Common)=0;
    pairscopy(2,Common)=0;
    pairscopy(3,Common)=0;
    pairscopy(4,Common)=0;

  else
        number1=find(Gaincopy==max(Gaincopy));
        H_1(:,i)=h(:,number1);
        Gaincopy(number1)=0; 
        number2=find(Gaincopy==max(Gaincopy));
        Gaincopy(number2)=0;
        H_weak1(:,i)=h(:,number1);
  end
    else
        rdnnumber=randperm(nUser,4);
H_1(:,1)=H_2(:,1);
H_weak1(:,1)=H_weak2(:,1);
H_1(:,2)=H_2(:,2);
H_weak1(:,2)=H_weak2(:,2);
break
    end
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
W_2=[];
W=H*(H'*H)^-1;
W_1=H_1*(H_1'*H_1)^-1;
W_2=H_2*(H_2'*H_2)^-1;
a1=0.85;
%RANDOM
SINR_RAD1=norm(H_2(:,1))^2*(1-a1)*P0/N0;
SINR_RAD2=norm(H_2(:,2))^2*(1-a1)*P0/N0;
SINR_RAD3=norm(H_weak2(:,1)'*W_2(:,1))^2*a1*P0/(norm(H_weak2(:,1)'*W_2(:,1))^2*(1-a1)*P0+norm(H_weak2(:,1)'*W_2(:,2))^2*P0+N0);
SINR_RAD4=norm(H_weak2(:,2)'*W_2(:,2))^2*a1*P0/(norm(H_weak2(:,2)'*W_2(:,2))^2*(1-a1)*P0+norm(H_weak2(:,2)'*W_2(:,1))^2*P0+N0);
%exhautive search 
SINR_exh1=norm(H_1(:,1))^2*(1-a1)*P0/N0;
SINR_exh2=norm(H_1(:,2))^2*(1-a1)*P0/N0;
SINR_exh3=norm(H_weak1(:,1)'*W_1(:,1))^2*a1*P0/(norm(H_weak1(:,1)'*W_1(:,1))^2*(1-a1)*P0+norm(H_weak1(:,1)'*W_1(:,2))^2*P0+N0);
SINR_exh4=norm(H_weak1(:,2)'*W_1(:,2))^2*a1*P0/(norm(H_weak1(:,2)'*W_1(:,2))^2*(1-a1)*P0+norm(H_weak1(:,2)'*W_1(:,1))^2*P0+N0);
%%%
SINR_str1=norm(H(:,1))^2*(1-a1)*P0/N0;
SINR_str2=norm(H(:,2))^2*(1-a1)*P0/N0;
SINR_W1=norm(H_weak(:,1)'*W(:,1))^2*a1*P0/(norm(H_weak(:,1)'*W(:,1))^2*(1-a1)*P0+norm(H_weak(:,1)'*W(:,2))^2*P0+N0);
SINR_W2=norm(H_weak(:,2)'*W(:,2))^2*a1*P0/(norm(H_weak(:,2)'*W(:,2))^2*(1-a1)*P0+norm(H_weak(:,2)'*W(:,1))^2*P0+N0);
SINR_Nconv1=norm(H_weak(:,1))^2*P0/N0;
SINR_Nconv2=norm(H_weak(:,2))^2*P0/N0;
%Conventional beamforming
SINR_conv1=norm(H(:,1))^2*P0/N0;
SINR_conv2=norm(H(:,2))^2*P0/N0;
SINR_conv3=norm(H_weak(:,1))^2*P0/N0;
SINR_conv4=norm(H_weak(:,2))^2*P0/N0;
switch sz2(2)
    case 0
        R_NOMA(m)=1/2*1/1e6*Bw*(log2(1+SINR_conv1)+log2(1+SINR_conv2)+log2(1+SINR_conv3)+log2(1+SINR_conv4));
    case 1
        R_NOMA(m)=1/1e6*Bw*(log2(1+SINR_str1)+log2(1+SINR_W1))+1/1e6*1/2*Bw*(log2(1+SINR_str2/(1-a1))+log2(1+SINR_Nconv2));
    case 2
        R_NOMA(m)=1/1e6*Bw*(log2(1+SINR_str1)+log2(1+SINR_str2)+log2(1+SINR_W1)+log2(1+SINR_W2));
end
R_EXH(m)=1/1e6*Bw*(log2(1+SINR_exh1)+log2(1+SINR_exh2)+log2(1+SINR_exh3)+log2(1+SINR_exh4));
R_Conventional(m)=1/1e6*Bw*(log2(1+SINR_RAD1)+log2(1+SINR_RAD2)+log2(1+SINR_RAD3)+log2(1+SINR_RAD4));
 if ite==1 
 R_NOMAtot(m)=R_NOMA(m);
 R_NEXH(m)=R_EXH(m);
 R_Convtot(m)=R_Conventional(m);
 else 
 R_NOMAtot(m)=(R_NOMAtot(m)*(ite-1)+R_NOMA(m))/ite;
 R_NEXH(m)=(R_NEXH(m)*(ite-1)+R_EXH(m))/ite;
 R_Convtot(m)=(R_Convtot(m)*(ite-1)+R_Conventional(m))/ite;
 end
 m=m+1;
end
end
if p==0.85
name1=plot([4,16,36,64,128],R_NOMAtot,'-or','LineWidth',1);
hold on

name3=plot([4,16,36,64,128],R_Convtot,'--ok','LineWidth',1);
hold on
name2=plot([4,16,36,64,128],R_NEXH,'-.ob','LineWidth',1);
hold on
end

legend([name1,name2,name3],["NOMA-BF(Exhautive)","NOMA-BF","NOMA-BF(Random)"]);
xlabel('Number of users');
ylabel('Sum Capacity [Mbps]');
end
ax = gca; 
ax.FontSize = 16; 