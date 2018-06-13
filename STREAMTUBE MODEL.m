% Reservoir Properties %
h=25;
phi=0.18; % Average Porosity
P =750; % Pressure drop(psi)
rw=1.5; % Effective Wellbore Radius (ft)

% Rock and Fluid Properties
kb=0.05; % Perm(darcy)
Sor=0.3;
Siw=0.3;
mu_o=1.62; % Oil Viscosity(cp)
mu_w=1; % Water Viscosity(cp)

alpha1=0.402; % Kro=alpha1*(1-Swd)^m
alpha2=0.248; % Krw=alpha2*Swd^n

m=2.06;
n=2.33;
A=40*43560; % Reservoir area, ft2

% Pore Volume
Vp=A*h*phi/5.615/8;


Sw=Siw:0.0002:(1-Sor);
Swd=(Sw-Siw)/(1-Siw-Sor);
Kro=alpha1*(1-Swd).^m;
Krw=alpha2*Swd.^n;
fw=1./(1+(mu_w/mu_o)*(Kro./Krw));

j=length(fw);
fw_dash=[0 1:(j-2) 0];
for i = 2:(j-1)
fw_dash(i)=(fw(i-1)-fw(i+1))/(Sw(i-1)-Sw(i+1));
end
slope=(fw-fw(1,1))./(Sw-Sw(1,1));

for k=j:-1:1
if (slope(k)-fw_dash(k))<=0
ind=k;
break
end
end

Swf=Sw(ind);
fw_Swf=fw(ind);
fw_dash_Swf=fw_dash(ind);
Swf_bar=Swf+(1-fw_Swf)/fw_dash_Swf;
Qi_bt=1/fw_dash_Swf;

ST1=[17.372, 1.506 .886  .631  .495 .397 .325 .286 .256 .231 .211 .196 .183 .173 .166 .161 .158 .155 .152 .150 .150 .152 .155 .158 .161 .166 .173 .183 .196 .211 .231 .256 .286 .325 .397 .495 .631 .886 1.506 17.372];
ST2=[17.880 1.558  .990  .639  .495 .397 .331 .282 .245 .220 .200 .185 .172 .162 .155 .150 .147 .144 .139 .124 .124 .139 .144 .147 .150 .155 .162 .172 .185 .200 .220 .245 .282 .331 .397 .495 .639 .900 1.558 17.880];
ST3=[19.620 1.696  .950  .670  .536 .440 .363 .306 .262 .229 .203 .182 .166 .152 .139 .128 .118 .110 .105 .101 .101 .105 .110 .118 .128 .139 .152 .166 .182 .203 .229 .262 .306 .363 .440 .536 .670 .950 1.696 19.620];
ST4=[36.582 2.920  1.848 1.197 .928 .768 .643 .535 .451 .378 .310 .253 .205 .167 .139 .121 .108 .096 .088 .082 .082 .088 .096 .108 .121 .139 .167 .205 .253 .310 .378 .451 .535 .643 .768 .928 1.197 1.848 2.920 36.582];
VF=[0.2237 0.2478 0.2909 0.2376];
tube=[ST1',ST2',ST3',ST4'];

for zz=1:4;
Qi=0:Qi_bt/40:Qi_bt;
lamdainv_Swf=1./(Kro(ind:end)/mu_o+Krw(ind:end)/mu_w);
lamdainv_ro=mu_o/Kro(1);
total=0;
for i=1:(length(lamdainv_Swf)-1)
total=total+(lamdainv_Swf(i)+lamdainv_Swf(i+1))*(fw_dash(ind+i)-fw_dash(ind+i-1))/2;
end
lamdainv_Swf_bar=total/(0-fw_dash_Swf);
lamdainv=lamdainv_ro+(lamdainv_Swf_bar-lamdainv_ro).*Qi(1:end)*fw_dash_Swf;

fw_dash_abt=fw_dash(ind:end);
Qi_abt=1./fw_dash_abt;

lamdainv_abt=1./(Kro(ind:ind+400)/mu_o+Krw(ind:ind+400)/mu_w);
lamdainv_bar_abt=1:length(Qi_abt);
for i=1:(length(Qi_abt))
    summ=0;
    for j=i:(length(fw_dash_abt)-1)
     summ=summ+(lamdainv_Swf(j+1)+lamdainv_Swf(j))*(fw_dash_abt(j+1)-fw_dash_abt(j))/2;
    end
    lamdainv_bar_abt(i)=summ/(fw_dash_abt(end)-fw_dash_abt(i));
end
Qi=[Qi(1:end),Qi_abt];
lamdainv=[lamdainv,lamdainv_bar_abt];

Qi_bbt=0:Qi_bt/40:Qi_bt;
qt_bbt=1:40;
for jj=1:40
    fw_dash_CELL=0:jj;
    for kk=2:(length(fw_dash_CELL))
fw_dash_CELL(kk)=1/Qi_bbt(jj+1)*(kk-1)/40;
    end
  indd_bbt=1:length(fw_dash_CELL);
 for i=1:length(fw_dash_CELL)
    mm=min(abs(fw_dash_CELL(i)-fw_dash_abt));
    for j=1:length(fw_dash_abt)
        if mm==min(abs(fw_dash_CELL(i)-fw_dash_abt(j)))
            indx_bbt=j;
        end
    end
    indd_bbt(i)=indx_bbt;
 end
ssss=0;
lamdainv_cell_bbt=1:(length(fw_dash_CELL)-1);
for j=1:(length(indd_bbt)-1)
ssss=0;
for i=indd_bbt(j+1):(indd_bbt(j)-1)
ssss=ssss+(lamdainv_Swf(i)+lamdainv_Swf(i+1))*(fw_dash_abt(i)-fw_dash_abt(i+1))/2;
end
lamdainv_cell_bbt(j)=ssss/(fw_dash_abt(indd_bbt(j+1))-fw_dash_abt(indd_bbt(j)));
end
lam=[];
lam(1:(40-length(lamdainv_cell_bbt)))=lamdainv_ro;
lambdainv_bbt_total=[lamdainv_cell_bbt,lam];
addition=tube(:,zz)'.*lambdainv_bbt_total;
qt_bbt(jj)=1.127*kb*P*h/sum(addition);
end
qt_bbt=[1.127*kb*h*P/sum(tube(:,zz)'*lamdainv_ro), qt_bbt];
t_bbt=0:(length(Qi_bbt)-1);
for i=2:(length(Qi_bbt))
     t_bbt(i)=t_bbt(i-1)+2*(Qi_bbt(i)-Qi_bbt(i-1))*(VF(zz)*Vp)/(qt_bbt(i)+qt_bbt(i-1))/1;
end

Qi_test=Qi_abt;
fw_dash_test=1./Qi_test;
qt=1:length(fw_dash_test);
for ii=1:length(fw_dash_test)
fw_dash_cell=linspace(fw_dash_test(ii)/40,fw_dash_test(ii),40);
 indd=1:40;
 for i=1:length(fw_dash_cell)
    nn=min(abs(fw_dash_cell(i)-fw_dash_abt));
    for j=1:length(fw_dash_abt)
        if nn==min(abs(fw_dash_cell(i)-fw_dash_abt(j)))
            indx=j;
        end
    end
    indd(i)=indx;
 end
 indd=[length(fw_dash_abt), indd];
sss=0;
lamdainv_cell=1:40;
for j=1:(length(indd)-1)
sss=0;
for i=indd(j+1):(indd(j)-1)
sss=sss+(lamdainv_Swf(i)+lamdainv_Swf(i+1))*(fw_dash_abt(i)-fw_dash_abt(i+1))/2;
end
lamdainv_cell(j)=sss/(fw_dash_abt(indd(j+1))-fw_dash_abt(indd(j)));
end
summation=tube(:,zz)'.*lamdainv_cell;
qt(ii)=1.127*kb*P*h/sum(summation);
end
t=[t_bbt(end), 2:length(Qi_abt)];
for i=2:length(Qi_abt)
     t(i)=t(i-1)+2*(Qi_abt(i)-Qi_abt(i-1))*(VF(zz)*Vp)/(qt(i)+qt(i-1))/1;
end
temp1=fw(ind:end).*qt;
temp2(1:length(Qi_bbt))=0;
qo_abt=(1-fw(ind:end)).*qt;
Np_bbt=Qi_bbt;
Np_abt=[Qi_bt,2:length(Qi_abt)];
for i=1:(length(Qi_abt)-1)
Np_abt(i+1)=Np_abt(i)+qo_abt(i)*(t(i+1)-t(i))/(VF(zz)*Vp);
end
Np_total(:,zz)=[Np_bbt, Np_abt(2:400)]';
qo_total(:,zz)=[qt_bbt,qo_abt(2:400)]';
t_total(:,zz)=[t_bbt, t(2:400)]';
Qi_total(:,zz)=[Qi_bbt, Qi_abt(2:400)]';
qt_total(:,zz)=[qt_bbt, qt(2:400)]';
qw_total(:,zz)=[temp2, temp1(2:400)]';
end

time=[0:50:8000]';

Qi_ST1=interp1(t_total(:,1),Qi_total(:,1),time);
Qi_ST2=interp1(t_total(:,2),Qi_total(:,2),time);
Qi_ST3=interp1(t_total(:,3),Qi_total(:,3),time);
Qi_ST4=interp1(t_total(:,4),Qi_total(:,4),time);

qt_ST1=interp1(t_total(:,1),qt_total(:,1),time);
qt_ST2=interp1(t_total(:,2),qt_total(:,2),time);
qt_ST3=interp1(t_total(:,3),qt_total(:,3),time);
qt_ST4=interp1(t_total(:,4),Qi_total(:,4),time);

qo_ST1=interp1(t_total(:,1),qo_total(:,1),time);
qo_ST2=interp1(t_total(:,2),qo_total(:,2),time);
qo_ST3=interp1(t_total(:,3),qo_total(:,3),time);
qo_ST4=interp1(t_total(:,4),Qi_total(:,4),time);

qw_ST1=interp1(t_total(:,1),qw_total(:,1),time);
qw_ST2=interp1(t_total(:,2),qw_total(:,2),time);
qw_ST3=interp1(t_total(:,3),qw_total(:,3),time);
qw_ST4=interp1(t_total(:,4),Qi_total(:,4),time);

Np_ST1=interp1(t_total(:,1),Np_total(:,1),time);
Np_ST2=interp1(t_total(:,2),Np_total(:,2),time);
Np_ST3=interp1(t_total(:,3),Np_total(:,3),time);
Np_ST4=interp1(t_total(:,4),Qi_total(:,4),time);

qt_res=(qt_ST1+qt_ST2+qt_ST3+qt_ST4)*8;
qo_res=(qo_ST1+qo_ST2+qo_ST3+qt_ST4)*8;
qw_res=(qw_ST1+qw_ST2+qw_ST3+qt_ST4)*8;
Np_res=(Np_ST1*VF(1)+Np_ST2*VF(2)+Np_ST3*VF(3)+Np_ST4*VF(4));

Qi_res=(Qi_ST1*VF(1)+Qi_ST2*VF(2)+Qi_ST3*VF(3)+Qi_ST4*VF(4));
[time, qt_res, Qi_res,qo_res, qw_res, Np_res]
figure
plot(time, qw_res, 'b',time, qo_res, 'r',time, qt_res, 'k')
title('qt (BLACK) , qo (RED) and qw (BLUE) vs Time')
xlabel('Time (days)')
ylabel('Flow Rate (Bbls/D)')
figure
plot(Qi(1:400),lamdainv(1:400))
title(' Average Apparent Viscosity  vs  PV Injected ')
xlabel('Qi')
ylabel('Average Apparent Viscosity (cp)')
figure
plot(Qi_res,Np_res)
title(' Cummulative Oil Recovery  vs  PV Injected ')
xlabel('Qi (PV)')
ylabel('Np (Barrels)')
plot(time,qt_res)
title('Water Injection Rate vs Time')
xlabel('Time (Days)')
ylabel('Water Injection rate (Bbls/D)')