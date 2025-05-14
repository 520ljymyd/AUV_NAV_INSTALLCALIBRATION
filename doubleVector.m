glvs
trj = trjfile('trj10ms.mat'); %获得avp, imu, avp0, wat, ts, repeats  

[nn, ts, nts] = nnts(2, trj.ts);
%inst = [0;0;0.5]*glv.min;  kod = 1;  qe = 0; dT = 0;  % od parameters
inst = [0;0;1.5]*glv.deg;  kod = 1;  qe = 0; dT = 0;  % od parameters
%安装偏差角度设置  inst[俯仰，横滚，航向]安装偏差角度
% 刻度系数误差为0，量化误差为0，里程计与IMU时间延迟为0
trjod = odsimu(trj, inst, kod, qe, dT, 0);
%使用 odsimu函数 输入真实的imu和avp信息、里程仪的参数，产生里程仪的测量数据
%希望获得DVL系的速度，积分得到DVL系的位移，得到DVL系统的空间矩阵
imuerr = imuerrset(0.01, 50, 0.001, 5);
%imuerr = imuerrset(0.001, 1);  %陀螺三轴零偏，加速计三轴零偏
imu = imuadderr(trjod.imu, imuerr);

davp = avperrset([0;0;0], 0, 0);  %姿态attitude 速度velocity，位置position误差设置
%davp = avperrset([0;0;0], 0, 0);
%dinst = [15;0;10]*glv.min; dkod = 0.2; 
%dinst = [0;0;0]*glv.min; dkod = 0.25; 
dinst = [0;0;1.5]*glv.deg; dkod = 0.001; 

[attsb, qnb] = alignsb(imu, avp(7:9));

%phi = [aa2phi(attsb,[0;0;0]), [[-imuerr.db(2);imuerr.db(1)]/glv.g0;-imuerr.eb(1)/(cos(avp(7))*glv.wie)]];

%trjod = odsimu(trj, inst, kod, qe, dT, 0);%DVL仿真