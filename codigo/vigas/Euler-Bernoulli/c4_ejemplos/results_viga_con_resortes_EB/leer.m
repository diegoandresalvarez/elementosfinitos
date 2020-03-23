fid = fopen('xx_EB.txt');  x = str2num(fscanf(fid,'%c')); fclose(fid);
fid = fopen('Vx_EB.txt');  V = str2num(fscanf(fid,'%c')); fclose(fid);
fid = fopen('Mx_EB.txt');  M = str2num(fscanf(fid,'%c')); fclose(fid);
fid = fopen('tx_EB.txt');  t = str2num(fscanf(fid,'%c')); fclose(fid);
fid = fopen('vxx_EB.txt'); v = str2num(fscanf(fid,'%c')); fclose(fid);

figure(1)
subplot(2,1,1);
hold on;
plot(x, v, 'r.');
subplot(2,1,2);
hold on;
plot(x, t, 'r.');

figure(2)
subplot(2,1,1);
hold on;
plot(x, M, 'r.');
subplot(2,1,2);
hold on;
plot(x, V, 'r.');