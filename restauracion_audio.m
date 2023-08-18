close all
clear all
clc

pkg load signal
pkg load audio

function h = pasaalto(m,fc,Fs)

  wc = fc*2*pi/Fs;
  x1 = linspace(-m/2,m/2,m+1);
  x2 = linspace(0,m,m+1);

  h1 = sinc((pi-wc)*x1/pi)*(pi-wc)/pi;
  h2 = cos(pi*x1);
  w = 0.54-0.46*cos(2*pi*x2/m);

  h = h1.*h2.*w;
end

function c = conv_muestras(i,overlap)
   c=round((i-1)*overlap+1);
end
function [Np,ip,fp] = plantilla(inicio,fin,fin_ant,inicio_sig,a)
    ancho = round(fin-inicio);
    Np1 = round(a*ancho/2);
    Np2 = round(a*ancho/2);
    fp = fin-ancho/2+Np2-1;
    ip = inicio+ancho/2-Np1;
    if fp > inicio_sig #fin plantilla toma parte del sig click
      fp = inicio_sig-1;
      Np2 = ancho/2+inicio_sig-fin;
    endif
    if ip < fin_ant #inicio plantilla toma parte del ant click
        ip = fin_ant;
        Np1 = ancho/2+inicio-fin_ant;
    endif
    Np=Np1+Np2;
end

function ir =  posicion_reemplazo(y,xp,Np,Nw,ip,side)
      if side == 'l'
        N = -Nw;
      elseif side == 'r'
        N = Np;

      endif
      xw = y(ip+N+1:ip+N+Np);
      R_max = dot(xp,xw);
      ir = ip+N;

      for k = 2:Nw-Np+1
        x_w = y(ip+N+k+1:ip+N+k+Np);
        R = dot(xp,x_w);
        if R > R_max
         ir = ip+N+k;
         R_max = R;
        endif
      endfor
end

function vs = ventana_suave(Np)
  Np2 = round(Np/2);
  h = Np2 - 50;
  k = 0.5;
  x = linspace(-Np2,Np2,Np);
  vs = 1./ ((1+exp(-k*(x+h))).* (1+exp(k*(x-h))))' ;
end

# Lectura del audio
filename = 'audio_con_clicks.wav';
[y, Fs] = audioread(filename);
t = (0:length(y)-1)/Fs;
tf = (length(y)-1)/Fs;

## ANALISIS DE LA SEÑAL
# Figura audio wav
figure(1);
plot(t, y);
axis('tight');
xlabel('Tiempo (s)');
ylabel('Amplitud');
print('senial.png');

# Espectrograma
nfft = 1024;
ventana=hamming(nfft);
overlap = nfft/2;

figure(2);
specgram(y, nfft, Fs, ventana, overlap);
colormap jet;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
print('espectrograma.png');

##FILTRADO
m = 300;
fc = 4000;
h = pasaalto(m,fc,Fs);
H = fft(h);
figure(11);
zplane(h,1);
figure(3);
stem(linspace(0,length(h)-1,length(h)),h);
xlabel('Tiempo');
ylabel('Amplitud');
print('rta_impulso.png');

figure(4);
freqz(h);
print('rta_en_frecuencia.png');

y_filtrada = filter(h, 1, y);

[S_f,f_f,t_f] = specgram(y_filtrada, nfft, Fs, ventana, overlap);

figure(5);
imagesc(t_f, f_f, 20*log10(abs(S_f)));
axis xy;
#colorbar;
colormap jet;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
print('espectrograma_filtro.png');

# Calculo de enrgias
energias = sum(abs(S_f).^2);
energia_media = mean(energias);

figure(6); hold on;
xe = linspace(0-150/Fs,tf-150/Fs,length(t_f));
stem(xe,energias);
line ([0 tf], [energia_media energia_media]);
legend('Energia por intervalo E_i','Energia Media');
legend('Location', 'northwest');
axis('tight');
xlabel('Tiempo');
ylabel('Energia');
print('energias.png');

# RESTAURACION

# condiciones iniciales
if energias(1)>energia_media
  control=true;
  inicio_sig=1;
  reemplazo1=true;
else
  control=false;
  reemplazo1 = false;
endif

a = 1.3;
b = 1.5;
fin_ant = 1;
avance = 10;
cant_clicks=0;
for i = 1:length(t_f)
  if energias(i)>energia_media && control==false
    inicio_sig = conv_muestras(i,overlap)-m/2;
    control=true;
    if reemplazo1
      [Np,ip,fp] = plantilla(inicio,fin,fin_ant,inicio_sig,a);
      xp = y(ip:fp);
      fin_ant = fin;
      cant_clicks=cant_clicks+1;
      # Ventana izquierda de reemplazo
      Nw1 = round(b*Np);
      if ip >= Np #para primer click
        if ip < Nw1
          Nw1 = ip;
        endif
        ir_1 = posicion_reemplazo(y,xp,Np,Nw1,ip,'l');
        x_w1=y(ir_1:ir_1+Np-1);;
      else
        Nw1 = 0;
        x_w1 = zeros(Np,1);
      endif

      # Ventana derecha de reemplazo
      Nw2=round(b*Np);
      if inicio_sig-fp >= Np
        if inicio_sig-fp < Nw2
            Nw2=inicio_sig-fp; #achico ventana
        endif
        ir_2 = posicion_reemplazo(y,xp,Np,Nw2,ip,'r');
        x_w2=y(ir_2:ir_2+Np-1);
      else
          #Buscar en clicks siguientes
          inicio_sig2 = inicio_sig;
          Np_interv=ceil((Np-1)/overlap); #redondea arriba
          sig=false;
          j=i;
          control2=true;
          while ( sig==false && j<i+avance && j<=length(t_f))
            if energias(j)>energia_media && control2==false
              inicio_sig2 = conv_muestras(j,overlap)-m/2;
              control2=true;
            elseif energias(j)<=energia_media && control2==true
              fin2 = conv_muestras(j,overlap)-m/2;
              control2=false;
              inicio2 = inicio_sig2;
              mover_r=0;
            elseif energias(j)<=energia_media && control2==false
              mover_r=mover_r+1;
              if mover_r>=Np_interv
                Nw2=conv_muestras(mover_r,overlap)-m/2;
                ir_2 = posicion_reemplazo(y,xp,Np,Nw2,fin2-Np,'r');
                x_w2=y(ir_2:ir_2+Np-1);
                sig=true;
              endif
            endif
            j=j+1;
            if j>length(t_f) || j >= i+avance
              Nw2=0;
              x_w2=zeros(Np,1);
            endif
          endwhile
      endif

      ventana_empalme = ventana_suave(Np);
      if Nw1 != 0 && Nw2 != 0
        alpha = 0.5;
      else
        alpha = 1;
      endif
      ventana_reemplazo = xp.*(1-ventana_empalme)+alpha*(x_w1+x_w2).*ventana_empalme;
      y(ip:fp)=ventana_reemplazo;

    else
      reemplazo1=true;
    endif
  elseif energias(i)<=energia_media && control==true
    fin = conv_muestras(i,overlap)-m/2;
    control=false;
    inicio = inicio_sig;
  endif
endfor

cant_clicks

#final ultimo click
if energias(length(t_f))< energia_media
  [Np,ip,fp] = plantilla(inicio,fin,fin_ant,length(y),a);
  xp = y(ip:fp);
  fin_ant = fin;

  Nw1 = round(b*Np);
  ir_1 = posicion_reemplazo(y,xp,Np,Nw1,ip,'l');
  x_w1= y(ir_1:ir_1+Np-1);;
  Nw2 = round(b*Np);
  if length(y)-fp >= Np
    if length(y)-fp < Nw2
        Nw2=length(y)-fp; #achico ventana
    endif
    ir_2 = posicion_reemplazo(y,xp,Np,Nw2,ip,'r');
    x_w2=y(ir_2:ir_2+Np-1);
  else
      Nw2 = 0;
      x_w2 = zeros(Np,1);
  endif
  ventana_empalme = ventana_suave(Np);
  if Nw1 != 0 && Nw2 != 0
    alpha = 0.5;
  else
    alpha = 1;
  endif
  ventana_reemplazo = xp.*(1-ventana_empalme)+alpha*(x_w1+x_w2).*ventana_empalme;
  y(ip:fp)=ventana_reemplazo;
else #no hay fin del ultimo click
  Np = round(a*(length(y)-inicio_sig));
  Nw1 = round(b*Np);
  Nw2 = 0;
  xp = y(length(y)-Np+1:length(y));
  ir_1 = posicion_reemplazo(y,xp,Np,Nw1,length(y)-Np,'l');
  x_w1=y(ir_1:ir_1+Np-1);;
  x_w2 = zeros(Np,1);
  ventana_empalme = ventana_suave(Np);
  if Nw1 != 0 && Nw2 != 0
    alpha = 0.5;
  else
    alpha = 1;
  endif
  ventana_reemplazo = xp.*(1-ventana_empalme)+alpha*(x_w1+x_w2).*ventana_empalme;
  y(length(y)-Np+1:end)=ventana_reemplazo;
endif


figure(7);
specgram(y, nfft, Fs, ventana, overlap);
colormap jet;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
print('espectrograma_restaurado.png');

figure(8);
plot(t, y);
axis('tight');
xlabel('Tiempo (s)');
ylabel('Amplitud');
print('senial_restaurada.png');

filename_r = 'audio_mejorado.wav';
audiowrite(filename_r,y,Fs);

figure(9);
N = 512;
ventana_empalme = ventana_suave(N);
x3=linspace(1,N,N);
plot(x3,ventana_empalme);
axis('tight');
xlabel('Tiempo');
ylabel('Amplitud');
print('ventana_empalme.png');



