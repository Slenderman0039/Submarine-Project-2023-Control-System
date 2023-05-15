 close all; clc
%definizione variabili
kteta = 143;
kr = 168;
b = 1.2;
teta = pi/6;
B = 20;
F = -kteta/(2*b);
FL = -sqrt(3)/2 * kteta/b - 2*B;
kzeta = 133;
M = 580;
J = 560;

%matrici forma di stato normalizzato del modello linearizzato approssimato
A = [ 0 0 1 0; 0 0 0 1; 0 -(F*sin(teta)+FL*cos(teta))/M -(kzeta*cos(teta))/M 0; 0 -(kteta*cos(teta))/J 0 -kr/J]
B = [ 0; 0; cos(teta)/M; -b/J]
C = [ 1 -b*cos(teta) 0 0 ]
D = 0

% autovettori ed autovalori di A
% T = autovettori
% landa = autovalori
[T,landa] = eig(A)
% colonna i-esima della matrice T è l'autovettore i-esimo dell'autovalore
% i-esimo

% Sistema 
%sys = ss(A, B, C, D);
% Funzione di trasferimento
%G = tf(sys);

% Poli
% P = pole(G)

% Zeri
% Z = zero(G)

%calcolo scomposizione del numeratore
alfa = cos(teta) + b*b*cos(teta)*M/J
beta = (cos(teta)*kr/J + b*b*cos(teta)*cos(teta)*kzeta/J )/alfa
gamma = ((F*sin(teta)+FL*cos(teta))*b/J + kteta*cos(teta)*cos(teta)/J )/alfa

delta = beta^2-4*gamma
s1 = (-beta + sqrt(delta))/2
s2 = (-beta - sqrt(delta))/2

costanteBode = alfa/M

%Sistema linearizzato
sys = ss(A, B, C, D)
G = tf(sys)
[z, p, k ] = zpkdata(sys)
G = zpk(G)
%bode(G)
%Funzione di trasferimento in anello chiuso Gc(s) approssimato
%dopo aver esportato Ctrl, oppure lo stesso se si esporta Gc con r2y
%Funzione di trasferimento in forma zpk
%n è il numero di ordini del sistema
% il vettore [wn,zeta,p] resistuisce le pulsazioni, lo smorzamento ed i
% poli del sistema approssimato
Gcs=zpk(Gc)
n=order(Gcs)
[wn,zeta,p] = damp(Gcs)
% bode(Gc) per il plot dei diagrammi di bode della Gc(s)
% bode(C*G) per il plot dei diagrammi di bode di C(jw)*G(jw)
% stepinfo(Gc) informazioni relative alla funzione step plottata con
% step(Gc)