close all
clear

% Metode for ?? forutse svingninger i derivatmarkedet ved en
% fourierutvikling. Tidssteg b??r v??re definert basert p?? grafen vi bruker
% som input, f.eks hvert steg er 1 uke eller 1 min. Dette b??r brukes i
% output.

% Input plasseres her, beskrivelse underveis
antall_variabler = 10000; % Antall ledd vi bruker i fourierutviklingen
ekstra_tidsledd = 20; % Antall tidssteg vi ??nsker ?? se inn i fremtiden
ymax = 1000; % Toppunktet p?? grafen for reskalering

% 1)
% F??rst m?? en matrise defineres. Matrisen m?? beskrive et bilde av
% tidsutviklingen til kursen in question.

% Loader bildet
%Kurs = imread('kursutvikling.jpeg');
%Kurs = imread('110613GOLDUSD.jpeg');
%Kurs = imread('120613SILVERUSD.jpeg');
Kurs = imread('us2907.jpg');

% Visualiserer matrisen
%{
figure(1)
Kurs = Kurs(:,:,1);
imagesc(Kurs)
title('R??data')
%}

% Definerer st??rrelsen p?? matrisen
lx = size(Kurs,1);
ly = size(Kurs,2);

% Definerer vektoren
Y = zeros(ly,1);

% Reduserer matrisen til en vektor
for j=1:ly
    min = 100000000;
    for i=1:lx
        if (Kurs(i,j) < min)
            min = Kurs(i,j);
            pos = i;
        end
    end
    Y(j) = lx - pos;
end

% Visualisere input til tiln??rmern

figure(1)
ymaxtemp = max(Y);
plot(Y/max(Y)*ymax)
title('Input til tiln??rmern')


% 2)
% Y inneholder n?? relativ funksjonsverdi til grafen. N?? m?? vi tiln??rme
% dette med en evig rekke cos og sin slik Fourier beskriver

% Fourier sier:
% Y = a0 + a1*cos(x) + b1*sin(x)
%        + a2*cos(2x) + b2*sin(2x)
%        + a3*cos(3x) + b3*sin(3x)
%        + ...

% Vi begynner med en enkel antagelse at a0 = Y(1) = startpunktet
% Vi ??nsker s?? ?? tiln??rme med s?? mange variabler at vi finner en eksakt
% beskrivelse av funksjonen med minste kvadraters metode. Her kan vi
% inkludere en optimizer for ?? finne hvilke a_n og b_n vi ??nsker, med
% n=1:antall_variabler.

% Definerer variablene v??re
a = zeros(antall_variabler, 1);
b = zeros(antall_variabler, 1);

% Definerer a0
a0 = Y(1);

% Inkluderer a0
for i=1:ly
    Y(i) = Y(i) - a0;
end

varianter = 100;
% For antallet variabler vi kj??rer en optimering
for Z = 1:antall_variabler
    minimum = 10000000000000000000000000000000000;

    % Vi bruker minste kvadraters metode forel??pig
    % F??rst for a - leddene
    
    % F??rst cosinus leddet
    for V = 0:varianter
        % Definerer en midlertidig sum og ??nsker ?? minimere denne
        diff_sum1 = 0;
        temp_a = (varianter - 2 * V)/varianter * max(abs(Y));
        
        temp_sum = 0;
        for i=1:ly
            temp_sum = temp_sum + (Y(i) - temp_a*cos(Z * i / lx))^2;
        end
        
        if (temp_sum < minimum)
            minimum = temp_sum;
            a(Z) = temp_a;
        end
    end
    
    for (i=1:ly)
        Y(i) = Y(i)-a(Z)*cos(Z*i/lx);
    end
    
    % S?? sinus leddet
    for V = 0:varianter
        % Definerer en midlertidig sum og ??nsker ?? minimere denne
        diff_sum1 = 0;
        temp_a = (varianter - 2 * V)/varianter * max(abs(Y));
        
        temp_sum = 0;
        for i=1:ly
            temp_sum = temp_sum + (Y(i) - temp_a*sin(Z * i / lx))^2;
        end
        
        if (temp_sum < minimum)
            minimum = temp_sum;
            b(Z) = temp_a;
        end
    end
    
    for (i=1:ly)
        Y(i) = Y(i)-b(Z)*sin(Z*i/lx);
    end
    
end

figure(2)
plot(Y*ymax/ymaxtemp)
title('Resten. Det vi ikke klarte ?? tiln??rme')

%{
% Output, leddene v??re i fouriertiln??rmern
a0 = a0
a
b
%}

% Dette outputet brukes til ?? konstruere funksjonen v??r
nye_ledd = ekstra_tidsledd;
y = zeros(ly + nye_ledd,1);

for (i=1:ly+nye_ledd)
    y(i) = y(i) + a0;
    for (j=1:antall_variabler)
        y(i) = y(i) + a(j)*cos(j*i/lx);
        y(i) = y(i) + b(j)*sin(j*i/lx);
    end
end

y = y/(max(y));
y = y*ymax;

% Vi plotter tiln??rmingen v??r med leddene i fremtiden
xlim = [ly ly];
ylim = [0 ymax];
figure(3)
plot(y)
title('Tiln??rmingen v??r, med ekstra ledd')
hold on
plot(xlim, ylim, 'r')
hold off

FREMTID = zeros(nye_ledd, 1);
for i=1:ly+nye_ledd
    if (i > ly)
        FREMTID(i-ly) = y(i);
    end
end

% Vi plotter fremtidsutviklingen v??r
figure(4)
plot(FREMTID)
title('Fremtiden')












