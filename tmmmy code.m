clear
clc
degrees = pi/180; 
theta = 17 * degrees;      % elevation angle 
phi = 12 * degrees;        % azimuthal angle 
lamda = (300:1:700) * 1e-9;    % free space wavelength 
k0 = (2 * pi) ./ lamda;
ni = 1.0;                    % incident medium refractive index
Kx = k0 * ni * sin(theta) * cos(phi);
nh = 5;
nl = 3;
N = 10;
dif = nh - nl;
aver = (nh + nl) / 2;
r = dif / (2 * aver);
xh = 0.5 * (10^-6);
xl = 0.25E-6;
alphah = Kx * xh;
alphal = Kx * xl;
alpha = alphah + alphal;
delalpha = alphah - alphal;
G = (-2 * r * sin(alphah) * 1i) ./ (1 - r^2);
cg = conj(G);
F = (exp(alpha * 1i) - r * r .* exp(delalpha * 1i)) ./ (1 - r^2);
fr = real(F);
fi = imag(F);
apsi = acosh(-fr);

% Preallocate R and T arrays for efficiency
R = zeros(1, length(lamda));
T = zeros(1, length(lamda));

% Loop through each wavelength
for idx = 1:length(lamda)
    % Calculating the transfer matrix elements for each wavelength
    M = [0 0; 0 0]; % Initialize a 2x2 matrix
    
    M(1,1) = cosh(N * apsi(idx)) - sinh(N * apsi(idx)) * (fi(idx) / sqrt((fr(idx)^2) - 1)) * 1i;
    M(1,2) = -sinh(N * apsi(idx)) * (cg(idx) / sqrt((fr(idx)^2) - 1));
    M(2,1) = -sinh(N * apsi(idx)) * (G(idx) / sqrt((fr(idx)^2) - 1));
    M(2,2) = cosh(N * apsi(idx)) + sinh(N * apsi(idx)) * (fi(idx) / sqrt((fr(idx)^2) - 1)) * 1i; 

    % Calculating reflectance and transmittance
    R(idx) = abs(M(2,1) / M(1,1));
    T(idx) = abs(1 / M(1,1));
end

% Plot the results
plot(lamda, R, 'r', lamda, T, 'b');
xlabel('Incident Wavelength (m)');
ylabel('Reflectance or Transmittance');
legend('Reflectance', 'Transmittance');
