theta = 0:0.01:2*pi();

for kk = 1 :length(theta)
    eraw = 40385;
    n = 1.00277387;
    ideg = 14.2407;
    wdeg = 318.6668;
    Mdeg = 85.9519;
    t = 264.51782528;
    t0 = 0;
    omegadeg = 225.5871;

    mu = 3.98604419e14;

    e = eraw*1e-07;
    a = mu^(1/3) / (2*n*pi()/ 86400)^(2/3);
    i = ideg*(pi()/180);
    w = wdeg*(pi()/180);
    omega = omegadeg*(pi()/180);
    Mt = Mdeg*(pi()/180);
    E = theta(kk);

    Vt = 2 * atan2( sqrt( 1 + e ) * sin( E / 2 ) , sqrt(1 - e ) * cos( E / 2));

    rc = a * (1 - e * cos(E));

    vcoeff = sqrt(mu * a ) / rc;

    Ot = [rc * cos(Vt), rc * sin(Vt),0];
    Odott = [vcoeff * -sin(E), vcoeff * sqrt(1 - e^2) * cos(E),0];
    rx(kk) = (Ot(1) * (cos(w) * cos(omega) - sin(w) * cos(i) * sin(omega)) - Ot(2) * (sin(w) * cos(omega) + cos(w) * cos(i) * sin(omega)));
    ry(kk) = (Ot(1)* (cos(w) * sin(omega) + sin(w) * cos(i) * cos(omega)) + Ot(2) * (cos(w) * cos(i) * cos(omega) - sin(w) * sin(omega)));
    rz(kk) = (Ot(1) * (sin(w) * sin(i)) + Ot(2) * (cos(w) * sin(i)));

    rdotx(kk) = (Odott(1) * (cos(w) * cos(omega) - sin(w) * cos(i) * sin(omega)) - Odott(2) * (sin(w) * cos(omega) + cos(w) * cos(i) * sin(omega)));
    rdoty(kk) = (Odott(1) * (cos(w) * sin(omega) + sin(w) * cos(i) * cos(omega)) + Odott(2) * (cos(w) * cos(i) * cos(omega) - sin(w) * sin(omega)));
    rdotz(kk) = (Odott(1) * (sin(w) * sin(i)) + Odott(2) * (cos(w) * sin(i)));
end
scatter3(rx,ry,rz)

%scatter3(r_x,r_y,r_z,70)
%legend('Eccentricity Variant', 'Time Variant')