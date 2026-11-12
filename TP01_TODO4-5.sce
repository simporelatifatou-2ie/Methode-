function q = trapezes(t, y)
    n = length(t) - 1;
    if n < 1 | length(y) <> length(t) then
        error("Deux vecteurs de meme longueur sont requis.");
    end
    q = 0;
    for i = 1:n
        h = t(i+1) - t(i);
        if h <= 0 then error("Les abscisses doivent croitre."); end
        q = q + h * (y(i) + y(i+1)) / 2;
    end
endfunction

function q = simpson(t, y)
    n = length(t) - 1;
    if n < 2 | length(y) <> length(t) then
        error("Vecteurs incompatibles ou trop courts.");
    end
    if modulo(n, 2) <> 0 then error("Simpson exige n pair."); end
    h = t(2) - t(1);
    for i = 1:n
        pas = t(i+1) - t(i);
        if pas <= 0 then error("Les abscisses doivent croitre."); end
        if abs(pas-h) > 1e-10 * max(1, abs(h)) then
            error("Simpson exige des pas egaux.");
        end
    end
    q = y(1) + y(n+1);
    for i = 2:n
        if modulo(i, 2) == 0 then
            q = q + 4*y(i);
        else
            q = q + 2*y(i);
        end
    end
    q = h*q/3;
endfunction

function [t, y] = echantillonner(f, a, b, n)
    if n < 1 | b <= a then error("Intervalle ou maillage incorrect."); end
    t = zeros(1, n+1);
    y = zeros(1, n+1);
    for i = 1:n+1
        t(i) = a + (i-1)*(b-a)/n;
        y(i) = f(t(i));
    end
endfunction

function p = ordre_observe(E_n, E_2n)
    if E_n <= 0 | E_2n <= 0 then
        p = %nan;
    else
        p = log(E_n/E_2n)/log(2);
    end
endfunction

// ---------- Application ----------

// Cas continu : modele normalise d'une puissance variable.
function y = modele(x)
    y = exp(x);
endfunction
Iref = %e-1; // Integrale exacte sur [0,1].
ns = [4 8 16 32];
ET = zeros(1, length(ns)); ES = ET;
errT = ET; errS = ET;
for k = 1:length(ns)
    n = ns(k);
    [t, y] = echantillonner(modele, 0, 1, n);
    // TODO 1 : utiliser les deux quadratures du CM1.
    ET(k) = %nan;
    ES(k) = %nan;
    // TODO 2 : erreurs absolues par rapport a Iref.
    errT(k) = %nan;
    errS(k) = %nan;
    mprintf("n=%d T=%.9f S=%.9f ET=%.3e ES=%.3e\n", ..
            n, ET(k), ES(k), errT(k), errS(k));
end
for k = 1:length(ns)-1
    // TODO 3 : mesurer les ordres entre n et 2n.
    pT = %nan;
    pS = %nan;
    mprintf("n=%d -> %d : pT=%.3f pS=%.3f\n", ..
            ns(k), ns(k+1), pT, pS);
end

// Cas discret : releves de puissance d'une pompe, temps en minutes.
tmin = [0 10 20 30 40];
Pkw = [2.0 2.6 3.0 2.5 1.9];
th = zeros(1, length(tmin));
for i = 1:length(tmin)
    // TODO 4 : convertir les minutes en heures.
    th(i) = tmin(i) / 60;
end
// TODO 5 : energie par les deux quadratures (en kWh).
Etrap = trapezes(th, Pkw);
Esimp = simpson(th, Pkw);
mprintf("Pompe : trapezes=%.6f kWh ; Simpson=%.6f kWh\n", ..
        Etrap, Esimp);
mprintf("Ecart=%.6f kWh ; borne des lectures=%.6f kWh\n", ..
        abs(Etrap-Esimp), 0.1*40/60);
