(* the loop *)
iterativesolver[innercut_, outercut_, innertol_, integraltol_, 
  grid_, \[Phi]_, \[Nu]_, \[CapitalLambda]_, initialinputh_, 
  initialh0_] := Module[
  {integralerr, innererr, i, j
   , data, plots = {}
   , inputh = initialinputh, hfromP = initialinputh, h0 = initialh0
   , h0new
   , inputhovergrid, hfromPovergrid
   , P, Pressure, F
   , uz, uzovergrid, delta},
  
  integralerr = 2 integraltol; i = 0; plots = {{}};
  While[integralerr > integraltol,
   i++; j = 0; innererr = 2 innertol;
   While[innererr > innertol,
    (* increment iteration counter*)
    j++;
    Print["i ", i, " j ", j];
    
    (* update the current h*)
    
    inputhovergrid = (0.5 inputh[#] + 0.5 hfromP[#]) & /@ grid;
    inputh = 
     Interpolation[MapThread[{#1, #2} &, {grid, inputhovergrid}]];
    
    (* ok go! get the pressure*)
    
    Pressure = 
     NDSolve[{1/r D[r   (inputh[r]^3)/12 P'[r], r] == 
        1/\[Phi]^2 (-1)/inputh[r], P[outercut] == 0, 
       P'[innercut] == 0}, P, {r, innercut, outercut}];
    Pressure = (P /. Pressure)[[1]];
    
    (* from the pressure, get the elastic displacement*)
    
    uz[r_] := 
     4 (1 - \[Nu]^2)/\[Pi]  NIntegrate[
       t/(t + r) Pressure[t] EllipticK[ (4 t r)/(t + r)^2], {t, 
        innercut, outercut}, 
       Method -> {Automatic, "SymbolicProcessing" -> 0}];
    (* construct an interpolating function over out grid*)
    
    uzovergrid = uz@grid;
    
    (* the new height function*)
    
    hfromPovergrid = h0 + grid^2/2 + uzovergrid;
    delta = Interpolation[MapThread[{#1, #2} &, {grid, uzovergrid}]];
    hfromP = 
     Interpolation[MapThread[{#1, #2} &, {grid, hfromPovergrid}]];
    
    (* have we converged at the origin? L2 norm*)
    (*nom=Sqrt[
    NIntegrate[(hfromP[r]-inputh[r])^2,{r,innercut,outercut},
    Method\[Rule]{Automatic,"SymbolicProcessing"\[Rule]0}]];
    denom=Sqrt[NIntegrate[inputh[r]^2,{r,innercut,outercut},
    Method\[Rule]{Automatic,"SymbolicProcessing"\[Rule]0}]];
    innererr=nom/denom;*)
    
    innererr = 
     Abs[((hfromP[innercut] - inputh[innercut])/inputh[innercut])];
    Print["inner error", innererr];
    
    AppendTo[
     plots[[i]], {Pressure[r], inputh[r], hfromP[r], delta[r]}];
    ];
   AppendTo[plots, {}];
   
   (* with the converged h and pressure, 
   do the integral to get the overall load*)
   
   F = 2*\[Pi]*NIntegrate[r Pressure[r], {r, innercut, outercut}];
   (* have we converged?*)
   integralerr = Abs[(F - 1)];
   
   If[integralerr < integraltol,
    (* append to lists*)
    
    data = {\[Phi], h0, inputh, Pressure, delta};
    ];
   
   (* new h0 guess, 
   depending on whether the force is too big or too small*)
   
   h0new = (h0) (F^\[CapitalLambda]);
   
   (* for the new profile, 
   which is the initial guess for the inner loop, use the old shape, 
   but just shifted down with the new h0*)
   
   inputh = 
    FunctionInterpolation[
     inputh[x] + (h0new - h0), {x, innercut, outercut}];
   hfromP = inputh;
   h0 = h0new;
   
   Print[" integralrr: ", integralerr];
   Print[" previous h0: ", h0];
   Print[" F from previous h0: ", F];
   Print["origin height:", hfromP[innercut]];
   Print[" new h0 ", h0new];
   ];
  
  (* at the end, return the gathered data and the intermediate plots*)

    Return[{data, plots}]
  ]