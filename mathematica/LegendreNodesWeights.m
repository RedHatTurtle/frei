(* ::Package:: *)

(* ::Input:: *)
(*maxOrder=20*)


(* ::Input:: *)
(*TableForm[x/.Table[NSolve[LegendreP[n,x]==0, x, 37],{n,1,maxOrder}]]*)


(* ::Input:: *)
(*Export["./legendreAbscissa.txt",%]*)


(* ::Input:: *)
(*TableForm[   Table[  ((2*(1-x^2)/(n+1)^2)/LegendreP[n+1,x]^2)/. NSolve[LegendreP[n,x]==0, x, 37],{n,1,maxOrder}   ]   ]*)


(* ::Input:: *)
(*Export["./legendreWeights.txt",%]*)


(* ::Input:: *)
(**)
(**)
(**)
(**)
(**)


(* ::Input:: *)
(*TableForm[Table[Append[Prepend[x/. NSolve[D[LegendreP[n-1,x],x]==0, x, 37],-1] ,1],{n,3,maxOrder}]]*)


(* ::Input:: *)
(*Export["./DlegendreAbscissa.txt",%]*)


(* ::Input:: *)
(*TableForm[   Table[  (2/(n^2-n)/LegendreP[n-1,x]^2)/. Append[Prepend[ NSolve[D[LegendreP[n-1,x],x]==0, x, 37],{x->-1} ],{x->1}],{n,3,maxOrder}   ]   ]*)


(* ::Input:: *)
(*Export["./DlegendreWeights.txt",%]*)


(* ::Input:: *)
(*"DlegendreWeights.txt"*)
