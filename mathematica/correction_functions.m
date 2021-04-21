(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
ClearAll;
SetOptions[Plot, Filling->Axis, GridLines->Automatic, Frame->True];
SetOptions[Graphics,Axes->False, Frame->False,LabelStyle->Directive[Black,Medium]];

RadauP[p_,x_]  :=(-1)^p*(LegendreP[p,x]-LegendreP[p-1,x])/2

psi[c_,p_]:=(c*(2*p+1)*(Factorial[2*p]/(2^p *Factorial[p]))^2)/2
Ges[c_,p_,x_]:=(((-1)^p)/2)*(LegendreP[p,x]-((psi[c,p])*LegendreP[p-1,x]+LegendreP[p+1,x])/(psi[c,p]+1)) (*Enegy Stable Correction Function*)


For [order=2,order<=5,order++,
	(* Define custum functions for the relevant polynomials *)
	Jac[x_]  :=JacobiP[order,a,b,x];
	Leg[x_]  :=ColorData[97,4]ColorData[97,4][order,x];
	Legp[x_] :=LegendreP[order-1,x];
	Rad[x_]  :=RadauP[order,x];
	Radp[x_] :=RadauP[order-1,x];

	(*Define Correction Function*)
	Gdg[x_]:=Rad[x];
	G2[x_]:=(order-1)/(2*order-1)*Rad[x]+(order)/(2*order-1)*Radp[x];

	pts={x,0}/.Solve[Legp[x]==0,x,Reals];
	pts=AppendTo[pts,{1,0}];
	pts=AppendTo[pts,{-1,1}];
	Gga[x_]:=InterpolatingPolynomial[pts,x];
	
	(*Print the current solution interpolation order*)
	Print["Order ", order, " solution => Order", order+1, " Flux"];

	(*Plot the Correction Functions*)
	dgGraph = Plot[Gdg[x], {x,-1,1}, PlotRange->{{-1.1,1.1},{-0.6,1.1}}, AspectRatio->Automatic, PlotStyle->ColorData[97,1], PlotLabel->"G_dg"];
	gaGraph = Plot[Gga[x], {x,-1,1}, PlotRange->{{-1.1,1.1},{-0.6,1.1}}, AspectRatio->Automatic, PlotStyle->ColorData[97,2], PlotLabel->"G_ga"];
	g2Graph = Plot[ G2[x], {x,-1,1}, PlotRange->{{-1.1,1.1},{-0.6,1.1}}, AspectRatio->Automatic, PlotStyle->ColorData[97,3], PlotLabel->"G_2"];
	GraphicsRow[{Show[dgGraph, gaGraph, g2Graph, PlotLabel->"All"], dgGraph, gaGraph, g2Graph}, ImageSize->Full, Spacings->0]//Print;

	(*Plot the Derivatives of the Correction Functions*)
	dgGraphd = Plot[Gdg'[x], {x,-1,1}, PlotRange->{{-1.1,1.1}, Automatic}, AspectRatio->GoldenRatio, PlotStyle->ColorData[97,1], PlotLabel->"G_dg'"];
	gaGraphd = Plot[Gga'[x], {x,-1,1}, PlotRange->{{-1.1,1.1}, Automatic}, AspectRatio->GoldenRatio, PlotStyle->ColorData[97,2], PlotLabel->"G_ga'"];
	g2Graphd = Plot[ G2'[x], {x,-1,1}, PlotRange->{{-1.1,1.1}, Automatic}, AspectRatio->GoldenRatio, PlotStyle->ColorData[97,3], PlotLabel->"G_2'"];
	GraphicsRow[{Show[dgGraphd, gaGraphd, g2Graphd, PlotLabel->"All"], dgGraphd, gaGraphd, g2Graphd}, ImageSize->Full, Spacings->0]//Print;
	
	(*Print Analytical Form of the Correction Functions*)
	GraphicsRow[Gdg[x], Gga[x], G2[x]] //Print ;
	Print[Gdg[x]];
	Print[Gga[x]];
	Print[G2[x]];

	(*Calculate the Correction Functions and Derivatives at Reference Points*)
	Print["x =", PaddedForm[N[x/.Solve[Leg[x]==0,x,Reals]], {8,6}]];
	Print["x =", PaddedForm[N[Range[-1,1,2/(order+2)]], {8,6}]];
	Print["   G_dg =", PaddedForm[ N[ Map[ Gdg , Range[-1,1,2/(order+2)] ] ], {8,6}] ];
	Print["   G_dg'=", PaddedForm[ N[ Map[ Gdg', Range[-1,1,2/(order+2)] ] ], {8,6}] ];
	Print["   G_ga =", PaddedForm[ N[ Map[ Gga , Range[-1,1,2/(order+2)] ] ], {8,6}] ];
	Print["   G_ga'=", PaddedForm[ N[ Map[ Gga', Range[-1,1,2/(order+2)] ] ], {8,6}] ];
	Print["   G_2  =", PaddedForm[ N[ Map[ G2  , Range[-1,1,2/(order+2)] ] ], {8,6}] ];
	Print["   G_2' =", PaddedForm[ N[ Map[ G2' , Range[-1,1,2/(order+2)] ] ], {8,6}] ];
	Print[];
	
	(*Get the Roots of the Polynomials*)
	Print["G_dg  roots: ", PaddedForm[ N[ x/.Solve[Gdg [x]==0, x, Reals] ], {8,6}] ];
	Print["G_dg' roots: ", PaddedForm[ N[ x/.Solve[Gdg'[x]==0, x, Reals] ], {8,6}] ];
	Print["G_ga  roots: ", PaddedForm[ N[ x/.Solve[Gga [x]==0, x, Reals] ], {8,6}] ];
	Print["G_ga' roots: ", PaddedForm[ N[ x/.Solve[Gga'[x]==0, x, Reals] ], {8,6}] ];
	Print["G_2   roots: ", PaddedForm[ N[ x/.Solve[G2  [x]==0, x, Reals] ], {8,6}] ];
	Print["G_2'  roots: ", PaddedForm[ N[ x/.Solve[G2' [x]==0, x, Reals] ], {8,6}] ];
]



