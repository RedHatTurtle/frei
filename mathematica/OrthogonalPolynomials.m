(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
ClearAll;
SetOptions[Plot, Filling->Axis, GridLines->Automatic, Frame->True];


a = -1/2;
b = -1/2;


For [order=1,order<=9,order++,
	(* Define custum functions for the relevant polynomials *)
	Jac[x_]  :=JacobiP[order,a,b,x];
	Leg[x_]  :=LegendreP[order,x];
	Legp[x_] :=LegendreP[order-1,x];
	Rad[x_]  :=(-1)^order*(Leg[x]-Legp[x])/2;

	(*Plot the Orthogonal Polynomials*)
	jGraph = Plot[Jac[x], {x,-1,1}, PlotRange->{{-1.1,1.1},{-1.1,1.1}}, AspectRatio->1, PlotStyle->ColorData[97,1], PlotLabel->"Jacobi"];
	lGraph = Plot[Leg[x], {x,-1,1}, PlotRange->{{-1.1,1.1},{-1.1,1.1}}, AspectRatio->1, PlotStyle->ColorData[97,4], PlotLabel->"Legendre"];
	rGraph = Plot[Rad[x], {x,-1,1}, PlotRange->{{-1.1,1.1},{-1.1,1.1}}, AspectRatio->1, PlotStyle->ColorData[97,3], PlotLabel->"Radau"];
	GraphicsRow[{Show[jGraph,lGraph,rGraph,PlotLabel->"All"],jGraph,lGraph,rGraph},ImageSize->1150, Spacings->0]//Print;
	
	(*Print Analytical Form of the Polynomials*)
	Print[Jac[x]];
	Print[Leg[x]];
	Print[Rad[x]];

	(*Calculate the Polynomials at Reference Points*)
	Print["x =", PaddedForm[N[Range[-1,1,2/(order+2)]], {8,6}]];
	Print["   Jac  =", PaddedForm[ N[ Map[ Jac  , Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Jac' =", PaddedForm[ N[ Map[ Jac' , Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Jac''=", PaddedForm[ N[ Map[ Jac'', Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Leg  =", PaddedForm[ N[ Map[ Leg  , Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Leg' =", PaddedForm[ N[ Map[ Leg' , Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Leg''=", PaddedForm[ N[ Map[ Leg'', Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Rad  =", PaddedForm[ N[ Map[ Rad  , Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Rad' =", PaddedForm[ N[ Map[ Rad' , Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print["   Rad''=", PaddedForm[ N[ Map[ Rad'', Range[-1,1,2/(order+2)] ] ], {9,6}] ];
	Print[];
	
	(*Get the Roots of the Polynomials*)
	Print["Jac  roots: ", PaddedForm[ N[ x/.Solve[Jac [x]==0, x, Reals] ], {8,6}] ];
	Print["Jac' roots: ", PaddedForm[ N[ x/.Solve[Jac'[x]==0, x, Reals] ], {8,6}] ];
	Print["Leg  roots: ", PaddedForm[ N[ x/.Solve[Leg [x]==0, x, Reals] ], {8,6}] ];
	Print["Leg' roots: ", PaddedForm[ N[ x/.Solve[Leg'[x]==0, x, Reals] ], {8,6}] ];
	Print["Rad  roots: ", PaddedForm[ N[ x/.Solve[Rad [x]==0, x, Reals] ], {8,6}] ];
	Print["Rad' roots: ", PaddedForm[ N[ x/.Solve[Rad'[x]==0, x, Reals] ], {8,6}] ];
]



