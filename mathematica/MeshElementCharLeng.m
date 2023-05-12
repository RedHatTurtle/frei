(* ::Package:: *)

(* ::Title:: *)
(*Calculate the Min Height of a Mesh Element*)


ClearAll;
$MaxExtraPrecision=128;
ProjectPointToPlane[x_,{p0_,p1_,p2_}]:=
	Module[{y,nrm, n0},
		nrm=Normalize[Cross[p1-p0,p2-p0]];
		n0=Projection[p0,nrm];
		y=x-Projection[x,nrm]+n0
	]
ProjectPointToPlaneIntersection[x_,{p0_,p1_,p2_},{q0_,q1_,q2_}]:=
	Module[{y,line, vect,n0,v1,v2,v3},
	line={v1,v2,v3}/.Solve[{v1,v2,v3}\[Element] InfinitePlane[{p0,p1,p2}] &&
{v1,v2,v3}\[Element] InfinitePlane[{q0,q1,q2}]][[1]];
	vect=Normalize[(line/.v1->1)-(line/.v1->0)];
	y=Projection[x,vect]+(line/.v1->0)
	]


(* ::Subtitle:: *)
(*Define the vertices and enforce coplanarity*)


(*Fixed vertices*)
vert1={(0*237+9)/237,(0*727+8)/727,(0*461+1)/461};
vert2={(1*127+6)/127,(0*659+1)/659,(0*353-7)/353};
vert4={(0*337+8)/337,(2*911-6)/911,(0*241-4)/241};
vert5={(0*593-6)/593,(0*167-7)/167,(3*641+2)/641};

(*Adjustable vertices*)
vert3={(1*941+6)/941,(2*281-6)/281,(0*503+7)/503};
vert6={(1*281-8)/281,(0*821+9)/821,(3*379+2)/379};
vert8={(0*701+3)/701,(2*419-8)/419,(3*907+9)/907};
vert7={(1*827+1)/827,(2*467-2)/467,(3*733-4)/733};

(*Coplanarity enforcement*)
vert3=ProjectPointToPlane[vert3,{vert1,vert2,vert4}];
vert6=ProjectPointToPlane[vert6,{vert1,vert2,vert5}];
vert8=ProjectPointToPlaneIntersection[vert8,{vert1,vert4,vert5},{vert2,vert4,vert6}];
vert7={x,y,z}/.NSolve[{x,y,z}\[Element] InfinitePlane[{vert5,vert6,vert8}] &&{x,y,z}\[Element] InfinitePlane[{vert2,vert3,vert6}] &&{x,y,z}\[Element] InfinitePlane[{vert3,vert4,vert8}],WorkingPrecision->64][[1]];

(*Check the Prism inclined face and Hexa faces*)
If[And[ CoplanarPoints[{vert2,vert4,vert8,vert6}],
	CoplanarPoints[{vert1,vert4,vert3,vert2}], CoplanarPoints[{vert1,vert2,vert6,vert5}],
	CoplanarPoints[{vert1,vert5,vert8,vert4}], CoplanarPoints[{vert5,vert6,vert7,vert8}],
	CoplanarPoints[{vert4,vert8,vert7,vert3}], CoplanarPoints[{vert2,vert3,vert7,vert6}]
], Print["All nodes on quadrilateral faces are coplanar"], Print["Something is wrong"]];


(* ::Subtitle:: *)
(*Print the final verified vertices*)


Print[TextString[ScientificForm[N[vert1,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert2,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert3,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert4,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert5,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert6,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert7,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]
Print[TextString[ScientificForm[N[vert8,50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}],ListFormat->{"[ ",",\n  "," ]"}]]


(* ::Subtitle:: *)
(*Define infinite line/planes for each face*)


lineY0=InfiniteLine[{vert1,vert2}]; (*Quad Face 1*)
lineX1=InfiniteLine[{vert2,vert3}]; (*Quad Face 2*)
lineY1=InfiniteLine[{vert3,vert4}]; (*Quad Face 3*)
lineX0=InfiniteLine[{vert1,vert4}]; (*Quad Face 4*)
lineXY1=InfiniteLine[{vert2,vert4}]; (*Tria Face 2*)

(*Tetrahedra inclined plane*)
planeXYZ1=InfinitePlane[{vert2,vert4,vert5}];

(*Pyramid inclined planes*)
planeXZ1=InfinitePlane[{vert2,vert3,vert5}];
planeYZ1=InfinitePlane[{vert3,vert4,vert5}];

(*Prism inclined plane*)
If[CoplanarPoints[{vert2,vert4,vert6,vert8}],planeXY1=InfinitePlane[{vert2,vert4,vert6}],Print["Error xy1"]];

(*Cartesian/Hexahedra planes*)
If[CoplanarPoints[{vert1,vert4,vert5,vert8}],planeX0=InfinitePlane[{vert1,vert4,vert5}],Print["Error x0"]];
If[CoplanarPoints[{vert2,vert3,vert6,vert7}],planeX1=InfinitePlane[{vert2,vert3,vert6}],Print["Error x1"]];
If[CoplanarPoints[{vert1,vert2,vert3,vert4}],planeY0=InfinitePlane[{vert1,vert2,vert5}],Print["Error y0"]];
If[CoplanarPoints[{vert1,vert2,vert5,vert6}],planeY1=InfinitePlane[{vert4,vert3,vert8}],Print["Error y1"]];
If[CoplanarPoints[{vert4,vert3,vert7,vert8}],planeZ0=InfinitePlane[{vert1,vert2,vert4}],Print["Error z0"]];
If[CoplanarPoints[{vert5,vert6,vert7,vert8}],planeZ1=InfinitePlane[{vert5,vert6,vert8}],Print["Error z1"]];


(* ::Subtitle:: *)
(*Calculate the distance form nodes to opposing planes*)


Print["Line: ", TextString[ScientificForm[N[Min[
	EuclideanDistance[vert1,vert2]
],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Tria: ", TextString[ScientificForm[N[Min[
	RegionDistance[lineY0,vert4],
	RegionDistance[lineXY1,vert1],
	RegionDistance[lineX0,vert2]
],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Quad: ", TextString[ScientificForm[N[Min[
	Max[	RegionDistance[lineY0,vert3],	RegionDistance[lineY0,vert4]	],
	Max[	RegionDistance[lineX1,vert4],	RegionDistance[lineX1,vert1]	],
	Max[	RegionDistance[lineX1,vert1],	RegionDistance[lineY1,vert2]	],
	Max[	RegionDistance[lineX0,vert2],	RegionDistance[lineX0,vert3]	]
],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Tetr: ", TextString[ScientificForm[N[Min[
	RegionDistance[planeZ0,vert5],
	RegionDistance[planeY0,vert3],
	RegionDistance[planeXYZ1,vert1],
	RegionDistance[planeX0,vert2]
],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pyra: ", TextString[ScientificForm[N[Min[
	RegionDistance[planeZ0,vert5],
	Max[	RegionDistance[planeY0,vert3],	RegionDistance[planeY0,vert4]	],
	Max[	RegionDistance[planeXZ1,vert1],	RegionDistance[planeXZ1,vert4]	],
	Max[	RegionDistance[planeYZ1,vert1],	RegionDistance[planeYZ1,vert2]	],
	Max[	RegionDistance[planeX0,vert2],	RegionDistance[planeX0,vert3]	]
],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pris: ", TextString[ScientificForm[N[Min[
	Max[	RegionDistance[planeY0,vert4],	RegionDistance[planeY0,vert8]	],
	Max[	RegionDistance[planeXY1,vert1],	RegionDistance[planeXY1,vert5]	],
	Max[	RegionDistance[planeX0,vert2],	RegionDistance[planeX0,vert6]	],
	Max[	RegionDistance[planeZ0,vert5],	RegionDistance[planeZ0,vert6],	RegionDistance[planeZ0,vert8]	],
	Max[	RegionDistance[planeZ1,vert1],	RegionDistance[planeZ1,vert4],	RegionDistance[planeZ1,vert2]	]
],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Hexa: ", TextString[ScientificForm[N[Min[
	Max[	RegionDistance[planeZ0,vert5],	RegionDistance[planeZ0,vert6],	RegionDistance[planeZ0,vert7],	RegionDistance[planeZ0,vert8]	],
	Max[	RegionDistance[planeY0,vert4],	RegionDistance[planeY0,vert8],	RegionDistance[planeY0,vert7],	RegionDistance[planeY0,vert3]	],
	Max[	RegionDistance[planeX0,vert2],	RegionDistance[planeX0,vert3],	RegionDistance[planeX0,vert7],	RegionDistance[planeX0,vert6]	],
	Max[	RegionDistance[planeZ1,vert1],	RegionDistance[planeZ1,vert4],	RegionDistance[planeZ1,vert3],	RegionDistance[planeZ1,vert2]	],
	Max[	RegionDistance[planeY1,vert1],	RegionDistance[planeY1,vert2],	RegionDistance[planeY1,vert6],	RegionDistance[planeY1,vert5]	],
	Max[	RegionDistance[planeX1,vert1],	RegionDistance[planeX1,vert5],	RegionDistance[planeX1,vert8],	RegionDistance[planeX1,vert4]	]
],50],NumberFormat->(Row[{#1,"e ",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Title:: *)
(*Detailed  Calculations for Debugging*)


(* ::Subtitle:: *)
(*Line Minimum Height  :*)


Print["Line: ", TextString[ScientificForm[N[EuclideanDistance[vert1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Subtitle:: *)
(*Triangle Minimum Height :*)


Print["Tria_Face1: ", TextString[ScientificForm[N[RegionDistance[lineY0,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Tria_Face2: ", TextString[ScientificForm[N[RegionDistance[lineXY1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Tria_Face3: ", TextString[ScientificForm[N[RegionDistance[lineX0,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Subtitle:: *)
(*Quadrilateral Minimum Height :*)


Print["Quad_Face1.Node3: ", TextString[ScientificForm[N[RegionDistance[lineY0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face1.Node4: ", TextString[ScientificForm[N[RegionDistance[lineY0,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face2.Node4: ", TextString[ScientificForm[N[RegionDistance[lineX1,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face2.Node1: ", TextString[ScientificForm[N[RegionDistance[lineX1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face3.Node1: ", TextString[ScientificForm[N[RegionDistance[lineY1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face3.Node2: ", TextString[ScientificForm[N[RegionDistance[lineY1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face4.Node2: ", TextString[ScientificForm[N[RegionDistance[lineX0,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Quad_Face4.Node3: ", TextString[ScientificForm[N[RegionDistance[lineX0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Subtitle:: *)
(*Tetrahedra Minimum Height :*)


Print["Tetr_Face1: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Tetr_Face2: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Tetr_Face3: ", TextString[ScientificForm[N[RegionDistance[planeXYZ1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Tetr_Face4: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Subtitle:: *)
(*Pyramid Minimum Height :*)


Print["Pyra_Face1.Node5: ", TextString[ScientificForm[N[RegionDistance[planeXYZ1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pyra_Face2.Node3: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pyra_Face2.Node4: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pyra_Face3.Node1: ", TextString[ScientificForm[N[RegionDistance[planeXZ1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pyra_Face3.Node4: ", TextString[ScientificForm[N[RegionDistance[planeXZ1,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pyra_Face4.Node1: ", TextString[ScientificForm[N[RegionDistance[planeYZ1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pyra_Face4.Node2: ", TextString[ScientificForm[N[RegionDistance[planeYZ1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pyra_Face5.Node2: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pyra_Face5.Node3: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Subtitle:: *)
(*Prism Minimum Height :*)


Print["Pris_Face1.Node3: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face1.Node6: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert8],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pris_Face2.Node1: ", TextString[ScientificForm[N[RegionDistance[planeXY1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face2.Node4: ", TextString[ScientificForm[N[RegionDistance[planeXY1,vert5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pris_Face3.Node2: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face3.Node5: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert6],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pris_Face4.Node4: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face4.Node5: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert6],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face4.Node6: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert8],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Pris_Face5.Node1: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face5.Node3: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Pris_Face5.Node2: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]


(* ::Subtitle:: *)
(*Hexa Minimum Height:*)


Print["Hexa_Face1.Node5: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face1.Node6: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert6],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face1.Node7: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert7],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face1.Node8: ", TextString[ScientificForm[N[RegionDistance[planeZ0,vert8],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Hexa_Face2.Node4: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face2.Node8: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert8],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face2.Node7: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert7],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face2.Node3: ", TextString[ScientificForm[N[RegionDistance[planeY0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Hexa_Face3.Node2: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face3.Node3: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face3.Node7: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert7],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face3.Node6: ", TextString[ScientificForm[N[RegionDistance[planeX0,vert6],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Hexa_Face4.Node1: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face4.Node4: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face4.Node3: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face4.Node2: ", TextString[ScientificForm[N[RegionDistance[planeZ1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Hexa_Face5.Node1: ", TextString[ScientificForm[N[RegionDistance[planeY1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face5.Node2: ", TextString[ScientificForm[N[RegionDistance[planeY1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face5.Node6: ", TextString[ScientificForm[N[RegionDistance[planeY1,vert6],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face5.Node5: ", TextString[ScientificForm[N[RegionDistance[planeY1,vert5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]

Print["Hexa_Face6.Node1: ", TextString[ScientificForm[N[RegionDistance[planeX1,vert1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face6.Node5: ", TextString[ScientificForm[N[RegionDistance[planeX1,vert5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face6.Node8: ", TextString[ScientificForm[N[RegionDistance[planeX1,vert8],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]
Print["Hexa_Face6.Node4: ", TextString[ScientificForm[N[RegionDistance[planeX1,vert4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]]









