(* ::Package:: *)

(* ::Title:: *)
(*Calculate  Polygon/Polyhedra Area/Volume from node coordinates*)


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
(*Define the cell polygons/polyhedra and calculate their area/volume*)


tria=Triangle[{vert1,vert2,vert4}];
quad=Polygon[{vert1,vert2, vert3,vert4}];
tetr=Tetrahedron[{vert1,vert2,vert4,vert5}];
pyra=Pyramid[{vert1,vert2,vert3,vert4,vert5}];
	pyraTetr1=Tetrahedron[{vert1,vert2,vert3,vert5}];
	pyraTetr2=Tetrahedron[{vert1,vert3,vert4,vert5}];
pris=Prism[{vert1, vert2, vert4, vert5, vert6, vert8}];
	prisTetr1=Tetrahedron[{vert1,vert2,vert4,vert5}];
	prisTetr2=Tetrahedron[{vert8,vert4,vert6,vert5}];
	prisTetr3=Tetrahedron[{vert6,vert2,vert5,vert4}];
hexa=Hexahedron[{vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8}];
	hexaTetr1=Tetrahedron[{vert1,vert2,vert4,vert5}];
	hexaTetr2=Tetrahedron[{vert3,vert2,vert7,vert4}];
	hexaTetr3=Tetrahedron[{vert6,vert2,vert5,vert7}];
	hexaTetr4=Tetrahedron[{vert8,vert5,vert4,vert7}];
	hexaTetr5=Tetrahedron[{vert5,vert2,vert4,vert7}];

Print["Line Leng   = ", TextString[ScientificForm[N[EuclideanDistance[vert1,vert2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["Tria Area   = ", TextString[ScientificForm[N[  Area[tria],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["Quad Area   = ", TextString[ScientificForm[N[  Area[quad],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["Tetr Volume = ", TextString[ScientificForm[N[Volume[tetr],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["Pyra Volume = ", TextString[ScientificForm[N[Volume[pyra],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["Pris Volume = ", TextString[ScientificForm[N[Volume[pris],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["Hexa Volume = ", TextString[ScientificForm[N[Volume[hexa],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print[];
Print["Decomposed Pyramid:"]
Print["Pyra Volume = ", TextString[ScientificForm[N[Volume[pyra],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["        Sum = ", TextString[ScientificForm[N[Volume[pyraTetr1]+Volume[pyraTetr2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  PyraTetr1 = ", TextString[ScientificForm[N[Volume[pyraTetr1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  PyraTetr2 = ", TextString[ScientificForm[N[Volume[pyraTetr2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print[];
Print["Decomposed Prism:"]
Print["Pris Volume = ", TextString[ScientificForm[N[Volume[pris],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["        Sum = ", TextString[ScientificForm[N[Volume[prisTetr1]+Volume[prisTetr2]+Volume[prisTetr3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  PrisTetr1 = ", TextString[ScientificForm[N[Volume[prisTetr1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  PrisTetr2 = ", TextString[ScientificForm[N[Volume[prisTetr2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  PrisTetr3 = ", TextString[ScientificForm[N[Volume[prisTetr3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print[];
Print["Decomposed Hexahedron:"]
Print["Hexa Volume = ", TextString[ScientificForm[N[Volume[hexa],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["        Sum = ", TextString[ScientificForm[N[Volume[hexaTetr1]+Volume[hexaTetr2]+Volume[hexaTetr3]+Volume[hexaTetr4]+Volume[hexaTetr5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  HexaTetr1 = ", TextString[ScientificForm[N[Volume[hexaTetr1],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  HexaTetr2 = ", TextString[ScientificForm[N[Volume[hexaTetr2],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  HexaTetr3 = ", TextString[ScientificForm[N[Volume[hexaTetr3],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  HexaTetr4 = ", TextString[ScientificForm[N[Volume[hexaTetr4],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];
Print["  HexaTetr5 = ", TextString[ScientificForm[N[Volume[hexaTetr5],50],NumberFormat->(Row[{#1,"e",#3}]&),NumberSigns->{"-"," "}]]];



