var documenterSearchIndex = {"docs":
[{"location":"#JToric-–-toric-geometry-in-Julia","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"","category":"page"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"CurrentModule = JToric","category":"page"},{"location":"#Goal","page":"JToric – toric geometry in Julia","title":"Goal","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"The goal of the package JToric is to make computations on toric geometry accessible in the OSCAR computer algebra system. As such, we follow the book Toric Varieties by David A. Cox, John B. Little and Henry K. Schenck]. Currently, this project is work-in-progress. Our goal is to ensure that one can perform all computations of Appendix B.","category":"page"},{"location":"#Generalities","page":"JToric – toric geometry in Julia","title":"Generalities","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"version","category":"page"},{"location":"#JToric.version","page":"JToric – toric geometry in Julia","title":"JToric.version","text":"version\n\nThe version number of the loaded JToric. Please mention this number in any bug report.\n\n\n\n\n\n","category":"constant"},{"location":"#Toric-Varieties","page":"JToric – toric geometry in Julia","title":"Toric Varieties","text":"","category":"section"},{"location":"#Constructors","page":"JToric – toric geometry in Julia","title":"Constructors","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"AffineNormalToricVariety(C::Cone)\nNormalToricVariety(C::Cone)\nNormalToricVariety(PF::PolyhedralFan)\nNormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )\nprojective_space\nhirzebruch_surface\ndel_pezzo","category":"page"},{"location":"#JToric.AffineNormalToricVariety-Tuple{Oscar.Cone}","page":"JToric – toric geometry in Julia","title":"JToric.AffineNormalToricVariety","text":"AffineNormalToricVariety(C::Cone)\n\nConstruct the affine normal toric variety U_C corresponding to a polyhedral cone C.\n\nExamples\n\nSet C to be the positive orthant in two dimensions.\n\njulia> C = Oscar.positive_hull([1 0; 0 1])\nA polyhedral cone in ambient dimension 2\n\njulia> antv = AffineNormalToricVariety(C)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"method"},{"location":"#JToric.NormalToricVariety-Tuple{Oscar.Cone}","page":"JToric – toric geometry in Julia","title":"JToric.NormalToricVariety","text":"NormalToricVariety(C::Cone)\n\nConstruct the (affine) normal toric variety X_Sigma corresponding to a polyhedral fan Sigma = C consisting only of the cone C.\n\nExamples\n\nSet C to be the positive orthant in two dimensions.\n\njulia> C = Oscar.positive_hull([1 0; 0 1])\nA polyhedral cone in ambient dimension 2\njulia> ntv = NormalToricVariety(C)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"method"},{"location":"#JToric.NormalToricVariety-Tuple{Oscar.PolyhedralFan}","page":"JToric – toric geometry in Julia","title":"JToric.NormalToricVariety","text":"NormalToricVariety(PF::PolyhedralFan)\n\nConstruct the normal toric variety X_PF corresponding to a polyhedral fan PF.\n\nExamples\n\nTake PF to be the normal fan of the square.\n\njulia> square = Oscar.cube(2)\nA polyhedron in ambient dimension 2\n\njulia> nf = Oscar.normal_fan(square)\nA polyhedral fan in ambient dimension 2\n\njulia> ntv = NormalToricVariety(nf)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"method"},{"location":"#JToric.NormalToricVariety-Tuple{Matrix{Int64}, Vector{Vector{Int64}}}","page":"JToric – toric geometry in Julia","title":"JToric.NormalToricVariety","text":"NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )\n\nConstruct the normal toric variety whose fan has ray generators r and maximal cones c.\n\nExamples\n\njulia> NormalToricVariety( [-1 5; 0 1; 1 0; 0 -1], [[1,2],[2,3],[3,4],[4,1]] )\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"method"},{"location":"#Oscar.Geometry.projective_space","page":"JToric – toric geometry in Julia","title":"Oscar.Geometry.projective_space","text":"projective_space( d::Int )\n\nConstruct the projective space of dimension d.\n\nExamples\n\njulia> projective_space( 2 )\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"function"},{"location":"#JToric.hirzebruch_surface","page":"JToric – toric geometry in Julia","title":"JToric.hirzebruch_surface","text":"hirzebruch_surface( r::Int )\n\nConstructs the r-th Hirzebruch surface.\n\nExamples\n\njulia> hirzebruch_surface( 5 )\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"function"},{"location":"#JToric.del_pezzo","page":"JToric – toric geometry in Julia","title":"JToric.del_pezzo","text":"delPezzo( b::Int )\n\nConstructs the delPezzo surface with b blowups for b at most 3.\n\nExamples\n\njulia> del_pezzo( 3 )\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\n\n\n","category":"function"},{"location":"#Conversion-among-GAP-and-Polymake-toric-varieties","page":"JToric – toric geometry in Julia","title":"Conversion among GAP and Polymake toric varieties","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"ntv_gap2polymake\nntv_polymake2gap","category":"page"},{"location":"#Properties-of-toric-varieties","page":"JToric – toric geometry in Julia","title":"Properties of toric varieties","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"isnormal\nisaffine\nisprojective\nissmooth\niscomplete\nhas_torusfactor\nis_orbifold\nissimplicial\nis_isomorphic_to_projective_space\nis_direct_product_of_projective_spaces\nis_gorenstein\nis_q_gorenstein\nis_fano","category":"page"},{"location":"#Hecke.isnormal","page":"JToric – toric geometry in Julia","title":"Hecke.isnormal","text":" isnormal(K::AnticNumberField) -> Bool\n\nReturns true if K is a normal extension of mathbb Q, false otherwise.\n\n\n\nisnormal(C::ClassField) -> Bool\n\nFor a class field C defined over a normal base field k, decide if C is normal over Q.\n\n\n\nisnormal(a::AlgAssAbsOrdIdl) -> Bool\n\nReturns true if a is a normal ideal and false otherwise.\n\n\n\nisnormal(G::T, H::T) where T <: GAPGroup\n\nReturn whether the subgroup H is normal in G, i. e., H is invariant under conjugation with elements of G.\n\n\n\n\n\nisnormal(A::MPolyQuo)\n\nGiven an affine algebra A over a perfect field, return true if A is normal, false otherwise.\n\nCAVEAT: The function computes the normalization of A. This may take some time.\n\n\n\nisnormal(P::Polyhedron)\n\nCheck whether P is normal.\n\nExamples\n\nThe 3-cube is normal.\n\njulia> C = cube(3)\nA polyhedron in ambient dimension 3\n\njulia> isnormal(C)\ntrue\n\nBut this pyramid is not:\n\njulia> P = convex_hull([0 0 0; 0 1 1; 1 1 0; 1 0 1]);\n\njulia> isnormal(P)\nfalse\n\n\n\n\n\nisnormal( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is normal. (This function is somewhat tautological at this point.)\n\nExamples\n\njulia> isnormal(projective_space( 2 ))\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.isaffine","page":"JToric – toric geometry in Julia","title":"JToric.isaffine","text":"isaffine( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is affine.\n\nExamples\n\njulia> isaffine( projective_space( 2 ) )\nfalse\n\n\n\n","category":"function"},{"location":"#JToric.isprojective","page":"JToric – toric geometry in Julia","title":"JToric.isprojective","text":"isprojective( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is projective, i.e. if the fan of v is the the normal fan of a polytope.\n\nExamples\n\njulia> isprojective( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#Hecke.issmooth","page":"JToric – toric geometry in Julia","title":"Hecke.issmooth","text":"issmooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem\n\nThrow an error if P is not a point of C, return false if P is a singular point of C, and true if P is a smooth point of C.\n\nExample\n\njulia> R, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])\n\njulia> C = Oscar.AffinePlaneCurve(x^2*(x+y)*(y^3-x^2))\nAffine plane curve defined by -x^5 - x^4*y + x^3*y^3 + x^2*y^4\n\njulia> P = Oscar.Point([QQ(0), QQ(0)])\nPoint with coordinates fmpq[0, 0]\n\njulia> Oscar.issmooth(C, P)\nfalse\n\n\n\nissmooth(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem\n\nThrow an error if P is not a point of C, return false if P is a singular point of C, and true if P is a smooth point of C.\n\nExample\n\njulia> S, (x, y, z) = PolynomialRing(QQ, [\"x\", \"y\", \"z\"])\n(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])\n\njulia> T, _ = grade(S)\n(Multivariate Polynomial Ring in x, y, z over Rational Field graded by\n  x -> [1]\n  y -> [1]\n  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])\n\njulia> C = Oscar.ProjPlaneCurve(x^2*(x+y)*(y^3-x^2*z))\nProjective plane curve defined by -x^5*z - x^4*y*z + x^3*y^3 + x^2*y^4\n\n\njulia> PP = projective_space(QQ, 2)\n(Projective space of dim 2 over Rational Field\n, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])\n\njulia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])\n(0 : 0 : 1)\n\njulia> Oscar.issmooth(C, P)\nfalse\n\n\n\nissmooth(P::Polyhedron)\n\nCheck whether P is smooth.\n\nExamples\n\nA cube is always smooth.\n\njulia> C = cube(8);\n\njulia> issmooth(C)\ntrue\n\n\n\n\n\nissmooth(PF::PolyhedralFan)\n\nDetermine whether PF is smooth.\n\nExamples\n\nEven though the cones of this fan cover the positive orthant together, one of these und thus the whole fan is not smooth.\n\njulia> PF = PolyhedralFan([0 1; 2 1; 1 0], IncidenceMatrix([[1, 2], [2, 3]]));\n\njulia> issmooth(PF)\nfalse\n\n\n\n\n\nissmooth( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is smooth.\n\nExamples\n\njulia> issmooth( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#Oscar.iscomplete","page":"JToric – toric geometry in Julia","title":"Oscar.iscomplete","text":"iscomplete(PF::PolyhedralFan)\n\nDetermine whether PF is complete, i.e. its support, the set-theoretic union of its cones, covers the whole space.\n\nExamples\n\nNormal fans of polytopes are complete.\n\njulia> iscomplete(normal_fan(cube(3)))\ntrue\n\n\n\n\n\niscomplete( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is complete.\n\nExamples\n\njulia> iscomplete( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.has_torusfactor","page":"JToric – toric geometry in Julia","title":"JToric.has_torusfactor","text":"has_torusfactor( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v has a torus factor.\n\nExamples\n\njulia> has_torusfactor( projective_space( 2 ) )\nfalse\n\n\n\n","category":"function"},{"location":"#JToric.is_orbifold","page":"JToric – toric geometry in Julia","title":"JToric.is_orbifold","text":"is_orbifold( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is an orbifold.\n\nExamples\n\njulia> is_orbifold( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.issimplicial","page":"JToric – toric geometry in Julia","title":"JToric.issimplicial","text":"issimplicial( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is simplicial. Hence, this function works just as is_orbifold. It is implemented for user convenience.\n\nExamples\n\njulia> issimplicial( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.is_gorenstein","page":"JToric – toric geometry in Julia","title":"JToric.is_gorenstein","text":"is_gorenstein( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is Gorenstein.\n\nExamples\n\njulia> is_gorenstein( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.is_q_gorenstein","page":"JToric – toric geometry in Julia","title":"JToric.is_q_gorenstein","text":"is_q_gorenstein( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is Q-Gorenstein.\n\nExamples\n\njulia> is_q_gorenstein( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.is_fano","page":"JToric – toric geometry in Julia","title":"JToric.is_fano","text":"is_fano( v::AbstractNormalToricVariety )\n\nChecks if the normal toric variety v is fano.\n\nExamples\n\njulia> is_fano( projective_space( 2 ) )\ntrue\n\n\n\n","category":"function"},{"location":"#Attributes-of-toric-varieties","page":"JToric – toric geometry in Julia","title":"Attributes of toric varieties","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"affine_open_covering\ncox_ring\nlist_of_variables_of_cox_ring\nclass_group\ntorus_invariant_divisor_group\nmap_from_character_to_principal_divisor\nmap_from_weil_divisors_to_class_group\ndim\ndim_of_torusfactor\ncoordinate_ring_of_torus\nlist_of_variables_of_coordinate_ring_of_torus\nis_product_of\nfactors\ncharacter_lattice\ntorus_invariant_prime_divisors\nirrelevant_ideal\nstanley_reisner_ideal\nmorphism_from_cox_variety\ncox_variety\nfan_of_variety\nfan\ncartier_torus_invariant_divisor_group\npicard_group\nname_of_variety\nzariski_cotangent_sheaf\ncotangent_sheaf\neuler_characteristic\nweil_divisors_of_variety\nzariski_cotangent_sheaf_via_euler_sequence\nzariski_cotangent_sheaf_via_poincare_residue_map\nnef_cone\nmori_cone\ntoric_ideal_binomial_generators(antv::AffineNormalToricVariety)","category":"page"},{"location":"#Hecke.class_group","page":"JToric – toric geometry in Julia","title":"Hecke.class_group","text":"class_group(K::AnticNumberField) -> GrpAbFinGen, Map\n\nShortcut for class_group(maximal_order(K)): returns the class group as an abelian group and a map from this group to the set of ideals of the maximal order.\n\n\n\nclass_group(O::NfOrd; bound = -1, method = 3, redo = false, large = 1000) -> GrpAbFinGen, Map\n\nReturns a group A and a map f from A to the set of ideals of O. The inverse of the map is the projection onto the group of ideals modulo the group of principal ideals. redo allows to trigger a re-computation, thus avoiding the cache. bound, when given, is the bound for the factor base.\n\n\n\n","category":"function"},{"location":"#AbstractAlgebra.Generic.dim","page":"JToric – toric geometry in Julia","title":"AbstractAlgebra.Generic.dim","text":"dim(Y::YoungTableau) -> BigInt\n\nReturn the dimension (using hook-length formula) of the irreducible representation of permutation group S_n associated the partition Y.part.\n\nSince the computation overflows easily BigInt is returned. You may perform the computation of the dimension in different type by calling dim(Int, Y).\n\nExamples\n\njulia> dim(YoungTableau([4,3,1]))\n70\n\njulia> dim(YoungTableau([3,1])) # the regular representation of S_4\n3\n\n\n\ndim(M::FreeModule{T}) where T <: FieldElement\n\nReturn the dimension of the given vector space.\n\n\n\ndim(N::Submodule{T}) where T <: FieldElement\n\nReturn the dimension of the given vector subspace.\n\n\n\ndim(N::QuotientModule{T}) where T <: FieldElement\n\nReturn the dimension of the given vector quotient space.\n\n\n\ndim(V::AbsSpace) -> Int\n\nReturn the dimension of the space V.\n\n\n\ndim(S::ZpGenus) -> fmpz\n\nReturn the dimension of this genus.\n\n\n\ndim(G::ZGenus) -> Int\n\nReturn the dimension of this genus.\n\n\n\ndim(I::MPolyIdeal)\n\nReturn the Krull dimension of I.\n\nExamples\n\njulia> R, (x, y, z) = PolynomialRing(QQ, [\"x\", \"y\", \"z\"])\n(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])\n\njulia> I = ideal(R, [y-x^2, x-z^3])\nideal(-x^2 + y, x - z^3)\n\njulia> dim(I)\n1\n\n\n\ndim(a::MPolyQuoIdeal)\n\nReturn the Krull dimension of a.\n\nExamples\n\njulia> R, (x, y, z) = PolynomialRing(QQ, [\"x\", \"y\", \"z\"], ordering=:degrevlex)\n(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])\n\njulia> Q, q = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]))\n(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^3*y^2 - x^2*y^3, x*y^4 - x*y^2), Map from\nMultivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^3*y^2 - x^2*y^3, x*y^4 - x*y^2) defined by a julia-function with inverse\n)\n\njulia> a = ideal(Q,[x^3*y^4-x+y, x*y+y^2*x])\nideal(x^3*y^4 - x + y, x*y^2 + x*y)\n\njulia> dim(a)\n1\n\n\n\ndim(Q::MPolyQuo)\n\nReturn the Krull dimension of the quotient ring Q.\n\nExamples\n\njulia> R, (x, y, z) = PolynomialRing(QQ, [\"x\", \"y\", \"z\"])\n(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])\n\njulia> Q, _ = quo(R, ideal(R, [y-x^2, x-z^3]))\n(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, x - z^3), Map from\nMultivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, x - z^3) defined by a julia-function with inverse\n)\n\njulia> dim(Q)\n1\n\n\n\ndim(C::Cone)\n\nReturn the dimension of C.\n\nExamples\n\nThe cone C in this example is 2-dimensional within a 3-dimensional ambient space.\n\njulia> C = Cone([1 0 0; 1 1 0; 0 1 0]);\n\njulia> dim(C)\n2\n\n\n\n\n\ndim(P::Polyhedron)\n\nReturn the dimension of P.\n\nExamples\n\njulia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];\n\njulia> P = convex_hull(V);\n\njulia> dim(P)\n2\n\n\n\n\n\ndim(PF::PolyhedralFan)\n\nReturn the dimension of PF.\n\nExamples\n\nThis fan in the plane contains a 2-dimensional cone and is thus 2-dimensional itself.\n\njulia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));\n\njulia> dim(PF)\n2\n\n\n\n\n\ndim( v::AbstractNormalToricVariety )\n\nComputes the dimension of the normal toric variety v.\n\n\n\n","category":"function"},{"location":"#JToric.dim_of_torusfactor","page":"JToric – toric geometry in Julia","title":"JToric.dim_of_torusfactor","text":"dim_of_torusfactor( v::AbstractNormalToricVariety )\n\nComputes the dimension of the torus factor of the normal toric variety v.\n\n\n\n","category":"function"},{"location":"#Hecke.picard_group","page":"JToric – toric geometry in Julia","title":"Hecke.picard_group","text":"picard_group(O::NfOrd) -> GrpAbFinGen, MapClassGrp\n\nReturns the Picard group of O and a map from the group in the set of (invertible) ideals of O.\n\n\n\npicard_group(O::AlgAssAbsOrd, prepare_ref_disc_log::Bool = false)\n  -> GrpAbFinGen, MapPicardGroup\n\nGiven an order O in a commutative algebra over mathbb Q, this function returns the picard group of O. If prepare_ref_disc_log is true, then (possibly expensive) preparations for the computation of refined discrete logarithms in non maximal orders are done.\n\n\n\n","category":"function"},{"location":"#JToric.euler_characteristic","page":"JToric – toric geometry in Julia","title":"JToric.euler_characteristic","text":"euler_characteristic( v::AbstractNormalToricVariety )\n\nComputes the Euler characteristic of the normal toric variety v.\n\n\n\n","category":"function"},{"location":"#JToric.nef_cone","page":"JToric – toric geometry in Julia","title":"JToric.nef_cone","text":"nef_cone( v::NormalToricVariety )\n\nComputes the nef cone of the normal toric variety v.\n\n\n\n\n\n","category":"function"},{"location":"#JToric.mori_cone","page":"JToric – toric geometry in Julia","title":"JToric.mori_cone","text":"mori_cone( v::NormalToricVariety )\n\nComputes the mori cone of the normal toric variety v.\n\n\n\n\n\n","category":"function"},{"location":"#JToric.toric_ideal_binomial_generators-Tuple{AffineNormalToricVariety}","page":"JToric – toric geometry in Julia","title":"JToric.toric_ideal_binomial_generators","text":"toric_ideal_binomial_generators(antv::AffineNormalToricVariety)\n\nGet the exponent vectors corresponding to the generators of the toric ideal associated to the affine normal toric variety antv.\n\nExamples\n\nTake the cyclic quotient singularity corresponding to the pair of integers (2,5).\n\njulia> C = Oscar.positive_hull([-2 5; 1 0])\nA polyhedral cone in ambient dimension 2\n\njulia> antv = AffineNormalToricVariety(C)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> toric_ideal_binomial_generators(antv)\npm::Matrix<long>\n-1 -1 2 1\n-1 0 3 -1\n0 -1 -1 2\n\n\n\n","category":"method"},{"location":"#Methods-of-toric-varieties","page":"JToric – toric geometry in Julia","title":"Methods of toric varieties","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"set_name_of_variety\ncoordinate_ring_of_torus( v::NormalToricVariety, names::Vector{String} )\ncox_ring( v::NormalToricVariety, name::String )\nBase.:*( v::NormalToricVariety, w::NormalToricVariety )\ncharacter_to_rational_function( l::Vector{Int}, v::NormalToricVariety )\nblowup_on_ith_minimal_torus_orbit( v::NormalToricVariety, i::Int )\nith_betti_number( v::NormalToricVariety, i::Int )\nnr_of_q_rational_points( v::NormalToricVariety, i::Int )","category":"page"},{"location":"#Base.:*-Tuple{NormalToricVariety, NormalToricVariety}","page":"JToric – toric geometry in Julia","title":"Base.:*","text":"*(x, y...)\n\nMultiplication operator. x*y*z*... calls this function with all arguments, i.e. *(x, y, z, ...).\n\nExamples\n\njulia> 2 * 7 * 8\n112\n\njulia> *(2, 7, 8)\n112\n\n\n\n\n\n","category":"method"},{"location":"#JToric.ith_betti_number-Tuple{NormalToricVariety, Int64}","page":"JToric – toric geometry in Julia","title":"JToric.ith_betti_number","text":"ith_betti_number( v::AbstractNormalToricVariety, i::Int )\n\nCompute the i-th Betti number of the normal toric variety v.\n\n\n\n","category":"method"},{"location":"#Toric-Divisors","page":"JToric – toric geometry in Julia","title":"Toric Divisors","text":"","category":"section"},{"location":"#Constructors-2","page":"JToric – toric geometry in Julia","title":"Constructors","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"ToricDivisor\ndivisor_of_character\ndivisor_of_class","category":"page"},{"location":"#JToric.ToricDivisor","page":"JToric – toric geometry in Julia","title":"JToric.ToricDivisor","text":"ToricDivisor(coeffs::AbstractVector, v::AbstractNormalToricVariety )\n\nConstruct the torus invariant divisor on the normal toric variety v as linear combination of the torus invariant prime divisors of v. The coefficients of thi linear combination are passed as list of integers as first argument.\n\nExamples\n\njulia> show( ToricDivisor( [1,1,2], projective_space( 2 ) ) )\nA torus invariant divisor on a normal toric variety\n\n\n\n","category":"type"},{"location":"#Properties-of-toric-divisors","page":"JToric – toric geometry in Julia","title":"Properties of toric divisors","text":"","category":"section"},{"location":"","page":"JToric – toric geometry in Julia","title":"JToric – toric geometry in Julia","text":"isample\nisbasepoint_free\niscartier\niseffective(td::ToricDivisor)\nisintegral(td::ToricDivisor) \nisnef\nis_primedivisor\nisprincipal\nisq_cartier\nisvery_ample\npolyhedron_of_divisor(td::ToricDivisor)","category":"page"},{"location":"#JToric.isample","page":"JToric – toric geometry in Julia","title":"JToric.isample","text":"isample(td::ToricDivisor)\n\nDetermine whether the toric divisor td is ample.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isample(td)\nfalse\n\n\n\n","category":"function"},{"location":"#JToric.isbasepoint_free","page":"JToric – toric geometry in Julia","title":"JToric.isbasepoint_free","text":"isbasepoint_free(td::ToricDivisor)\n\nDetermine whether the toric divisor td is basepoint free.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isbasepoint_free(td)\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.iscartier","page":"JToric – toric geometry in Julia","title":"JToric.iscartier","text":"is_cartier( d::ToricDivisor )\n\nChecks if the divisor d is Cartier.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> iscartier(td)\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.iseffective-Tuple{ToricDivisor}","page":"JToric – toric geometry in Julia","title":"JToric.iseffective","text":"iseffective(td::ToricDivisor)\n\nDetermine whether the toric divisor td is effective.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> iseffective(td)\ntrue\n\n\n\n","category":"method"},{"location":"#Hecke.isintegral-Tuple{ToricDivisor}","page":"JToric – toric geometry in Julia","title":"Hecke.isintegral","text":"isintegral(td::ToricDivisor)\n\nDetermine whether the toric divisor td is integral.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isintegral(td)\ntrue\n\n\n\n","category":"method"},{"location":"#JToric.isnef","page":"JToric – toric geometry in Julia","title":"JToric.isnef","text":"isnef(td::ToricDivisor)\n\nDetermine whether the toric divisor td is nef.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isnef(td)\ntrue\n\n\n\n","category":"function"},{"location":"#Hecke.isprincipal","page":"JToric – toric geometry in Julia","title":"Hecke.isprincipal","text":"isprincipal(A::NfOrdIdl) -> Bool, NfOrdElem\nisprincipal(A::NfOrdFracIdl) -> Bool, NfOrdElem\n\nTests if A is principal and returns (mathtttrue alpha) if A = langle alpharangle or (mathttfalse 1) otherwise.\n\n\n\nisprincipal(D::ProjCurveDivisor{S}) where S <: FieldElem\n\nReturn true if the divisor D is principal, and false otherwise\n\nExample\n\njulia> S, (x, y, z) = PolynomialRing(QQ, [\"x\", \"y\", \"z\"])\n(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])\n\njulia> T, _ = grade(S)\n(Multivariate Polynomial Ring in x, y, z over Rational Field graded by\n  x -> [1]\n  y -> [1]\n  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])\n\njulia> C = Oscar.ProjPlaneCurve(T(y^2*z - x*(x-z)*(x+3*z)))\nProjective plane curve defined by -x^3 - 2*x^2*z + 3*x*z^2 + y^2*z\n\n\njulia> PP = projective_space(QQ, 2)\n(Projective space of dim 2 over Rational Field\n, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])\n\njulia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])\n(0 : 1 : 0)\n\njulia> E = Oscar.ProjCurveDivisor(C, P)\n(0 : 1 : 0)\n\njulia> Oscar.isprincipal(E)\nfalse\n\n\n\nisprincipal(td::ToricDivisor)\n\nDetermine whether the toric divisor td is principal.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isprincipal(td)\nfalse\n\n\n\n","category":"function"},{"location":"#JToric.isq_cartier","page":"JToric – toric geometry in Julia","title":"JToric.isq_cartier","text":"isq_cartier(td::ToricDivisor)\n\nDetermine whether the toric divisor td is Q-Cartier.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isq_cartier(td)\ntrue\n\n\n\n","category":"function"},{"location":"#JToric.isvery_ample","page":"JToric – toric geometry in Julia","title":"JToric.isvery_ample","text":"isvery_ample(td::ToricDivisor)\n\nDetermine whether the toric divisor td is very ample.\n\nExamples\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isvery_ample(td)\nfalse\n\n\n\n","category":"function"},{"location":"#JToric.polyhedron_of_divisor-Tuple{ToricDivisor}","page":"JToric – toric geometry in Julia","title":"JToric.polyhedron_of_divisor","text":"polyhedron_of_divisor(td::ToricDivisor)\n\nConstruct the polyhedron P_D of a torus invariant divisor D=td as in 4.3.2 of CLS11. The lattice points of this polyhedron correspond to the global sections of the divisor.\n\nExamples\n\nThe polyhedron of the divisor with all coefficients equal to zero is a point, if the ambient variety is complete. Changing the coefficients corresponds to moving hyperplanes. One direction moves the hyperplane away from the origin, the other moves it across. In the latter case there are no global sections anymore and the polyhedron becomes empty.\n\njulia> H = hirzebruch_surface(4)\nA normal toric variety corresponding to a polyhedral fan in ambient dimension 2\n\njulia> td0 = ToricDivisor([0,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isfeasible(polyhedron_of_divisor(td0))\ntrue\n\njulia> dim(polyhedron_of_divisor(td0))\n0\n\njulia> td1 = ToricDivisor([1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isfeasible(polyhedron_of_divisor(td1))\ntrue\n\njulia> td2 = ToricDivisor([-1,0,0,0], H)\nA torus invariant divisor on a normal toric variety\n\njulia> isfeasible(polyhedron_of_divisor(td2))\nfalse\n\n\n\n","category":"method"}]
}
