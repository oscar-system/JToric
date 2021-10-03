######################
# 1: Attributes of ToricVarieties
######################

@doc Markdown.doc"""
    dim( v::AbstractNormalToricVariety )

Computes the dimension of the normal toric variety `v`.
"""
function dim( v::AbstractNormalToricVariety )
    return v.polymakeNTV.FAN_DIM::Int
end
export dim


@doc Markdown.doc"""
    dim_of_torusfactor( v::AbstractNormalToricVariety )

Computes the dimension of the torus factor of the normal toric variety `v`.
"""
function dim_of_torusfactor( v::AbstractNormalToricVariety )

    if has_torusfactor( v ) == false
        return 0
    end
    
    dimension_of_fan = v.polymakeNTV.FAN_DIM::Int
    ambient_dimension = v.polymakeNTV.FAN_AMBIENT_DIM::Int
    return ambient_dimension - dimension_of_fan
end
export dim_of_torusfactor


@doc Markdown.doc"""
    euler_characteristic( v::AbstractNormalToricVariety )

Computes the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic( v::AbstractNormalToricVariety )
    f_vector = Vector{Int}(v.polymakeNTV.F_VECTOR)
    return f_vector[ dim( v ) ]
end
export euler_characteristic


@doc Markdown.doc"""
    weil_divisors_of_variety( v::AbstractNormalToricVariety )

Compute the Weil divisors of the normal toric variety `v`.
"""
function weil_divisors_of_variety( v::AbstractNormalToricVariety )
    gap_divisors = GAP.Globals.WeilDivisorsOfVariety( v.GapNTV )
    return [ToricDivisor(extract_gap_divisor_coeffs(d), v) for d in gap_divisors ]
end
export weil_divisors_of_variety


struct ZariskiCotangentSheafViaEulerSequence
           GapZariskiCotangentSheafViaEulerSequence::GapObj
end
export ZariskiCotangentSheafViaEulerSequence


@doc Markdown.doc"""
    zariski_cotangent_sheaf_via_euler_sequence( v::AbstractNormalToricVariety )

Computes the Zariski cotangent sheaf of the normal toric variety `v` via the Euler sequence.
"""
function zariski_cotangent_sheaf_via_euler_sequence( v::AbstractNormalToricVariety )
    gap_ZariskiCotangentSheafViaEulerSequence = GAP.Globals.ZariskiCotangentSheafViaEulerSequence( v.GapNTV )
    return ZariskiCotangentSheafViaEulerSequence( gap_ZariskiCotangentSheafViaEulerSequence )
end
export zariski_cotangent_sheaf_via_euler_sequence


struct ZariskiCotangentSheafViaPoincareResidueMap
           GapZariskiCotangentSheafViaPoincareResidueMap::GapObj
end
export ZariskiCotangentSheafViaPoincareResidueMap


@doc Markdown.doc"""
    zariski_cotangent_sheaf_via_poincare_residue_map( v::AbstractNormalToricVariety )

Computes the Zariski cotangent sheaf of the normal toric variety `v` via the Poincare residue map.
"""
function zariski_cotangent_sheaf_via_poincare_residue_map( v::AbstractNormalToricVariety )
    gap_ZariskiCotangentSheafViaPoincareResidueMap = GAP.Globals.ZariskiCotangentSheafViaPoincareResidueMap( v.GapNTV )
    return ZariskiCotangentSheafViaPoincareResidueMap( gap_ZariskiCotangentSheafViaPoincareResidueMap )
end
export zariski_cotangent_sheaf_via_poincare_residue_map


#struct UnderlyingSheaf
#          GapUnderlyingSheaf::GapObj
#end
#export UnderlyingSheaf

#function underlying_sheaf( v::AbstractNormalToricVariety )
#   gap_Underlying = GAP.Globals.UnderlyingSheaf( v.GapNTV )
#    return UnderlyingSheaf( gap_UnderlyingSheaf )
#end
#export underlying_sheaf


"""
    nef_cone( v::NormalToricVariety )

Computes the nef cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> pp = projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> nef = nef_cone(pp)
A polyhedral cone in ambient dimension 1

julia> using Oscar

julia> rays(nef)
1-element VectorIterator{RayVector{Polymake.Rational}}:
 [-1]
```
"""
nef_cone( v::NormalToricVariety ) = Cone( v.polymakeNTV.NEF_CONE )
export nef_cone



"""
    mori_cone( v::NormalToricVariety )

Computes the mori cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> pp = projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> mori = mori_cone(pp)
A polyhedral cone in ambient dimension 1

julia> using Oscar

julia> rays(mori)
1-element VectorIterator{RayVector{Polymake.Rational}}:
 [-1]
```
"""
mori_cone( v::NormalToricVariety ) = Cone( v.polymakeNTV.MORI_CONE )
export mori_cone


######################
# 2: Methods of ToricVarieties
######################


@doc Markdown.doc"""
    ith_betti_number( v::AbstractNormalToricVariety, i::Int )

Compute the i-th Betti number of the normal toric variety `v`.
"""
function ith_betti_number( v::AbstractNormalToricVariety, i::Int )
    if isodd(i)
        return 0
    end
    k = div(i, 2)
    f_vector = Vector{Int}(v.polymakeNTV.F_VECTOR)
    pushfirst!(f_vector, 1)
    betti_number = sum( (-1)^(i-k) * binomial(i,k) * f_vector[ dim( v ) - i + 1 ] for i=k:dim( v ))
    return betti_number
    
end
export ith_betti_number


@doc Markdown.doc"""
    toric_ideal_binomial_generators(antv::AffineNormalToricVariety)

Get the exponent vectors corresponding to the generators of the toric ideal
associated to the affine normal toric variety `antv`.

# Examples
Take the cyclic quotient singularity corresponding to the pair of integers
`(2,5)`.
```jldoctest
julia> C = Oscar.positive_hull([-2 5; 1 0])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> toric_ideal_binomial_generators(antv)
pm::Matrix<long>
-1 -1 2 1
-1 0 3 -1
0 -1 -1 2
```
"""
function toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
    pmntv = pm_ntv(antv)
    result = pmntv.TORIC_IDEAL.BINOMIAL_GENERATORS
    return result
end
export toric_ideal_binomial_generators
