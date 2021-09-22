######################
# 1: The Julia type for ToricVarieties
######################

struct normalToricVariety
           GapNTV::GapObj
           polymakeNTV::Polymake.BigObject
end
export normalToricVariety


######################
# 2: Generic constructors
######################

"""
    NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )

Construct the normal toric variety whose fan has ray generators `r` and maximal cones `c`.

# Examples
```julia-repl
julia> NormalToricVariety( [-1 5; 0 1; 1 0; 0 -1], [[1,2],[2,3],[3,4],[4,1]])
NormalToricVariety(GAP: <A toric variety of dimension 2>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x0000561365a4f3d0))
```
"""
function NormalToricVariety( rays::Matrix{Int}, cones::Vector{Vector{Int}} )
    # construct the toric variety in GAP
    gap_rays = GapObj( rays, recursive = true )
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )

    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} GAP.Globals.POLYMAKE_JULIA_MODULE.common.rows(Incidence.pm_incidencematrix)

    pmntv = GAP.Globals.POLYMAKE_JULIA_MODULE.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )

    # wrap it into a struct and return
    return normalToricVariety( variety, pmntv )
end
export NormalToricVariety

"""
    NormalToricVariety( v::GapObj )

Construct the Julia wrapper for a `GAP` toric variety `v`.
"""
function NormalToricVariety(GapNTV::GapObj)
   pmNTV = ntv_gap2polymake(GapNTV)
   return normalToricVariety(GapNTV, pmNTV)
end
export NormalToricVariety

"""
    ntv_gap2polymake( v::GapObj )

Convert a `GAP` toric variety `v` into a `Polymake` toric variety.
"""
function ntv_gap2polymake(GapNTV::GapObj)
    ff = GAP.Globals.Fan(GapNTV)
    R = GAP.Globals.RayGenerators(ff)
    MC = GAP.Globals.RaysInMaximalCones(ff)
    rays = Matrix{Int}(R)
    cones = [findall(x->x!=0, vec) for vec in [[e for e in mc] for mc in MC]]
    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} GAP.Globals.POLYMAKE_JULIA_MODULE.common.rows(Incidence.pm_incidencematrix)
    pmntv = GAP.Globals.POLYMAKE_JULIA_MODULE.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )
    return pmntv
end
export ntv_gap2polymake

"""
    ntv_polymake2gap( v::Polymake.BigObject )

Convert a `Polymake` toric variety `v` into a `GAP` toric variety.
"""
function ntv_polymake2gap(polymakeNTV::Polymake.BigObject)
    rays = Matrix{Int}(polymakeNTV.RAYS)
    gap_rays = GapObj( rays, recursive = true )
    cones = [findall(x->x!=0, v) for v in eachrow(polymakeNTV.MAXIMAL_CONES)]
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )
    return variety
end
export ntv_polymake2gap


######################
# 3: Special constructors
######################

"""
    projective_space( d::Int )

Construct the projective space of dimension `d`.

# Examples
```julia-repl
julia> projective_space(  2 )
NormalToricVariety(GAP: <A projective toric variety of dimension 2>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x0000562fbcc21b70))
```
"""
function projective_space( d::Int )
    # construct the projective space in gap
    variety = GAP.Globals.ProjectiveSpace( d )
    f = GAP.Globals.POLYMAKE_JULIA_MODULE.fan.normal_fan(GAP.Globals.POLYMAKE_JULIA_MODULE.polytope.simplex(d))
    pmntv = GAP.Globals.POLYMAKE_JULIA_MODULE.fulton.NormalToricVariety(f)
    
    # wrap it and return
    return normalToricVariety( variety, pmntv )
end
export projective_space


######################
# 4: Properties
######################


"""
    is_normal_variety( v::normalToricVariety )

Checks if the normal toric variety `v` is normal. (This function is somewhat tautological at this point.)

# Examples
```julia-repl
julia> is_normal_variety( projective_space(  2 ) )
true
```
"""
function is_normal_variety( v::normalToricVariety )
    return true
end
export is_normal_variety


"""
    is_affine( v::normalToricVariety )

Checks if the normal toric variety `v` is affine.

# Examples
```julia-repl
julia> is_normal_variety( projective_space(  2 ) )
false
```
"""
function is_affine( v::normalToricVariety )
    return v.polymakeNTV.AFFINE
end
export is_affine


"""
    is_projective( v::normalToricVariety )

Checks if the normal toric variety `v` is projective, i.e. if the fan of `v` is the the normal fan of a polytope.

# Examples
```julia-repl
julia> is_projective( projective_space(  2 ) )
true
```
"""
function is_projective( v::normalToricVariety )
    return v.polymakeNTV.PROJECTIVE
end
export is_projective


"""
    is_smooth( v::normalToricVariety )

Checks if the normal toric variety `v` is smooth.

# Examples
```julia-repl
julia> is_smooth( projective_space(  2 ) )
true
```
"""
function is_smooth( v::normalToricVariety )
    return v.polymakeNTV.SMOOTH
end
export is_smooth


"""
    is_complete( v::normalToricVariety )

Checks if the normal toric variety `v` is complete.

# Examples
```julia-repl
julia> is_complete( projective_space(  2 ) )
true
```
"""
function is_complete( v::normalToricVariety )
    return v.polymakeNTV.COMPLETE
end
export is_complete


"""
    has_torusfactor( v::normalToricVariety )

Checks if the normal toric variety `v` has a torus factor.

# Examples
```julia-repl
julia> has_torusfactor( projective_space(  2 ) )
false
```
"""
function has_torusfactor( v::normalToricVariety )
    return GAP.Globals.HasTorusfactor( v.GapNTV )::Bool
end
export has_torusfactor


"""
    is_orbifold( v::normalToricVariety )

Checks if the normal toric variety `v` is an orbifold.

# Examples
```julia-repl
julia> is_orbifold( projective_space(  2 ) )
true
```
"""
function is_orbifold( v::normalToricVariety )
    return GAP.Globals.IsOrbifold( v.GapNTV )::Bool
end
export is_orbifold


"""
    is_simplicial( v::normalToricVariety )

Checks if the normal toric variety `v` is simplicial. Hence, this function works just as `is_orbifold`. It is implemented for user convenience.

# Examples
```julia-repl
julia> is_simplicial( projective_space(  2 ) )
true
```
"""
function is_simplicial( v::normalToricVariety )
    return v.polymakeNTV.SIMPLICIAL
end
export is_simplicial


"""
    is_isomorphic_to_projective_space( v::normalToricVariety )

Checks if the normal toric variety `v` is isomorphic to projective space.

# Examples
```julia-repl
julia> is_isomorphic_to_projective_space( projective_space(  2 ) )
true
```
"""
function is_isomorphic_to_projective_space( v::normalToricVariety )
    return GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapNTV )::Bool
end
export is_isomorphic_to_projective_space


"""
    is_direct_product_of_projective_spaces( v::normalToricVariety )

Checks if the normal toric variety `v` is isomorphic to a direct product of projective space.

# Examples
```julia-repl
julia> is_direct_product_of_projective_spaces( projective_space(  2 ) )
true
```
"""
function is_direct_product_of_projective_spaces( v::normalToricVariety )
    return GAP.Globals.IsDirectProductOfPNs( v.GapNTV )::Bool
end
export is_direct_product_of_projective_spaces
