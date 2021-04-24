#############################################################################
##
##  ExternalPolymakePolytope.gd  Convex package
##                               Martin Bies
##
##  Copyright 2021               University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
##  Chapter Polytopes in Polymake
##
#############################################################################


##############################################################################################
##
##  Section GAP category of PolymakePolytopes
##
##############################################################################################

DeclareRepresentation( "IsPolymakePolytopeRep", IsPolymakePolytope and IsAttributeStoringRep, [ ] );

BindGlobal( "TheFamilyOfPolymakePolytopes", NewFamily( "TheFamilyOfPolymakePolytopes" ) );

BindGlobal( "TheTypeOfPolymakePolytope", NewType( TheFamilyOfPolymakePolytopes, IsPolymakePolytopeRep ) );


##############################################################################################
##
##  Constructors for PolymakePolytopes
##
##############################################################################################


InstallGlobalFunction( Polymake_PolytopeByGenerators,
  function( arg )
    local poly, i, matrix, temp, dim;
    
    if Length( arg )= 0 or ForAll( arg, IsEmpty ) then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[1] ) then
        
        return Polymake_PolytopeByGenerators( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if ( not IsEmpty( arg[ 1 ] ) ) and not ( IsMatrix( arg[ 1 ] ) ) then
            Error( "Wronge input: The first argument should be a Gap matrix!" );
        fi;
        
        if ( not IsEmpty( arg[ 2 ] ) ) and not ( IsMatrix( arg[ 2 ] ) ) then
            Error( "Wronge input: The second argument should be a Gap matrix!" );
        fi;
        
        poly := rec( vertices := arg[ 1 ],
                     lineality := arg[ 2 ],
                     number_type := "rational",
                     rep_type := "V-rep" );
        ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
        return Polymake_CanonicalPolytopeByGenerators( poly );
        
    fi;
    
end );


InstallGlobalFunction( Polymake_PolytopeFromInequalities,
  function( arg )
    local poly, i, temp, matrix, dim;
    
    if Length( arg ) = 0 or ForAll( arg, IsEmpty ) then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[ 1 ] ) then
        
        return Polymake_PolytopeFromInequalities( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if ( not IsEmpty( arg[ 1 ] ) ) and not ( IsMatrix( arg[ 1 ] ) ) then
            Error( "Wronge input: The first argument should be a Gap matrix!" );
        fi;
        
        if ( not IsEmpty( arg[ 2 ] ) ) and not ( IsMatrix( arg[ 2 ] ) ) then
            Error( "Wronge input: The second argument should be a Gap matrix!" );
        fi;
        
        poly := rec( inequalities := arg[ 1 ],
                     equalities := arg[ 2 ],
                     number_type := "rational",
                     rep_type := "H-rep" );
        ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
        return Polymake_CanonicalPolytopeFromInequalities( poly );
        
    fi;
    
end );


##############################################################################################
##
##  Canonicalize polytopes
##
##############################################################################################

InstallMethod( Polymake_CanonicalPolytopeByGenerators,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, vertices, scaled_vertices, i, scale, lineality, scaled_lineality, new_poly;
    
    if poly!.rep_type = "H-rep" then
        
        return fail;
        
    else
        
        # compute vertices
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".VERTICES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        vertices := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational vertices - we turn them into integral vectors
        scaled_vertices := [];
        for i in [ 1 .. Length( vertices ) ] do
            scale := Lcm( List( vertices[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_vertices, [ scale * vertices[ i ] ] );
        od;
        
        # extract lineality
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational lineality - we turn them into integral vectors
        scaled_lineality := [];
        for i in [ 1 .. Length( lineality ) ] do
            scale := Lcm( List( lineality[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_lineality, [ scale * lineality[ i ] ] );
        od;
        
        # construct the new poly
        new_poly := rec( vertices := scaled_vertices,
                         lineality := scaled_lineality,
                         number_type := "rational",
                         rep_type := "V-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakeCone );
        return new_poly;
        
    fi;
    
end );

InstallMethod( Polymake_CanonicalPolytopeFromInequalities,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, ineqs, scaled_ineqs, i, scale, eqs, scaled_eqs, new_poly;
    
    if poly!.rep_type = "V-rep" then
        
        return fail;
        
    else
        
        # compute facets
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational facets - we turn them into integral vectors
        scaled_ineqs := [];
        for i in [ 1 .. Length( ineqs ) ] do
            scale := Lcm( List( ineqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_ineqs, [ scale * ineqs[ i ] ] );
        od;
        
        # compute affine hull
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".AFFINE_HULL" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational affine hulls - we turn them into integral vectors
        scaled_eqs := [];
        for i in [ 1 .. Length( eqs ) ] do
            scale := Lcm( List( eqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_eqs, [ scale * eqs[ i ] ] );
        od;
        
        # construct the new poly
        new_poly := rec( inequalities := scaled_ineqs,
                         equalities := scaled_eqs,
                         number_type := "rational",
                         rep_type := "H-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakeCone );
        return new_poly;
        
    fi;
    
end );


##############################################################################################
##
##  Conversion of polytopes
##
##############################################################################################

InstallMethod( Polymake_V_Rep,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, vertices, lineality, scaled_lineality, new_poly;
    
    if poly!.rep_type = "V-rep" then
        
        return poly;
        
    else
        
        # compute vertices
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".VERTICES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        vertices := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational vertices - we turn them into integral vectors
        scaled_vertices := [];
        for i in [ 1 .. Length( vertices ) ] do
            scale := Lcm( List( vertices[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_vertices, [ scale * vertices[ i ] ] );
        od;
        
        # compute lineality
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational lineality - we turn them into integral vectors
        scaled_lineality := [];
        for i in [ 1 .. Length( lineality ) ] do
            scale := Lcm( List( lineality[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_lineality, [ scale * lineality[ i ] ] );
        od;
        
        # construct the new poly
        new_poly := rec( vertices := scaled_vertices,
                         lineality := scaled_lineality,
                         number_type := "rational",
                         rep_type := "V-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakeCone );
        return new_poly;
        
    fi;
    
end );


InstallMethod( Polymake_H_Rep,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, ineqs, scaled_ineqs, eqs, scaled_eqs, new_poly;
    
    if poly!.rep_type = "H-rep" then
        
        return poly;
        
    else
        
        if poly!.rep_type = "V-rep" and poly!.matrix = [] then
            return Polymake_PolytopeFromInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # compute inequalities
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational facets - we turn them into integral vectors
        scaled_ineqs := [];
        for i in [ 1 .. Length( ineqs ) ] do
            scale := Lcm( List( ineqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_ineqs, [ scale * ineqs[ i ] ] );
        od;
        
        # compute equalities
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".AFFINE_HULL" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational affine hulls - we turn them into integral vectors
        scaled_eqs := [];
        for i in [ 1 .. Length( eqs ) ] do
            scale := Lcm( List( eqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_eqs, [ scale * eqs[ i ] ] );
        od;
        
        # construct the new poly
        new_poly := rec( inequalities := scaled_ineqs,
                         equalities := scaled_eqs,
                         number_type := "rational",
                         rep_type := "H-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakeCone );
        return new_poly;
        
    fi;
    
end );


##############################################################################################
##
##  Attributes of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_AmbientSpaceDimension,
              "finding the dimension of the ambient space of the poly",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Length( Polymake_H_Rep( poly )!.matrix[1] )-1;
    
end );


InstallMethod( Polymake_Dimension,
              " returns the dimension of the poly",
            [ IsPolymakePolytope ],
  function( poly )
    local help_poly, rays, lin, command_string, s;
    
    if Polymake_IsEmpty( poly ) then 
        return -1;
    fi;
    
    # compute v-representation
    help_poly := Polymake_V_Rep( poly );
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( help_poly ), ".CONE_DIM" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s ) - 1;
    
end );


InstallMethod( Polymake_GeneratingVertices,
              " return the list of generating vertices",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( Polymake_V_Rep( poly )!.generating_vertices );
    
end );


InstallMethod( Polymake_GeneratingRays,
              " return the list of generating vertices",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( Polymake_V_Rep( poly )!.generating_rays );
    
end );


InstallMethod( Polymake_Equalities,
              " return the list of equalities of a poly",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( ( Polymake_H_Rep( poly ) )!.equalities );
    
end );


InstallMethod( Polymake_Inequalities,
              " return the list of inequalities of a poly",
              [ IsPolymakePolytope ],
  function( poly )
    local equ, ineq;
    
    equ := ( Polymake_H_Rep( poly ) )!.equalities;
    ineq := ( Polymake_H_Rep( poly ) )!.inequalities;
    
    return Set( Concatenation( equ, (-1) * equ, ineq ) );
    
end );


InstallMethod( Polymake_LatticePoints,
              " return the list of the lattice points of poly",
              [ IsPolymakePolytope ],
  function( poly )
    local help_poly, command_string, s, res_string;
    
    # compute v-representation
    help_poly := Polymake_V_Rep( poly );
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( help_poly ), ".LATTICE_POINTS_GENERATORS" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    res_string := SplitString( s, '\n' );
    res_string := List( [ 2 .. Length( res_string ) - 3 ], i -> Concatenation( "[", ReplacedString( ReplacedString( res_string[ i ], " ", "," ), "<", "" ), "]" ) );
    return EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
    
end );


##############################################################################################
##
##  Properties of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_IsEmpty,
               "finding if the poly empty is or not",
               [ IsPolymakePolytope ],
  function( poly )
    
    return Length( Polymake_V_Rep( poly )!.matrix ) = 0;
    
end );


InstallMethod( Polymake_IsPointed,
               "finding if the poly is pointed or not",
               [ IsPolymakePolytope ],
  function( poly )
    local help_poly, rays, lin, command_string, s;
    
    # compute V-representation
    help_poly := Polymake_V_Rep( poly );
    
    # parse the rays into format recognized by Polymake
    command_string := Concatenation( Polymake_V_Rep_command_string( help_poly ), ".POINTED" );
    
    # issue command in Julia and fetch result as string
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


InstallMethod( Polymake_IsBounded,
              " returns if the polytope is bounded or not",
              [ IsPolymakePolytope ],
  function( poly )
    local help_poly, command_string, s;
    
    help_poly := Polymake_H_Rep( poly );
    command_string := Concatenation( Polymake_H_Rep_command_string( help_poly ), ".BOUNDED" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


##############################################################################################
##
##  Command strings
##
##############################################################################################

InstallMethod( Polymake_H_Rep_command_string,
               "construct command string for H-Representation of polytope in Julia",
               [ IsPolymakePolytope ],
  function( poly )
    local ineqs, eqs, command_string;
    
        # check if the given poly is a V-rep
        if not ( poly!.rep_type = "H-rep" ) then
            return fail;
        fi;
        
        # prepare string with inequalities
        ineqs := poly!.inequalities;
        ineqs := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "PolytopeByGAP4PackageConvex", " = Julia.Polymake.polytope.Polytope( INEQUALITIES = [ ", JoinStringsWithSeparator( ineqs, "; " ), " ] " );
        
        # check if we also need equalities
        eqs := poly!.equalities;
        if ( Length( eqs ) > 0 ) then
            eqs := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, " EQUATIONS = [ ", JoinStringsWithSeparator( eqs, "; " ), " ] " );
        fi;
        
        # append closing bracket
        command_string := Concatenation( command_string, ")" );
        
        # add return
        return command_string;
        
end );

InstallMethod( Polymake_V_Rep_command_string,
               "construct command string for V-Representation of polytope in Julia",
               [ IsPolymakePolytope ],
  function( poly )
    local vertices, lin, command_string;
        
        # check if the given poly is a V-rep
        if not ( poly!.rep_type = "V-rep" ) then
            return "fail";
        fi;
        
        # prepare string with vertices -- as polymake considers them affine, we have to add a 1 at the beginning
        vertices := poly!.matrix;
        vertices := List( [ 1 .. Length( vertices ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( vertices[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "PolytopeByGAP4PackageConvex", " = Julia.Polymake.polytope.Polytope( POINTS = [ ", JoinStringsWithSeparator( vertices, "; " ), "] " );
        
        # see if we need lineality
        lin := poly!.lineality;
        if ( Length( lin ) > 0 ) then
            lin := List( [ 1 .. Length( lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, " INPUT_LINEALITY = [ ", JoinStringsWithSeparator( lin, "; " ), " ] " );
        fi;
        
        # append closing bracket
        command_string := Concatenation( command_string, ")" );
        
        # add return
        return command_string;
        
end );
