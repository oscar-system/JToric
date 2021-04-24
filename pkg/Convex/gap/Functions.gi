#############################################################################
##
##  Functions.gi        Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
## Chapter: Functionality
##
#############################################################################

## Section: Availability and loading of Polymake

##
InstallMethod( PolymakeAvailable, [  ],
  function( )
    local available;
    available := false;
    
    # Check if Polymake and the gap-julia interface are available
    if ( TestPackageAvailability( "JuliaInterface", ">= 0.5.2" ) <> fail ) then
        if IsPackageMarkedForLoading( "JuliaInterface", ">= 0.5.2" ) then
            if JuliaImportPackage("Polymake") then
                available := true;
                ImportJuliaModuleIntoGAP( "Polymake" );
            fi;
        fi;
    fi;
    
    return available;
    
end );

##
InstallMethod( CddInterfaceAvailable, [  ],
  function( )
    local available;
    available := false;
    
    if ( TestPackageAvailability( "CddInterface", ">= 2020.06.24" ) <> fail ) then
        if IsPackageMarkedForLoading( "CddInterface", ">= 2020.06.24" ) then
            available := true;
        fi;
    fi;
    
    return available;
    
end );

##
InstallMethod( NormalizInterfaceAvailable, [  ],
  function( )
    local available;
    available := false;
    
    if ( TestPackageAvailability( "NormalizInterface", ">= 1.2.0" ) <> fail ) then
        if IsPackageMarkedForLoading( "NormalizInterface", ">= 1.2.0" ) then
            available := true;
        fi;
    fi;
    
    return available;
    
end );
