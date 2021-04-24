#############################################################################
##
##  Functions.gd        Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
#! @Chapter Functionality
##
#############################################################################

#! @Section Availability of Polymake

#! @Arguments
#! @Returns a boolean
#! @Description
#! Checks if the polymake functionality is available in Julia.
DeclareOperation( "PolymakeAvailability", [ ] );

#! @Arguments
#! @Returns a boolean
#! @Description
#! Checks if the polymake functionality is available and if yes loads it. For success, it returns true and false otherwise.
DeclareOperation( "LoadPolymake", [ ] );
