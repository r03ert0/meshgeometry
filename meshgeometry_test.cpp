#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "meshgeometry.h"

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}

TEST_CASE( "Test linear algebra functions", "[linalg]") {
    SECTION("dot product") {
        REQUIRE( dot3D((float3D){1,2,3}, (float3D){1,1,1}) == 6 );
    }
    SECTION("addition") {
        float3D s = add3D((float3D){1,2,3}, (float3D){1,1,1});
        REQUIRE( (s.x == 2 && s.y == 3 && s.z == 4) );
    }
}
