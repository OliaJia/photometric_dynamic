                                                                                                                             //This is a test for generating data in the project of visibility subspace
//By Yaoguang Jia, Feb. 15th



#include "colors.inc"
#include "stones.inc"
#include "textures.inc"
#include "shapes.inc"
#include "glass.inc"
#include "metals.inc"
#include "woods.inc"

camera {
    location <0, 0, -3>
    look_at  <0, 0,  0>
  }
    
 
 sphere {
    <-1, 0.5, 0>, 0.5
    texture {
      pigment{White}
    }
  }
  sphere {
    <1, 0.5, 0>, 0.5
    texture {
      pigment{White}
    }
  } 
  
  
box {
    <-2, -0.5,   0>,  // Near lower left corner
    < 2,0, 4 >   // Far upper right corner
    texture {
      pigment{Tom_Wood}
                    // directions
    }
    rotate y*0     // Equivalent to "rotate <0,20,0>"
  }   
  
  
  
  


light_source { <0,4000, 0> color White}