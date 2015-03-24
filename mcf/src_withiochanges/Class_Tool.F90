      MODULE Class_Tool
      
        USE mcf_header
        
        IMPLICIT NONE
        SAVE
        
        TYPE Tool
           INTEGER      :: flag
        END TYPE Tool
        
        
        INTERFACE tool_new
           MODULE PROCEDURE tool_init           
        END INTERFACE
        
      CONTAINS
        
#include "src/tool/tool_new.F90"
#include "src/tool/tool_uppercase.F90"
#include "src/tool/tool_cross_product.F90"
#include "src/tool/tool_rotation_matrix.F90" 
#include "src/tool/tool_barycentric_coordinate.F90"
#include "src/tool/tool_arbitrary_point_picking.F90"
#include "src/tool/tool_quadrilateral_point_picking.F90"
#include "src/tool/tool_wedge_point_picking.F90"
#include "src/tool/tool_triangle_point_picking.F90"
#include "src/tool/tool_tetrahedron_point_picking.F90"
        
      END MODULE Class_Tool
