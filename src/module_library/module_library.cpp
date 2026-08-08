#include "module_library.h"
#include "../framework/module_creator.h"  // for create_mc

// Include all the header files that define the modules.
#include "example_module.h"
#include "FvCB.h"
#include "c3_canopy.h"
#include "c3_leaf_photosynthesis.h"
#include "c3photo_module.h"
#include "multilayer_canopy_properties.h"
#include "multilayer_c3_canopy.h"
#include "multilayer_canopy_integrator.h"

creator_map BMLePhoto::module_library::library_entries =
{
    {"example_module", &create_mc<example_module>},
    {"FvCB", &create_mc<FvCB>},
    {"c3_canopy", &create_mc<c3_canopy>},
    {"c3_leaf_photosynthesis", &create_mc<c3_leaf_photosynthesis>},
    {"c3photo", &create_mc<c3photo_module>},
    {"ten_layer_canopy_properties", &create_mc<ten_layer_canopy_properties>},
    {"ten_layer_c3_canopy", &create_mc<ten_layer_c3_canopy>},
    {"ten_layer_canopy_integrator", &create_mc<ten_layer_canopy_integrator>}
};
