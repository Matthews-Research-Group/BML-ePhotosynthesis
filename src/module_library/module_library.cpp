#include "module_library.h"
#include "../framework/module_creator.h"  // for create_mc

// Include all the header files that define the modules.
#include "example_module.h"
#include "c3_leaf_photosynthesis.h"
#include "multilayer_canopy_properties.h"
#include "multilayer_c3_canopy.h"
#include "multilayer_canopy_integrator.h"

creator_map BMLePhoto::module_library::library_entries =
{
    {"example_module", &create_mc<example_module>},
    {"c3_leaf_photosynthesis", &create_mc<c3_leaf_photosynthesis>},
    {"ten_layer_canopy_properties", &create_mc<ten_layer_canopy_properties>},
    {"ten_layer_c3_canopy", &create_mc<ten_layer_c3_canopy>},
    {"ten_layer_canopy_integrator", &create_mc<ten_layer_canopy_integrator>}
};
