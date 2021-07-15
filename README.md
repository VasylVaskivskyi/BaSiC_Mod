## BaSiC Mod

Modified version of ImageJ/fiji plugin BaSiC (Background and Shading Correction) that can run in headless mode.


The original package is here:
https://github.com/marrlab/BaSiC

### Installation

1. Copy jars from `dependencies` folder (or `jars` folder in the repo) into fiji `jars` folder.
2. Copy `BaSiC_Mod.jar` into fiji `plugins` folder.


### Example .ijm macro
```
run("BaSiC Mod",
    "input_stack={path_to_stack}" +
    " flat-field_image_path=[]" +
    " dark-field_image_path=[]" +
    " output_dir={out_dir}" +
    " shading_estimation=[Estimate shading profiles]" +
    " shading_model=[Estimate flat-field only (ignore dark-field)]" +
    " setting_regularisation_parameters=Automatic" +
    " temporal_drift=Ignore" +
    " correction_options=[Compute shading only]" +
    " lambda_flat=0.500" +
    " lambda_dark=0.500");
    
run("Quit");
eval("script", "System.exit(0);");
```    