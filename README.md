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
    " output_suffix=_img_name" +
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

`input_stack` accepts either an image stack 
e.g. `input_stack=C:\\img\\dir\\ImageStack.tif`
or a comma separated list of image paths,
 e.g. `input_stack='C:\\img\\dir\\1_00001_Z007_CH1.tif,C:\\img\\dir\\1_00002_Z007_CH1.tif'`
 
Depeding on the parameters the output is saved into the following directories:
```
out_dir
 |-- corrected
 |    `-- corrected_img_name.tif
 |-- flatfield
 |    `-- flatfield_img_name.tif
 `-- darkfield
     `-- darkfield_img_name.tif
```
