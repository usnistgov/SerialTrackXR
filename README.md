# SerialTrack-X-ray 2D

## Purpose
This repository contains scripts for running a custom implementation of ("SerialTrack",<https://github.com/FranckLab/SerialTrack>) designed for 2D X-ray projection images taken during in situ micro-X-ray computed tomography (XCT) experiments. The particle tracking-based technique allows for surface displacement and strain fields to be measured during, for example, mini-tensile test motions. This provides quantitative visualization of non-uniformities in the strain field (i.e., necking, shear banding), direct measurement of applied strain to supplement crosshead-based data, and quantification of system drift, among other uses. 

## Requirements and Content

### Dependencies
This is a Matlab-based code. It was developed using Releases 2021b - 2023b under Windows 10. The following toolboxes are used, although it may be straightforward to adapt to code to reduce these dependencies we have not tested this:
- Curve Fitting Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox
- System Identification Toolbox
- Wavelet Toolbox

### Directory structure
Since this is a derivative of the SerialTrack package, it follows similar workflow and usage and has similar content. The codebase has been streamlined and updated for the specific needs of the measurement task. The folder structure is as follows:

- `function` contains the library of supporting function scripts that SerialTrack relies on for pre-processing, running the algorithm, and post-processing. Licenses for 3rd party functions are included as .txt file.
- `imgFolder` contains subdirectories with images sequences that are the input data for the algorithm. These can be elsewhere, but since the user is usually prompted to select a file or folder and the selection window starts at the top level direction this is a convenient organizational method.
- `results` output folder for results, typically for automatically saved `.mat` files, `.fig` files saved by the user, etc. 
- `src` core source scripts to run SerialTrack.
- The `main_hardpar_inc_xray_particles.m` file is the script used to input experiment specifications and is the run-script to execute the full process. 

## Running the SerialTrack code
See also the `HowToRun_SerialTrack-XR.txt` document. 

Briefly, the major steps are as follows:
1. Prepare images in a directory you can find. Usually this will be an image sequence named in alphanumerical order of the time-series in which they were acquired. Typically image loading method "0" works well for this type of data.
2. If the images have not been brightfield corrected save the brightfield image in a directory you can find.
3. If you wish to use an image mask (although it is not typically required) prepare the mask according to the "Image binary mask file" mode you indent to use.
4. In the `main` script"
 1. Update the `MPTPara.xstep` variable with the micro-to-pixel ratio from your instrument.
 2. Update any tracking parameters according to your problem: usually the defaults are a good starting point to get a sense of the tracking performance, but, for example, for different particle densities `n_neighborsMax` and `n_neighborsMin` might be changed, or for lower resolution images `edge_width` and/or `grid_spaceing` might be reduced.
5. Run the script and follow the prompts.
6. Save output, visualize results.

## Other links and info
For the data used in the development and calibration of the model see: 
> Kafka, Orion L., Landauer, Alexander K., Benzing, Jake T., Moser, Newell H., Mansfield, Elisabeth, Garboczi, Edward J. (2023), Dataset for: A technique for in-situ displacement and strain measurement with laboratory-scale X-ray Computed Tomography, National Institute of Standards and Technology, https://doi.org/10.18434/mds2-3127

For a static, archival version of the code at version 1.0.0 see: https://doi.org/10.18434/mds2-3118

Please cite this code as:
> Landauer, Alexander K, Kafka, Orion L (2023), Code for SerialTrack-XR particle tracking, National Institute of Standards and Technology, https://doi.org/10.18434/mds2-3118


## Contact and support
For questions, please open a new entry in the "Issues" tab. If needed, you can also find authors' contact information via the associated paper (see above). 

The corresponding author is Orion Kafka (NIST MML Applied Chemicals and Materials Division). Alex Landauer (NIST MML Materials Measurement Science Division) is the primary developers of the code.

## License

NIST Software Licensing Statement

NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.
