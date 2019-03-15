# VectorEntrainment3D
A vector-based 3D rolling motion sediment entrainment model was developed with application to X-ray computed tomography (XCT) scanned riverbed grains.  VectorEntrainment3D extracts grain characteristics and calculates grain-to-grain contact points in an effort to estimate threshold of entrainment critical shear stress.  A vector-based 3D moment balance is used to calculate the entrainment threshold of individual surface grains. A cohesive force model is used to estimate resistance forces associated with coarse grain contact with a fine-grain matrix.  Once critical shear stress is determined, three angles are calculated, which describe the orientation of the grain's centre of mass and its two contact points: the bearing and tilt angles for the 3D rotation plane and the pivot angle contained therein.  Images of riverbed samples must be processed to separate grains and to extract the fine-grain matrix.  Two sets of image samples used in the following paper are also provided.

[Voepel, H., J. Leyland, R. Hodge, S. Ahmed, and D. Sear (2019),
Development of a vector-based 3D grain entrainment model with
application to X-ray computed tomography (XCT)scanned riverbed
sediment, *Earth Surface Processes and Landforms* doi: 10.1002/esp.4608](https://doi.org/10.1002/esp.4608)

## Getting Started

Clone this GitHub project to a local folder on a computer that has MATLAB installed.  Please note the following Prerequisites.

### Prerequisites

To run the 3D entrainment code, you will need MATLAB R2018a or higher installed with the following MATLAB Toolboxes:

* Image Processing Toolbox
* Mapping Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox

There are a suite of functions and subroutines that are run from the `main.m` script.  Follow the instructions in the `main.m` header.

## License

This project is licensed under the GNU General Public License version 3 - see the [LICENSE](LICENSE) file for details
