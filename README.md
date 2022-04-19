# ML for nanopores

### Tailoring Nanoporous Graphene via Machine Learning: Predicting Probabilities and Formation Times of Arbitrary Nanopore Shapes
> Cite our paper here

## Contents

This repository contains the necessary MATLAB and Python scripts to manually generate the data and train a two-stage machine learning (ML) model, respectively, for predicting the probabilities and formation times of graphene nanopores.
The first stage in the two-stage model aims at segregating the nanopores (corresponding to both probability and formation time prediction) into majority/minority classes based on a fixed threshold for probability and formation time. Further, the second stage aims at predicting the actual probability/formation times of the nanopores based on the seperation of nanopores in stage 1. The complexity in the majority datasets have been captured using the gradient-boosted ML models and the linearity in the minority datasets have been captured using linear regression.

The dataset was generated using the DFT-KMC-Graph Theory approach proposed by Govind Rajan et al.[1] The codes for generating the datasets have been written in MATLAB. The MATLAB scripts to generate the datasets in the form of CSV files can be found in the folder named "generate_features." A total of 10000 KMC runs were carried out to obtain a database of nanopores for each value of vacancy size (N) between 4 and 22 (inclusive of both) considered. The simulated probabilities of the nanopores, for all the different N values were calculated after counting identical isomers obtained from the KMC simulations. The probability of a given shape is the number of times the isomer was seen in the KMC runs divided by 10,000. The formation time was also an output of the KMC simulation. 

![alt text](https://github.com/agrgroup/MLforNanopores/blob/main/TOC_image.png)

The sequence recommended for executing the MATLAB files using MATLAB's command window is:
* <pre>generate_isomers</pre> Generates nanopore isomers stochastically based on the KMC algorithm and stores their propoerties in the form of CSV files, which will further be used to generate the corresponding physical descriptors. The N values corresponding to the nanopores and the number of KMC runs must be pre-specified in order to run the simulation.


* <pre>generateCSV</pre> Generates the physical features for nanopores whose N values are specified based on the previously obtained KMC data. 

<b>Note</b>: Please refer to the souce code in [1] for a more detailed description of the used MATLAB scripts.


The sequence recommended for executing the Python files using the Windows command prompt is:
* <pre>python load.py</pre> Concatenates the individual CSV files in the "physicalFeatures" folder into one dataframe and splits it into training and test CSV files. (Executing this script is <b> not </b> necessary if the required train.csv and test.csv files are already present in the "input" folder in the "probability" and "formation_time" directories).


* <pre>python cv.py</pre> Cross validates the two-stage model using any the KFold technique. This helps in knowing if the model is learning stably or suffering from any bias or variance-related problems.


* <pre>python train.py</pre> Trains the two-stage model (see model.py) using the training set and saves the model as a pickle file (PKL/.pkl). Optionally, the tuners (see tuners.py) can also be enabled by uncommenting lines 41 and 42 in the "train.py" file corresponding to probability prediction, and lines 35 and 36 in the same file but corresponding to formation time prediction. The final two-stage model is trained on the complete training set and isn't an ensemble of models over the different folds trained during cross validation.

Note that the above suggested sequence for executing Python files must be followed for both targets, i.e., the nanopore probability and formation time, separately.

## References

<a id="1">[1]</a>  Govind Rajan, A.; Silmore, K.; Swett, J.; Robertson, A. W.; Warner, J. H.; Blankschtein, D.; Strano, M. S. Addressing the isomer cataloguing problem for nanopores in two-dimensional materials. Nature Materials, 18, pp 129–135 (2019). Source Code: https://github.com/agrgroup/nanopore_isomers.
