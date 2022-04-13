## Tailoring Nanoporous Graphene via Machine Learning: Predicting Probabilities and Formation Times of Arbitrary Nanopore Shapes

This repository contains the necessary datasets and Python files to train a two-stage machine learning (ML) model for predicting the probabilities and formation times of graphene nanopores.
The first stage in the two-stage model aims at segregating the nanopores (corresponding to both probability and formation time prediction) into majority/minority classes based on a fixed threshold for probability and formation time. And,
the second stage aims at predicting the actual probability/formation times of the nanopores based on the seperation of nanopores in stage 1. The complexity in the datasets have been captured using the gradient boosted ML models and
the linearity in the datasets have also been captured using linear regression.

![alt text](https://github.com/agrgroup/MLforNanopores/blob/main/TOC_image.png)

<pre>
The sequence recommended for executing the files using the command prompt is:
    1. python load.py - Concatenates the individual CSV files in the "physicalFeatures" folder into one dataframe and splits it 
    into training and test CSV files.
    [Optional if the required train.py and test.py files are already present in the "input" folder in the "probability" and 
    "formation_time" directories].
    
    2. python cv.py - Cross validates the two-stage model using any the KFold technique. This helps in knowing if the model is learning stably
    or suffering from any of bias or variance-related problems.
    
    3. python train.py - Trains the two-stage model (see model.py) using the training set and saves the model as a pickle file 
    (PKL/.pkl). Optionally, the tuners (see tuners.py) can also be enabled by uncommenting lines 41 and 42 in the "train.py" file corresponding 
    to probability prediction, and lines 35 and 36 in the same file but corresponding to formation time prediction. Npte that the two-stage model
    is trained on the complete training set and isn't an ensemble of models over different folds trained during cross validation.
</pre>

If you use our code for your research, please cite: 
