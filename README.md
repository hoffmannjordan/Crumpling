# Machine Learning in a data-limited regime: Augmenting experiments with synthetic data uncovers order in crumpled sheets
Jordan Hoffmann
Yohai Bar-Sinai
Lisa Lee
Jovana Andrejevic
Shruti Mishra
Shmuel M. Rubinstein
Chris H. Rycroft

# Research Article
Preprint Available At: https://arxiv.org/abs/1807.01437

# Code Availability
Mathematica code for the neural network. Put everything in the same directory and unzip the zip folder. For more training data, contact Yohai or Shmuel.

Trained weights can be downloaded here:
https://www.dropbox.com/s/8ph0ll2l2u7ao0g/Trained_Weights_NI.wlnet?dl=0 


Flat folding code written in C++ by Chris H. Rycroft (chr@seas.harvard.edu). 
Uses voro++_2d€” see: http://math.lbl.gov/voro++/

Chris H. Rycroft, Voro++: A three-dimensional Voronoi cell library in C++, Chaos 19, 041111 (2009).

Code compiles on Ubuntu. Once compiled, generates `flatfold_gen` executable. Take four command line arguments:
`./flatfold_gen NUM_FOLDS RAND_SEED FRACTION_RADIAL_FOLD ALL_SAME_BOOL`
Where `NUM_FOLDS` is the number of folds, the `RAND_SEED` is the random seed, `FRACTION_RADIAL_FOLD` is a value between 0 and 100 for the fraction of folds that are inward folds. `ALL_SAME_BOOL` is whether or not to always fold in the same direction. 

Running the executable generates two .dat files and a `fold.gnuplot` script. The `fold.gnuplot` script, when run, generates a fold_fig.png that is used by the Mathematica code. A python (`run.py`) script allows you to customize options and generate figures. 
```bash
> mkdir images
> python run.py > run.sh
> sh run.sh
```
# Article Description
In the paper, we used tools from machine learning to try to understand the physics of crumpled sheets. Below, we show an image of a rigid, flat-folded 
sheet along with a crumple sheet. Below, we show curvature maps from a flat-folded and crumpled sheet. The two images are not form the same sheets.
![Crumpled Sheet](../master/ims/sheets.png)
In the case of rigid flat folding, there are specific geometric rules that govern each vertex. Therefore, given either the ridges or the valleys, the location of
 the other can be inferred with near certainty. We asked to what degree this is true in crumpled sheets.
 Our goal is outlined below, where we try to train a network that is given the valleys to predict the distance from the nearest ridge.
 ![Crumpled Sheet](../master/ims/goal.png) 
 _to be continued_