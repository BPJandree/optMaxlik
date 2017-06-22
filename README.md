
# optMaxlik
R package for minimization of Kullback-Leibler divergence of custom likelihood functions.

    library(devtools)
    install_github("BPJandree/optMaxlik")

Example code is contained in the [Vignette](https://github.com/BPJandree/optMaxlik/blob/master/optMaxlik.pdf)

You may also like my other package [AutoGLM](https://github.com/BPJandree/AutoGLM) for automated General to Specific modelling of generalized linear models suitable for very large datasets.

Feel free to contact me with suggestions or comments. This is the first version of the package, limited testing has been done.

# First update!

- Added an option to select parameters that will not be droppped;
- Fixed an issue where the algorithm crashes if nothing has to be dropped in the first iteration.
- If the numerical optimizer returns with non-zero exit status, op.maxLik will return the first optimization object.
