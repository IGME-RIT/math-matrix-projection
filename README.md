# Matrix Projection and Rejection

Previously, vector projection was introduced and explored. However, this is only good to project onto one vector.

Now that we have matrices, the question must be asked "Can projection be represented as a matrix?"

What if we instead want to project onto a subspace spanned by multiple vectors? E.g. what if we want to project onto a 2D plane in R3?

Indeed, we will see that vector projection can be represented as a matrix, and the process can be extended.

# Setup

You will need to have CMake installed on your computer, and properly added to your path. In order to setup, run the following in a shell, then open the project in your preferred editor. Windows setup has been configured for use with Visual Studio.

Windows:
```
cd path/to/folder
setup.cmd
```
Linux:
```
cd path/to/folder
./setup
```
