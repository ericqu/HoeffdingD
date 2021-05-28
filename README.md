# HoeffdingD.jl
HoeffdingD.jl implements in pure Julia the Hoeffding measure of dependence as described in the original paper: [A Non-Parametric Test of Independence](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-19/issue-4/A-Non-Parametric-Test-of-Independence/10.1214/aoms/1177730150.full) in particular chapter 5. The package also implements the D-test of independence described in chapter 9 of the same paper.

The advantage of this statistic is to detect nonlinear relationships that Pearson's correlation or Spearman's rank correlation are unable to detect.

## Installation
Enter the Pkg REPL by pressing ```]``` from the Julia REPL. Then install the package with: ```pkg> add https://github.com/ericqu/HoeffdingD```

## Usage

Here we demonstrate the classic example of detecting linear and quadratic relationships with Hoeffding measure contrasted with Perason Correlation and Spearman's rank correlations.

### Data generation
'''julia 
x = -2:0.1:2
linear_f(x) = 2x ; quad_f(x) = x^2
y_linear = linear_f.(x)
y_quad = quad_f.(x)
'''



## Reference
links