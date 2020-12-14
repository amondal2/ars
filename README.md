# ars
R package implementing an adaptive rejection sampler. 

Sample usage:
```
> devtools::load_all("./")
Loading ars
> ars(dnorm, 10)
 [1]  0.89259554 -1.40340976  0.04997463 -0.66866277 -0.43522863 -0.86795337  0.62318480 -0.05521466 -0.05610585
[10]  0.12313170
```

To run tests:
```
> testthat::test_package('ars')
```
