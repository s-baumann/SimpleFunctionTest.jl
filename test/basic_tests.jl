using SimpleFunctionTest

f(x) = (5-x)^2
grid = -5:0.01:5
# This should test there are no infs, missings, or NaNs but not test any convexity/monotonicity
ftr = test_univariate_function(f, grid) # This should test there are no infs, missings, or NaNs but not test any convexity/monotonicity
ftr.pass_missing_test
ftr.pass_nan_test
ftr.pass_infs_test
all_pass(ftr)
# Now doing a quasiconvex and convex tests.
ftr = test_univariate_function(f, grid, monotonicity = :quasiconvex, curvature = :convex)
ftr.monotonicity_test
ftr.curvature_test
# In this case there is also decreasing (as the domain ends at one). But it is not convex.
ftr = test_univariate_function(f, grid, monotonicity = :decreasing, curvature = :concave)
ftr.monotonicity_test
!ftr.curvature_test

# Now we should get an inf from 1/0.
g(x) = 1/x
ftr = test_univariate_function(g, grid)
!ftr.pass_infs_test
# This is not decreasing or increasing or concave or convex
ftr = test_univariate_function(g, grid, monotonicity = :decreasing, curvature = :concave)
!ftr.monotonicity_test
!ftr.curvature_test
ftr = test_univariate_function(g, grid, monotonicity = :increasing, curvature = :convex)
!ftr.monotonicity_test
!ftr.curvature_test
