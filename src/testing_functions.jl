const infs_vals    = Set([:no, :only_at_edges, :yes])
const missing_vals = Set([:no, :only_at_edges, :yes])
const NaN_vals     = Set([:no, :only_at_edges, :yes])
const monotonicity_vals = Set([:constant, :increasing, :decreasing, :nondecreasing, :nonincreasing, :quasiconcave, :quasiconvex, :any])
const curvature_vals    = Set([:concave, :convex, :linear, :any])

struct FunctionTestResult{T<:Real}
    func::Function
    clean_grid::Array{T,1}
    clean_evals::Array{T,1}
    first_diffs::Array{T,1}
    gradients::Array{T,1}
    pass_missing_test::Union{Missing,Bool}
    pass_nan_test::Union{Missing,Bool}
    pass_infs_test::Union{Missing,Bool}
    monotonicity_test::Union{Missing,Bool}
    curvature_test::Union{Missing,Bool}
end

function test_univariate_function(f::Function, grid::Union{Array,StepRangeLen,UnitRange};
            monotonicity::Symbol = :any, curvature::Symbol = :any, evals = f.(grid), tol = eps(),
            print_reports::Bool = false, infs::Symbol = :no, missings::Symbol = :no, NaNs::Symbol = :no, linearity_test_tol = 1e-2)
    if !(monotonicity in monotonicity_vals) error("Not a valid value for monotonicity in the test_univariate_function. You must input a value from the following list ", monotonicity_vals) end
    if !(curvature    in curvature_vals)    error("Not a valid value for curvature in the test_univariate_function. You must input a value from the following list ", curvature_vals) end
    if (sum(ismissing.(grid)) > 0) || (sum(isnan.(grid)) > 0) || (sum(isinf.(grid)) > 0) error("Grid can not be input containing missing, NaN or Inf values.") end
    grid = collect(grid)
    if !all(sort(grid) .== grid) error("A strictly increasing grid must be input") end # So now we are dealing strictly with a sorted array.
    clean_grid, clean_evals, pass_missing_test, pass_nan_test, pass_infs_test = test_existence(grid, evals, infs, missings, NaNs, print_reports)
    len   = length(clean_grid)
    first_diffs       = clean_evals[2:len] .- clean_evals[1:(len-1)]
    x_diffs           = clean_grid[ 2:len] .- clean_grid[ 1:(len-1)]
    gradients         = first_diffs ./ x_diffs
    monotonicity_test = test_monotonicity(clean_evals, first_diffs, monotonicity, tol)
    curvature_test    = test_curvature(gradients, curvature, tol, linearity_test_tol)
    return FunctionTestResult(f, clean_grid, clean_evals, first_diffs, gradients, pass_missing_test, pass_nan_test, pass_infs_test, monotonicity_test, curvature_test)
end

function test_existence(grid::Array, evals::Array, infs::Symbol, missings::Symbol, NaNs::Symbol, print_reports::Bool)
    N = length(evals)
    # Missings
    pass_missing_test = true
    where_missing = findall(ismissing.(evals))
    if length(where_missing) > 0
        if (missings == :no) || ((missings == :only_at_edges) && (N > 2) && (length(findall(ismissing.(evals[2:(N-1)]))) > 0))
            pass_missing_test = false
            if print_reports println("Failed missing test. Inputting the following values yields a missing ", print(grid[where_missing])) end
        end
    end
    # NaNs
    pass_nan_test = true
    where_nan = findall(isnan.(evals))
    if length(where_nan) > 0
        if (nans == :no) || ((nans == :only_at_edges) && (N > 2) && (length(findall(isnan.(evals[2:(N-1)]))) > 0))
            pass_nan_test = false
            if print_reports println("Failed NaN test. Inputting the following values yields a NaN ", print(grid[where_nan])) end
        end
    end
    # Infs
    pass_infs_test = true
    where_inf = findall(isinf.(evals))
    if length(where_inf) > 0
        if (infs == :no) || ((infs == :only_at_edges) && (N > 2) && (length(findall(isinf.(evals[2:(N-1)]))) > 0))
            pass_infs_test = false
            if print_reports println("Failed inf test. Inputting the following values yields an inf ", print(grid[where_inf])) end
        end
    end
    reduced_list = sort(unique([where_missing..., where_nan..., where_inf...]))
    obs_to_keep  = in.(1:N, Ref(reduced_list)) .== false
    return grid[obs_to_keep], evals[obs_to_keep], pass_missing_test, pass_nan_test, pass_infs_test
end

function test_monotonicity(evals::AbstractArray, first_diffs::AbstractArray, monotonicity::Symbol, tol::Real)
    N = length(evals)
    if N == 1 return true end
    monotonicity_test = true
    if monotonicity == :constant
        if maximum(evals) - minimum(evals) > tol monotonicity_test = false end
    elseif monotonicity == :increasing
        if minimum(first_diffs) < 0    monotonicity_test = false end
    elseif monotonicity == :decreasing
        if maximum(first_diffs) > 0    monotonicity_test = false end
    elseif monotonicity == :nondecreasing
        if minimum(first_diffs) < -tol monotonicity_test = false end
    elseif monotonicity == :nonincreasing
        if maximum(first_diffs) > tol monotonicity_test = false end
    elseif monotonicity == :quasiconcave
        max_index = findlast(evals .== maximum(evals))
        if sum(first_diffs[1:(max_index-1)] .< 0.0) > 0 monotonicity_test = false end
        if sum(first_diffs[max_index:(N-1)] .> 0.0) > 0 monotonicity_test = false end
    elseif monotonicity == :quasiconvex
        min_index = findlast(evals .== minimum(evals))
        if sum(first_diffs[1:(min_index-1)] .> 0.0) > 0 monotonicity_test = false end
        if sum(first_diffs[min_index:(N-1)] .< 0.0) > 0 monotonicity_test = false end
    elseif monotonicity == :any
        # Nothing to do, monotonicity test passes.
    else
        error("Unreachable code")
    end
    return monotonicity_test
end

function test_curvature(gradients::AbstractArray, curvature::Symbol, tol::Real, linearity_test_tol::Real)
    N = length(gradients)
    if N < 3 return true end
    if curvature == :any return true end
    curvature_test = true
    diff_gradients = gradients[2:N] .- gradients[1:(N-1)]
    if curvature == :concave
        if maximum(diff_gradients) > tol curvature_test = false end
    elseif curvature == :convex
        if minimum(diff_gradients) < -tol curvature_test = false end
    elseif curvature == :linear
        problem = abs(maximum(diff_gradients) - minimum(diff_gradients)) > linearity_test_tol
        if problem curvature_test = false end
    else
        error("Unreachable code")
    end
    return curvature_test
end

import Plots.plot
function plot(ftr::FunctionTestResult, num::Integer = 0)
    if num < 0 | num > 1 error("The input number must be 0 or 1") end
    if num == 0 return plot(ftr.clean_grid, ftr.clean_evals) end
    if num == 1
        len = length(ftr.gradients)
        return plot(ftr.clean_grid[1:len], ftr.gradients[1:len])
    end
end
function all_pass(ftr::FunctionTestResult)
    result = true
    if !(ismissing(ftr.pass_missing_test) || ftr.pass_missing_test)
        result = false
        println("Failings missing test")
    end
    if !(ismissing(ftr.pass_nan_test) || ftr.pass_nan_test)
        result = false
        println("Failings NaN test")
    end
    if !(ismissing(ftr.pass_infs_test) || ftr.pass_infs_test)
        result = false
        println("Failings Infs test")
    end
    if !(ismissing(ftr.monotonicity_test) || ftr.monotonicity_test)
        result = false
        println("Failings Monotonicity test")
    end
    if !(ismissing(ftr.curvature_test) || ftr.curvature_test)
        result = false
        println("Failings Curvature test")
    end
    return result
end
