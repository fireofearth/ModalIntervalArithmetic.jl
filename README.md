# ModalIntervalArithmetic.jl

Modal interval arithmetic using algorithms defined in Chapter 5 of the book Modal Interval Analysis: New Tools for Numerical Information, by Miguel Sainz.

# Dependencies

[ModalIntervalArithmetic.jl](https://github.com/fireofearth/ModalIntervalArithmetic.jl)  

# Usage

Instantiate modal intervals of the form `ModalInterval{T<:Real} <: Number` using the following constructors

`ModalInterval(a::Real, b::Real)` constructs a modal interval `[a, b]`  
`ModalInterval{T}(a::Real)` constructs a modal interval `[a, a]` representing a real number, cast to type `T`  
`ModalInterval(a::Real)` constructs a modal interval `[a, a]` representing a real number  
`ModalInterval(X::Interval)` constructs a proper modal interval from classical interval `X`  
`ModalInterval{T}(X::Interval)` constructs a proper modal interval from classical interval `X`, with endpoints cast to type `T`  
`ModalInterval{T}(X::ModalInterval)`  

Julia operations are overloaded to support `ModalInterval`.

## TODOs

- Extend modal intervals to additional unary/binary operations (trig., exp.).
- Implement AE-extensions of real functions.
- Automatic inner and outer roundings of intervals.

## Related Works

Miguel Sainz et al. 2014. *Modal Interval Analysis: New Tools for Numerical Information*. Springer Lecture Notes in Mathematics, New York, NY.
PDF:<https://www.springer.com/gp/book/9783319017204>
