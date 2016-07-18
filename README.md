# Pi.jl
Various ways to compute π in Julia

I collected these for fun after a remark on the Julia mailing list. 

## Usage

```julia
include("src/pil.jl")
```

In the following, `n` is the requested number of digits of π. 

- `mpfr_pi(n)`: divert the problem to the `mpfr` library
- `atan_pi(n)`: compute π as $4 \arc \tan 1$ (probably the same as above)
- `spigot_pi(n)`: compute π one digit at the time
- `borwein2_pi(n)`: an algorithm from Borwein & Borwein
- `borwein_1984_pi(n)`: another algorithm from Borwein
- `chudnovsky2_1989_pi(n): an algorithm from Chudnovsky
- `gaussLegendre_pi(n)`: the classical algorithm
- `bbp_pi(n, start)`: compute π in hexidecimal digits, starting at an arbitrary offset, withou using arbitrary precision math

## Suggestions, improvements, etc

Are very welcome, please submit a pull request. 
