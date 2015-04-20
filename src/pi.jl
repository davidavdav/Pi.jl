## pil.jl (c) 2015 David A. van Leeuwen and others
## This file is licensed under the MIT software license.

## Various algorithms to compute π in many (hexa)decimal digits.

## We count digits after the period, e.g., for ndigits=2, π ≈ 3.14

nbits(ndigits::Int) = iceil(Int, log(2,10) * (ndigits+1))

## Ivar Nesje, diverting problem to mpfr 
function mpfr_pi(ndigits::Int)
    with_bigfloat_precision(nbits(ndigits)) do
        big(π)
    end
end

## seems to give identical result
function atan_pi(ndigits::Int)
    with_bigfloat_precision(nbits(ndigits)) do
        4atan(BigFloat(1))
    end
end

## spigot.  This computes π one decimal digit at the time.  
## Adapted from Hans W Borgers:
function spigot_pi(n)
    p = 0
    Pi = fill(uint8('0'), n+6)
    no_nines = 0
    d = n + 2
    N = ifloor(10*d/3.0 + 1)
    a = fill(2, N+1)
    ci = 1
    for l in 1:d
        for i=1:length(a)
            a[i] *= 10
        end
        for i in (N+1):-1:2
            j = 2*i - 1
            q, a[i] = divrem(a[i], j)
            a[i-1] += q * (i-1)
        end
        q, a[1] = divrem(a[1], 10)
        if q < 9
            Pi[ci] = uint8('0')+p
            for i = 1:no_nines
                Pi[ci+i] = uint8('9')
            end
            ci += no_nines+1
            p = q
            no_nines = 0
        elseif q == 9
            no_nines += 1
        elseif q == 10
            p += 1
            Pi[ci] = uint('0')+p
            for i = 1:no_nines
                Pi[ci+i] = uint8('0')
            end
            ci += no_nines+1
            p = 0
            no_nines = 0
        else
            error("spigot_pi: algorithm error!")
        end
    end
    Pi[1], Pi[2] = Pi[2], uint8('.') # remove first 0 and insert decimal point
    ASCIIString(Pi[1:n+2])
end

## from numerical recipes, Borwein & Borwein
function borwein2_pi(ndigits::Int)
    with_bigfloat_precision(nbits(ndigits)) do
        ɛ = BigFloat(10) ^ -ndigits
        x = √ BigFloat(2)
        p = 2 + x
        y = sx = √x
        while true
            lastp = p
            x = 0.5(sx + 1/sx)
            p *= (x+1)/(y+1)
            if abs(p-lastp) < ɛ break end
            sx = √x
            y = (y*sx + 1/sx) / (y+1)
        end
        p
    end
end

## from Wikipedia http://en.wikipedia.org/wiki/Borwein%27s_algorithm
function borwein_1984_pi(ndigits::Integer)
    with_bigfloat_precision(nbits(ndigits)) do
        ɛ = BigFloat(10) ^ -ndigits
        a = √ BigFloat(2)
        b = BigFloat(0)
        p = a + 2
        while true
            lastp = p
            sa = √a
            b = (b+1)sa / (a+b)
            a = 0.5(sa + 1/sa)
            p *= (1+a)b / (1+b)
            if abs(p-lastp) < ɛ break end           
        end
        p
    end
end

function chudnovsky2_1989_pi(ndigits::Int)
    with_bigfloat_precision(nbits(ndigits)) do
        ɛ = BigFloat(10) ^ -ndigits
        threek = sixk = k = s = BigFloat(0)
        sign = k!3 = threek! = sixk! = BigFloat(1)
        denfact = BigFloat(640320) ^ 3
        den = √denfact
        while true
            lasts = s
            s += sign * sixk! * (13591409 + 545140134k) / (threek! * k!3 * den)
            if abs(lasts - s) < ɛ break end
            k += 1
            for i=1:3
                threek += 1
                threek! *= threek
            end
            for i=1:6
                sixk += 1
                sixk! *= sixk
            end
            k!3 *= k*k*k
            den *= denfact
            sign = -sign
        end
        1 / 12s
    end
end

## Johan Sigfrids
function gaussLegendre_pi(ndigits::Integer)
    with_bigfloat_precision(nbits(ndigits)) do
        n = iceil(log2(ndigits))
        a = BigFloat(1)
        b = a / √ BigFloat(2)
        t = a / BigFloat(4)
        x = a
        while n > 0
            y = a
            a = (a + b) / 2
            b = √(b * y)
            t -= x * (y - a)^2
            x *= 2 
            n -= 1
        end
        (a + b)^2 / (4 * t)
    end
end


## helper function for BBP alorithm
function modpow(b, n::Integer, c)
    ## wikipedia
    r = 1 % c
    b %= c
    while n>0
        if isodd(n)
            r = (r*b) % c
        end
        n >>= 1
        b = (b*b) % c
    end
    return r
    ## Bailey, Borwein, Borwein, Plouffe (incorrect?)
    t = nextpow2(n+1) >> 1
    r = 1
    while true
        if n ≥ t
            r = (b*r) % c
            n -= t
        end
        t = div(t,2)
        if t < 1
            break
        end
        r = (r*r) % c
    end
    return r
    ## Straightforward (slow)
    r = 1
    while n>0
        r = (r*b) % c
        n -= 1
    end
    r
end

## compute 1 digit according to the Borwein, Borwein, and Plouffe algorithm.
## http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/pi-quest.pdf
function bbp_pi_digit(n::Int)
    if n==0
        return 3
    else
        n -= 1
    end
    const o = [1, 4, 5, 6]
    const w = [4, -2, -1, -1]
    frac = 0.
    for k=0:n
        for i=1:4
            den = 8k+o[i]
            frac += (w[i] * modpow(16, n-k, den)) / den
        end
    end
    ifloor(mod(frac, 1) * 16)::Int     
end

## The Borwein, Borwein and Plouffe formula of a hex digit of pi.
## This can be run in parallel, and can start at any digit.
## cf. table 3 from the publication: `bbp_pi(14, 10^6)`
function bbp_pi(n::Int, start::Int=0)
    if start==0
        range = 0:n
    else
        range = start + (0:n-1)
    end
    digits = @parallel (vcat) for i in range
        digit = bbp_pi_digit(i)
        if digit < 10
            d = '0'+digit
        else
            d = 'a'+digit-10
        end
        if i==0
            uint8(vcat(d,'.'))
        else
            uint8(d)
        end
    end
    ASCIIString(vcat(digits))
end

function bbp_pi_float(n::Int)
    digits = @parallel (vcat) for i in 0:n
        digit = bbp_pi_digit(i)
        nibble = BitArray(4)
        for i = 1:4
            nibble[5-i] = isodd(digit)
            digit >>= 1
        end
        nibble
    end
    bigfloat(digits, 2)
end

## This function interprets a bitvector as mantissa for a BigFloat number.
## Admittedy, somewhat clumsy. 
function bigfloat(x::BitVector, exp=1)
    ## remove leading falses and reverse
    nb = length(x)
    x = reverse(x)
    while length(x)>1 && !last(x)
        pop!(x)
    end
    if length(x)==1
        return with_bigfloat_precision(()->big(0.), nb)
    end
    ## shift bits until a leading true is found
    nb = length(x)
    firstbit = reinterpret(Uint, typemin(Int)) ## 0x8000000000000000
    while last(x.chunks) & firstbit == 0
        insert!(x, 1, false)
    end
    with_bigfloat_precision(nb) do
        f = big(1.)
        f.exp = exp
        for i=1:length(x.chunks)
            unsafe_store!(f.d, x.chunks[i], i)
        end
        return f
    end
end

# Helper to check bbp algorithm, converts float to hexadecimal string. 
function base16(x, n::Int)
    digits = Uint8[]
    nbefore = max(0,ifloor(log(16,x)))+1
    for i=1:nbefore
        x *= 16
    end
    for i = -nbefore:n
        if i==0
            push!(digits, '.')
            continue
        end
        d, r = divrem(x, 16)
        if d<10 
            push!(digits, '0'+d)
        else
            push!(digits, 'a'+d-10)
        end
        x = 16r
    end
    ASCIIString(digits)
end
