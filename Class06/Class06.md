# Class6:R functions
Zixuan Zeng (A16142927)

- [Our first (silly) function](#our-first-silly-function)
- [A second function](#a-second-function)
- [a protein generating function](#a-protein-generating-function)

All functions in R have at least 3 things:

- A **name**, we oucj this and use it to call our function, -Input
  **arguments** (there can be multiple) The **body** lines o R code that
  do the work

## Our first (silly) function

write a function to add some numbers

``` r
add <- function(x,y=1) {
  x + y
}
```

Now we can call the function:

``` r
add(c(10,10),100)
```

    [1] 110 110

``` r
add(10,100)
```

    [1] 110

## A second function

Write a function to generate random nucleotide sequences of a user
specified length:

THe `sample()` function can be helpful here.

``` r
V <- sample(c("A","C","G","T"),size=50, replace = TRUE)
```

I want the a 1 element long character vector that looks “TCATTG” not “T”
“C” “A” “T” “T” “G”

``` r
V <- sample(c("A","C","G","T"),size=50, replace = TRUE)
paste(V,collapse= "")
```

    [1] "GCATTTATCCAACGTACTTAGAGGCGACTCCTCCGCCGGTAACAGCTTCT"

Turn this into generate_DNA

``` r
generate_DNA <- function(size = 50){
V <- sample(c("A","C","G","T"),size=size, replace = TRUE)
paste(V,collapse= "")
}
```

Test it:

``` r
generate_DNA(60)
```

    [1] "GTTATAGATGACTAGTAACAATGCTATCGGAAGTCAAAGTCATAAGGACTACCAATTACT"

``` r
fasta <- FALSE
if(fasta){
  cat("HELLO You!")
} else{
  cat("No you dont!")
}
```

    No you dont!

Add the ability to return a multi-element vector or a single element
fasta like vector.

``` r
generate_DNA <- function(size = 50, fasta=TRUE){
V <- sample(c("A","C","G","T"),size=size, replace = TRUE)
s <- paste(V,collapse= "")

if(fasta){
  return(s)
} else{
  return(V)
}
}
```

Test:

``` r
generate_DNA(60,fasta=FALSE)
```

     [1] "A" "C" "C" "G" "T" "A" "T" "G" "C" "A" "A" "A" "G" "A" "G" "A" "T" "G" "G"
    [20] "C" "A" "T" "C" "C" "T" "G" "G" "A" "G" "T" "T" "C" "T" "G" "G" "A" "C" "C"
    [39] "C" "T" "G" "C" "G" "T" "G" "A" "C" "T" "G" "C" "C" "A" "A" "C" "T" "T" "C"
    [58] "G" "A" "C"

## a protein generating function

``` r
generate_protein <- function(size = 50, fasta=TRUE){
V <- sample(c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),size=size, replace = TRUE)
s <- paste(V,collapse= "")

if(fasta){
  return(s)
} else{
  return(V)
}
}
generate_protein(6)
```

    [1] "DMSHMY"

Use our new `generate_protein`function to make random protein sequences
of length 6 to 12 (i.e. one length 6, one length 7, up to length 12)

brute force:

``` r
generate_protein(6)
```

    [1] "YWYHIQ"

``` r
generate_protein(7)
```

    [1] "AQESSIG"

``` r
generate_protein(8)
```

    [1] "SITSRFGL"

A second way is to use a `for()` loop:

``` r
lengths <- 6:12
lengths
```

    [1]  6  7  8  9 10 11 12

``` r
for(i in lengths){
  cat(">",i,"\n",sep="")
  aa <- generate_protein(i)
  cat(aa)
  cat("\n")
}
```

    >6
    GQDEGP
    >7
    SEMQFQM
    >8
    RYQFVQPW
    >9
    VYKYWYLTM
    >10
    NNHKTRRATS
    >11
    GHVRNAKKWWA
    >12
    AWEQHMMDFHDF

A thrid, and better, way to solve this is to use the `apply()` family of
functions,specifically the `sapply()` function in this case

``` r
sapply(6:12,generate_protein)
```

    [1] "NRGYQS"       "HLYHKGC"      "DHEISTYH"     "DRSRQTSRR"    "LDMWWPHQLP"  
    [6] "GTTQMCGSHRM"  "CPILTLGRRFST"
