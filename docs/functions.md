The ar() function
=================

Other than constructors for classes, the `ar()` function is intended as the
primary API for library. While it _is_ possible to directly call the underlying
operations (wedge, dot, division etc) it very quickly leads to noisy code and
becomes more of an exercise in programming than one of mathematics and physics.
As a matter of fact, it is also possible to use the `ar()` function to create
instances of Alphas and symbolic Pairs with a more concise syntax:

```Python
ai  # Equivalent to Alpha(i)
pi  # Equivalent to Pair(i)
```

Under the hood, `ar()` is actually a callable class that contains a lexer and
parser for a lightweight DSL that closely resembles the mathematical syntax used
in the published papers so far. The following basic operations are supported:

```Python
a ^ b  # wedge product of a and b
a . b  # dot product of a and b
a * b  # full (or geometric) product of a and b
a / b  # a divided by b
a \ b  # a divided into b

<a>n   # The grade(n) projection of a.
       # This selects out the elements of a that have an alpha
       # of the grade(n) (think of n as the length of the index).
       # If used on a Pair or Alpha, this returns the argument if
       # is of grade(n).
```

### Example usage
The following are several short examples of how to use the `ar()` function to
carry out computations in arpy. Note that variables defined in a REPL session
(or in the enclosing scope within a function) are useable within and `ar()`
function call.

#### Finding the product and quotient of Alphas and Pairs
```Python
>>> ar('a12 ^ a23')
α31

>>> ar('a012 ^ p31')
(-α023, +ξ31)

>>> my_alpha = ar('a12 ^ a23')
>>> ar('my_alpha / a123')
-α2
```

### Projecting and adding MultiVectors
```Python
>>> mvec = MultiVector(['1', '2', '3', '123'])
>>> mvec
{
  α123 +ξ123
  α1 +ξ1
  α2 +ξ2
  α3 +ξ3
}

>>> ar('<mvec>1')
{
  α1 +ξ1
  α2 +ξ2
  α3 +ξ3
}

>>> ar('<mvec>2')  # No grade(2) elements so we have an empty MultiVector
{

}

>>> ar('<mvec>3')
{
  α123 +ξ123
}

# The returned MultiVectors from a projection can be added back together to
# recreate the original.
>>> mvec == ar('<mvec>1') + ar('<mvec>3')
True
```


# Differential Operators
By default, arpy provides two differential operators - `Dmu` and `DG` - which
are the operators discussed in John's papers. There is a helper function that
allows you to define your own operators on the fly along with a lower level
`AR_differential` function that allows you to perform one of calculations as you
try out different components.

```Python
>>> mvec = MultiVector(['1', '2', '3', '123'])  # The spatial subgroup: XiA
>>> mvec == XiA
True

>>> Dmu(XiA)
{
  αp (+∂1ξ1, +∂2ξ2, +∂3ξ3)
  α23 (+∂1ξ123, -∂3ξ2, +∂2ξ3)
  α31 (+∂2ξ123, +∂3ξ1, -∂1ξ3)
  α12 (+∂3ξ123, -∂2ξ1, +∂1ξ2)
  α0123 -∂0ξ123
  α10 +∂0ξ1
  α20 +∂0ξ2
  α30 +∂0ξ3
}

>>> DE = differential_operator(['10', '20', '30', '0123'])  # wrt XiE
>>> DE(XiA)
{
  α0 (+∂0123ξ123, -∂10ξ1, -∂20ξ2, -∂30ξ3)
  α023 (+∂10ξ123, -∂0123ξ1, -∂30ξ2, +∂20ξ3)
  α031 (+∂20ξ123, +∂30ξ1, +∂0123ξ2, -∂10ξ3)
  α012 (+∂30ξ123, -∂20ξ1, +∂10ξ2, -∂0123ξ3)
}

# This is equivalent to the following:
>>> AR_differential(XiA, ['10', '20', '30', '0123'])
{
  α0 (+∂0123ξ123, -∂10ξ1, -∂20ξ2, -∂30ξ3)
  α023 (+∂10ξ123, -∂0123ξ1, -∂30ξ2, +∂20ξ3)
  α031 (+∂20ξ123, +∂30ξ1, +∂0123ξ2, -∂10ξ3)
  α012 (+∂30ξ123, -∂20ξ1, +∂10ξ2, -∂0123ξ3)
}
```
