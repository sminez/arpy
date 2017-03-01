Available Classes and their APIs
================================

# Alpha
This is the sole representation of unit elements (αμνρ) within the algebra. All
other classes (other than the Xi elements) and functions require _Alpha_
instances in order to function in accordance with the principle of absolute
relativity.
### Constructor
A string representing the desired α value: `1`, `-012`, `p` etc.
It is also possible to pass the sign as a second argument of +- 1 but this is
primarily intended for programatic use within the library itself.
### Methods
None! Use of Alphas is via the library functions and `ar()` interface.

# Xi
A real value.
### Constructor
Any python object. If an Alpha is passed then its index is extracted to create a
symbolic Xi value denoted `ξ<index>`.
### Methods
None.

# Pair
A tuple of an `Alpha` and a `Xi`.
### Constructor
As with the Xi class, any Python object can be passed to the constructor as the
second argument (the first must be an Alpha). If a Xi is passed then it is used
as is; otherwise the second argument will be implicitly converted to a Xi
object.
### Methods
None.

# MultiVector
An arbitrary linear combination of Pairs.
### Constructor
Any Python collection type that supports the iterator proctocol can be passed so
long as its components are Pairs, Alphas or strings.
In the case of Alphas and strings, the components will be converted to symbolic
Pairs as detailed above before storage.
### Methods
Multivectors support iteration (all Pairs returned in MTAE order), Dict style
lookup and the `.get()` method where the argument is either a string or an
Alpha. The return is a list of all Pair that with that Alpha and a linear sum of
all components in the MultiVector with that Alpha value.

It is also possible to add two MultiVectors which will create a new MultiVector
from the linear sum of all components of each MultiVector.

`BTAE_grouped()` will pretty print a version of the MultiVector that has been
collected into the BTAE α notation.

`collected_terms()` will attempt to factorise the Xi values and pretty print the
output.
