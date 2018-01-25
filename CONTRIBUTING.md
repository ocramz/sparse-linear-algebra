# Contribution guidelines

Chatroom: https://gitter.im/ocramz/sparse-linear-algebra
Issues: https://github.com/ocramz/sparse-linear-algebra/issues

`sparse-linear-algebra` is a young and experimental library, but it aims to provide high-performance _and_ high-level functionality.

A combination of approaches is necessary to meet these ambitious goals:

* Performance
* Ergonomics
* Applications

## Performance

This is the trickiest bit (especially in Haskell). The initial version was definitely not performant, which motivated a
complete rewrite and reorganization of the backend (which will be split into `sparse-linear-algebra-accelerate` and possibly also 
`sparse-linear-algebra-vector`).
Got an idea for code generation of a certain method? Would you like to have a `repa` backend or perhaps have something to say
about how a method is implemented? Don't esitate to join the discussion.

## Ergonomics

One of the first aims of this library was to bridge the gap between scientific computing textbooks and application code.
We often think in math, but must talk to computers; the translation doesn't need to be a traumatic experience!
To all prospective users and contributors: do try the library out! Is there something missing, or does something feel wrong or clunky? 
Please get in touch and share all your thoughts !

## Applications/Libraries

This library aims to being used in the "real world", so all feedback from developers of applications and wrapping libraries is more than welcome.
Anything from performance regressions, or un-exposed functionality such as exceptions etc.
Shoot us a line, or even better a pull request!
