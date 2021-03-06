	0.3.2
	* Introduced `logging-effect` as a more versatile alternative to debug logging in IO

	0.3.1
	* Changed `SpMatrix` to use `foldlWithKey'` for efficiency (Joshua Moerman)
	
	* Bumped LTS to 11.3 (GHC 8.2.2)
	
	* Removed unneeded dependencies from stack.yaml

	0.3
	* Fixed a number of instances, un-commented tests (Joshua Moerman)
	
	* Documented issues with complex number support (Joshua Moerman)
	
	0.2.9.9
	* Moved to IntMap.Strict (Gregory Schwartz)
	
	* Stackage LTS bump to 10.4 (GHC 8.2)
	
	0.2.9.7
	Improved pretty printer:

	* Fixed display precision (e.g. 2 decimal digits), fixed column width output for vectors and matrices
	
	* Small and large values (wrt fixed precision) switch to scientific notation
	
	0.2.9.4
	Exceptions constructors are exported by Numeric.LinearAlgebra.Sparse


	0.2.9.1
	* Uses classes from `vector-space` such as AdditiveGroup, VectorSpace and InnerSpace
	
	* QuickCheck tests for algebraic properties, such as matrix-vector products and soon more abstract ones e.g. positive semi-definite matrices
	
	* Getting rid of `error` in favor of MonadThrow exceptions for high-level operations such as matrix algorithms
