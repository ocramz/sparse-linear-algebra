
	0.2.9.7
	------

	Improved pretty printer:

	* Fixed display precision (e.g. 2 decimal digits), fixed column width output for vectors and matrices
	
	* Small and large values (wrt fixed precision) switch to scientific notation
	
	0.2.9.4
	-------

	Exceptions constructors are exported by Numeric.LinearAlgebra.Sparse


	0.2.9.1
	-------

	* Uses classes from `vector-space` such as AdditiveGroup, VectorSpace and InnerSpace
	* QuickCheck tests for algebraic properties, such as matrix-vector products and soon more abstract ones e.g. positive semi-definite matrices
	
	* Getting rid of `error` in favor of MonadThrow exceptions for high-level operations such as matrix algorithms
