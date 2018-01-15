module Data.Sparse.Internal.Vector.Utils where

import qualified Data.Vector as V 
import qualified Data.Vector.Mutable as VM

import Control.Monad.ST

-- | Given a number of rows(resp. columns) `n` and a _sorted_ Vector of Integers in increasing order (containing the row(col) indices of nonzero entries), return the cumulative vector of nonzero entries of length `n + 1` (the "row(col) pointer" of the CSR(CSC) format). NB: Fused count-and-accumulate
-- E.g.:
-- > csPtrV (==) 4 (V.fromList [1,1,2,3])
-- [0,0,2,3,4]
csPtrV :: (a -> Int -> Bool) -> Int -> V.Vector a -> V.Vector Int
csPtrV eqf n xs = V.create createf where
  createf :: ST s (VM.MVector s Int)
  createf = do
    let c = 0
    vm <- VM.new (n + 1)
    VM.write vm 0 0  -- write `0` at position 0
    let loop v ll i count | i == n = return ()
                          | otherwise = do
                              let lp = V.length $ V.takeWhile (`eqf` i) ll
                                  count' = count + lp
                              VM.write v (i + 1) count'
                              loop v (V.drop lp ll) (succ i) count'
    loop vm xs 0 c
    return vm
