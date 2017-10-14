{-# language FlexibleContexts #-}
module Data.Array.Accelerate.Sparse.SMatrix where

import Data.Typeable
import Data.Int (Int32, Int64)


import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Vector, Segments, DIM1, DIM2, Exp, Any(Any), All(All), Z(Z), (:.)((:.)) )

import Data.Array.Accelerate.Sparse.SVector







-- | Sparse matrix, coordinate representation (COO), array-of-struct encoding
data SMatrixCOO e = SMatrixCOO {
    smcooNrows :: Int64
  , smcooNcols :: Int64 
  , smcooEntries :: Array DIM1 (COOElem e)
                           } deriving (Eq, Show)




fromList :: (Show e, Typeable e, A.Elt (COOElem' Int64 e), Foldable t) =>
                  Int64 -> Int64 -> t (Int64, Int64, e) -> SMatrixCOO e
fromList m n ll = SMatrixCOO m n mtx
  where
    mtx = A.fromList dim $ foldr insert [] ll 
    dim = Z :. length ll
    insert (i, j, x) acc = CooE' (i, j, x) : acc 

-- -- identity ::
-- --    (A.IsNum a, A.Elt a) =>
-- --    Exp DIM2 -> Acc (Array DIM2 a)
-- identity sh =
--    A.generate sh
--       (withMatrixIndex $
--        \(_ :. r :. c) -> A.fromIntegral $ A.boolToInt (r A.== c))

-- -- withMatrixIndex :: (A.Lift c a, A.Unlift d b) =>
-- --      (a -> b) -> c (A.Plain a) -> d (A.Plain b)
-- withMatrixIndex f = A.lift . f . A.unlift      





-- | Matrix element in COO form
newtype COOElem' i e = CooE' (i, i, e) deriving (Eq, Show)
type COOElem e = COOElem' Int64 e

-- getRow, getCol :: COOElem i e -> i
getRow (CooE' (i, _, _)) = i
getCol (CooE' (_, j, _)) = j

-- getElem :: COOElem i e -> e
getElem (CooE' (_, _, e)) = e


-- | Lexicographic ordering of matrix elements, rows-first
instance (Ord i, Eq e) => Ord (COOElem' i e) where
  (CooE' (i, j, _)) <= (CooE' (i', j', _))
    | i < i' = True
    | i == i' && j <= j' = True
    | otherwise = False


-- -- No instance for (A.Elt (COOElem' Int64 Double))
-- -- NB: EltRepr is not exported. This seems to mean that users are not supposed to create new array element types. Hm
instance (Show i, Show e, Typeable i, Typeable e) => A.Elt (COOElem' i e) where






-- | CSR


-- | Compressed-sparse row. The definition appearing in the 2011 Accelerate paper. Efficient for matrix-vector action but not for matrix-matrix product
data SMatrixCSR i e = SMatrixCSR {
    smNrows :: Int
  , smNcols :: Int
  , smSegments :: Segments i
  , smEntries :: SVector e } deriving (Eq, Show)
