{-# language FlexibleContexts #-}
module Data.Array.Accelerate.Sparse.SMatrix where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Vector, Segments, DIM1, DIM2, Exp, Any(Any), All(All), Z(Z), (:.)((:.)))

import Data.Array.Accelerate.Sparse.SVector







-- | Sparse matrix, coordinate representation (COO), array-of-struct encoding
data SMatrixCOO i e = SMatrixCOO {
    smcooNrows :: Int
  , smcooNcols :: Int   
  , smcooEntries :: Array DIM1 (COOElem i e)
                           } deriving (Eq, Show)





fromList :: (A.Elt (COOElem i e), Foldable t) =>
                  Int -> Int -> t (i, i, e) -> SMatrixCOO i e
fromList m n ll = SMatrixCOO m n mtx
  where
    mtx = A.fromList dim $ foldr insert [] ll 
    dim = Z :. length ll
    insert (i, j, x) acc = CooE (i, j, x) : acc 

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
newtype COOElem i e = CooE (i, i, e) deriving (Eq, Show)

getRow, getCol :: COOElem i e -> i
getRow (CooE (i, _, _)) = i
getCol (CooE (_, j, _)) = j

getElem :: COOElem i e -> e
getElem (CooE (_, _, e)) = e


-- | Lexicographic ordering of matrix elements, rows-first
instance (Eq e, Ord i) => Ord (COOElem i e ) where
  (CooE (i, j, _)) <= (CooE (i', j', _))
    | i < i' = True
    | i == i' && j <= j' = True
    | otherwise = False



-- | CSR


-- | Compressed-sparse row. The definition appearing in the 2011 Accelerate paper. Efficient for matrix-vector action but not for matrix-matrix product
data SMatrixCSR i e = SMatrixCSR {
    smNrows :: Int
  , smNcols :: Int
  , smSegments :: Segments i
  , smEntries :: SVector e } deriving (Eq, Show)
