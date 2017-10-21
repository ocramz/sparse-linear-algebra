{-# language FlexibleContexts, TypeFamilies #-}
{-# language ScopedTypeVariables #-}
{-# language UndecidableInstances #-}
module Data.Array.Accelerate.Sparse.SMatrix where

import Data.Typeable
import Data.Int (Int32, Int64)

-- import Data.Array.Accelerate.Lift
import Data.Array.Accelerate.Product
import Data.Array.Accelerate.Smart
import Data.Array.Accelerate.Type
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Array.Data

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Vector, Segments, DIM1, DIM2, Exp, Any(Any), All(All), Z(Z), (:.)((:.)) )

import Data.Array.Accelerate.Sparse.COOElem
import Data.Array.Accelerate.Sparse.SVector


-- | Sparse matrix, coordinate representation (COO), array-of-struct encoding
data SMatrixCOO i e = SMatrixCOO {
    smcooNrows :: i
  , smcooNcols :: i
  , smcooEntries :: Array DIM1 (COOElem i e)
                           } deriving (Eq, Show)


fromListCOO :: (Elt e, Elt i, Foldable t) => i -> i -> t (i, i, e) -> SMatrixCOO i e
fromListCOO m n ll = SMatrixCOO m n mtx
  where
    mtx = A.fromList d $ foldr ins [] ll 
    d = Z :. length ll
    ins (i, j, x) acc = CooE (i, j, x) : acc 

toListCOO :: SMatrixCOO i e -> [COOElem i e]
toListCOO (SMatrixCOO _ _ ll) = A.toList ll





-- | CSR


-- | Compressed-sparse row. The definition appearing in the 2011 Accelerate paper. Efficient for matrix-vector action but not for matrix-matrix product
data SMatrixCSR i e = SMatrixCSR {
    smNrows :: Int
  , smNcols :: Int
  , smSegments :: Segments i
  , smEntries :: SVector e } deriving (Eq, Show)



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

