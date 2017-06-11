{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses #-}
module Data.Sparse.Internal.SpMatrix_Lex where

import qualified Data.Vector as V 

import Data.Sparse.Internal.Utils
import Data.Sparse.Types
import Numeric.LinearAlgebra.Class

-- * SpMatrix1 : entry coordinates are stored in lexicographic encoding

data SpMatrix1 a = SM1 { smNrows :: Rows,
                         smNcols :: Cols,
                         smLexOrd :: LexOrd,
                         smPtr :: V.Vector Int,
                         smDat :: V.Vector (LexIx, a) } deriving (Eq, Show)

instance Functor SpMatrix1 where
  fmap f (SM1 r c et p x) = SM1 r c et p (fmap g x) where
    g (i, v) = (i, f v)

instance FiniteDim SpMatrix1 where
  type FDSize SpMatrix1 = (Rows, Cols)
  dim x = (smNrows x, smNcols x)

instance HasData SpMatrix1 a where
  type HDData SpMatrix1 a = V.Vector (LexIx, a)
  dat = smDat
  nnz = V.length . dat

instance Sparse SpMatrix1 a where
  spy mm = fromIntegral (nnz mm) / fromIntegral (m * n) where (m, n) = dim mm

-- instance Num a => SpContainer SpMatrix1 a where
--   type ScIx SpMatrix1 = LexIx
--   scToList = V.toList . fmap snd . dat
--   scLookup = smLookup
--   mm @@ i = fromMaybe 0 (scLookup mm i)

-- smLookup :: SpMatrix1 a -> LexIx -> Maybe a
-- smLookup mm i = snd <$> V.find ((== i) . fst) (dat mm)






-- ** Encoding/decoding between coordinate representation and lexicographic index
-- | Coordinates to lexicographic 
encode :: LexOrd -> Rows -> Cols -> (IxRow, IxCol) -> LexIx
encode t 
  | t == RowsFirst = encodeR
  | otherwise = encodeC where
      encodeR m _ (i, j) = m * j + i
      encodeC _ n (i, j) = n * i + j

-- | Lexicographic  to coordinates
decode :: LexOrd -> Rows -> Cols -> LexIx -> (IxRow, IxCol)
decode t
  | t == RowsFirst = decodeR
  | otherwise = decodeC where
      decodeR m _ ii = (i, j) where
        (j, i) = divMod ii m
      decodeC _ n jj = divMod jj n


-- instance Num a => SparseMatrix SpMatrix1 a where
--   smFromVector et (m, n) v = SM1 m n et p' v' where
--     enc (i, j, x) = (encode et m n (i, j), x) -- encode indices
--     v' = sortByIx $ V.map enc v
--     p' = csrPtrVTup et m n v'
--   encodeIx mm o = encode o m n where (m, n) = dim mm
--   decodeIx mm o = decode o m n where (m, n) = dim mm
-- --   smTranspose ColsFirst = transposeSM1


-- | Generate the row/column row pointer according to the respective encoding
csrPtrVTup ::
  LexOrd -> Rows -> Cols -> V.Vector (LexIx, t) -> V.Vector Int
csrPtrVTup et m n v
  | et == RowsFirst = csPtrV (\(i, _) i0 -> fst (dec i) == i0) m v
  | otherwise       = csPtrV (\(i, _) i0 -> snd (dec i) == i0) n v 
    where
      dec = decode et m n
      




-- NB : binary operations on SpMatrix1 are defined only if dimensions AND lexical order encoding are equal
-- instance Num a => AdditiveGroup (SpMatrix1 a) where
--   -- zero = SM1 0 0 ColsFirst V.empty V.empty
--   negated = fmap negate

-- stream merge -> SpMatrix1 binary lift





-- -- | Transpose: remap keys ( O(N) ) followed by merge sort ( O(log N) )
-- transposeSM1 :: Num a => SpMatrix1 a -> SpMatrix1 a
-- transposeSM1 mm = SM1 n m p' d' where
--   fi (i, j) = (j, i)
--   (m, n) = dim mm
--   d' = sortByIx $ fmapIx (remapCoords fi mm) (dat mm)
--   p' = csrPtrV' (\(i, x) i0 -> i==i0) n d'

-- remapCoords ::
--   SparseMatrix m e => ((IxRow, IxCol) -> (IxRow, IxCol)) -> m e -> LexIx -> LexIx
-- remapCoords fi mm = encodeIx mm ColsFirst . fi . decodeIx mm ColsFirst


-- -- | Map a function over the indices (internal use only)
-- fmapIx :: Functor f => (a -> b) -> f (a, t) -> f (b, t)
-- fmapIx f tt = fmap g tt where g (i, x) = (f i, x)

-- -- | Sort according to index 
-- sortByIx :: V.Vector (LexIx, b) -> V.Vector (LexIx, b)
-- sortByIx = V.modify (VA.sortBy f) where
--   f x y = compare (fst x) (fst y)
