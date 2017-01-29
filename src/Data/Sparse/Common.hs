{-# language TypeFamilies #-}
{-# language TypeOperators, GADTs #-}
{-# language FlexibleInstances, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------
module Data.Sparse.Common
       ( module X,
         insertRowWith, insertRow, insertColWith, insertCol,
         diagonalSM,
         outerProdSV, (><), toSV, svToSM,
         lookupRowSM, 
         extractCol, extractRow,
         extractVectorDenseWith, extractRowDense, extractColDense,
         extractDiagDense,
         extractSubRow, extractSubCol,
         extractSubRow_RK, extractSubCol_RK,
         fromCols) where

import Control.Exception
import Control.Exception.Common
import Control.Monad.Catch

import Data.Sparse.Utils as X
import Data.Sparse.PPrint as X
import Data.Sparse.Types as X
import Data.Sparse.Internal.IntMap2 -- as X
import Data.Sparse.SpMatrix as X
import Data.Sparse.SpVector as X
-- import Data.Sparse.Internal.CSR as X

import Numeric.Eps as X
import Numeric.LinearAlgebra.Class as X

import qualified Data.IntMap as IM

import Data.Complex

-- import Control.Applicative
-- import Data.Traversable

import Data.Maybe (fromMaybe, maybe)
import qualified Data.Vector as V

import Data.VectorSpace

-- withBoundsSM m ij e f
--   | isValidIxSM m ij = f m ij
--   | otherwise = error e

-- | modify the size of a SpVector. Do not use directly
resizeSV :: Int -> SpVector a -> SpVector a
resizeSV d2 (SV _ sv) = SV d2 sv


mapKeysSV fk (SV d sv) = SV d $ IM.mapKeys fk sv

-- * Insert row/column vector in matrix

-- | Insert row , using the provided row index transformation function
insertRowWith :: (IxCol -> IxCol) -> SpMatrix a -> SpVector a -> IM.Key -> SpMatrix a
insertRowWith fj (SM (m,n) im) (SV d sv) i
  | not (inBounds0 m i) = error "insertRowSM : index out of bounds"
  | n >= d = SM (m,n) $ IM.insert i (insertOrUnion i sv' im) im
  | otherwise = error $ "insertRowSM : incompatible dimensions " ++ show (n, d)
    where sv' = IM.mapKeys fj sv
          insertOrUnion i sv im = maybe sv (IM.union sv) (IM.lookup i im)

-- | Insert row          
insertRow :: SpMatrix a -> SpVector a -> IM.Key -> SpMatrix a
insertRow = insertRowWith id

-- | Insert column, using the provided row index transformation function
insertColWith :: (IxRow -> IxRow) -> SpMatrix a -> SpVector a -> IxCol -> SpMatrix a
insertColWith fi smm sv j
  | not (inBounds0 n j) = error "insertColSM : index out of bounds"
  | m >= mv = insIM2 smm vl j
  | otherwise = error $ "insertColSM : incompatible dimensions " ++ show (m,mv) where
      (m, n) = dim smm
      mv = dim sv
      vl = toListSV sv
      insIM2 im2 ((i,x):xs) j = insIM2 (insertSpMatrix (fi i) j x im2) xs j
      insIM2 im2 [] _ = im2

-- | Insert column
insertCol :: SpMatrix a -> SpVector a -> IxCol -> SpMatrix a
insertCol = insertColWith id



-- * Outer vector product

-- | Outer product (all-with-all matrix)
outerProdSV, (><) :: Num a => SpVector a -> SpVector a -> SpMatrix a
outerProdSV v1 v2 = fromListSM (m, n) ixy where
  m = dim v1
  n = dim v2
  ixy = [(i,j, x * y) | (i,x) <- toListSV v1 , (j, y) <- toListSV v2]

(><) = outerProdSV



-- * Diagonal matrix

-- | Fill the diagonal of a SpMatrix with the components of a SpVector
diagonalSM :: SpVector a -> SpMatrix a
diagonalSM sv = ifoldSV iins (zeroSM n n) sv where
  n = dim sv
  iins i = insertSpMatrix i i



-- * Matrix-vector conversions

-- | promote a SV to SM
svToSM :: SpVector a -> SpMatrix a
svToSM (SV n d) = SM (n, 1) $ IM.singleton 0 d

-- |Demote (n x 1) or (1 x n) SpMatrix to SpVector
toSV :: SpMatrix a -> SpVector a
toSV (SM (m,n) im) = SV d (ff im) where
  ff | m > n = IM.map g  -- column case
     | otherwise = g
  g = snd . head . IM.toList
  d | m==1 && n==1 = 1
    | m==1 && n>1 = n 
    | n==1 && m>1 = m
    | otherwise = error $ "toSV : incompatible matrix dimension " ++ show (m,n)



-- | Lookup a row in a SpMatrix; returns an SpVector with the row, if this is non-empty
lookupRowSM :: SpMatrix a -> IxRow -> Maybe (SpVector a)
lookupRowSM sm i = SV (ncols sm) <$> IM.lookup i (dat sm) 


-- * Extract a SpVector from an SpMatrix
-- ** Sparse extract

-- |Extract ith row
extractRow :: SpMatrix a -> IxRow -> SpVector a
extractRow m i
  | inBounds0 (nrows m) i = fromMaybe (zeroSV (ncols m)) (lookupRowSM m i)
  | otherwise = error $ unwords ["extractRow : index",show i,"out of bounds"]

-- |Extract jth column
extractCol :: SpMatrix a -> IxCol -> SpVector a
extractCol m j = toSV $ extractColSM m j      



-- ** Dense extract (default == 0)
-- | Generic extraction function
extractVectorDenseWith ::
  Num a => (Int -> (IxRow, IxCol)) -> SpMatrix a -> SpVector a
extractVectorDenseWith f mm = fromListDenseSV n $ foldr ins [] ll  where
  ll = [0 .. n - 1]
  (m, n) = dim mm
  ins i acc = mm @@ f i : acc

-- | Extract ith row (dense)
extractRowDense :: Num a => SpMatrix a -> IxRow -> SpVector a
extractRowDense mm iref = extractVectorDenseWith (\j -> (iref, j)) mm

-- | Extract jth column
extractColDense :: Num a => SpMatrix a -> IxCol -> SpVector a
extractColDense mm jref = extractVectorDenseWith (\i -> (i, jref)) mm

-- | Extract the diagonal
extractDiagDense :: Num a => SpMatrix a -> SpVector a
extractDiagDense = extractVectorDenseWith (\i -> (i, i))



-- | extract row interval (all entries between columns j1 and j2, INCLUDED, are returned)
-- extractSubRow :: SpMatrix a -> IxRow -> (IxCol, IxCol) -> SpVector a
-- extractSubRow m i (j1, j2) = case lookupRowSM m i of
--   Nothing -> zeroSV (ncols m)
--   Just rv -> ifilterSV (\j _ -> j >= j1 && j <= j2) rv

-- |", returning in Maybe
-- extractSubRow :: SpMatrix a -> IxRow -> (Int, Int) -> Maybe (SpVector a)
-- extractSubRow m i (j1, j2) =
--   resizeSV (j2 - j1) . ifilterSV (\j _ -> j >= j1 && j <= j2) <$> lookupRowSM m i

-- | Extract an interval of SpVector components, changing accordingly the resulting SpVector size. Keys are _not_ rebalanced, i.e. components are still labeled according with respect to the source matrix.
extractSubRow :: SpMatrix a -> IxRow -> (Int, Int) -> SpVector a
extractSubRow m i (j1, j2) = fromMaybe (zeroSV deltaj) vfilt where
  deltaj = j2 - j1 + 1
  vfilt = resizeSV deltaj .
          ifilterSV (\j _ -> j >= j1 && j <= j2) <$> lookupRowSM m i

-- | extract row interval, rebalance keys by subtracting lowest one
extractSubRow_RK :: SpMatrix a -> IxRow -> (IxCol, IxCol) -> SpVector a
extractSubRow_RK m i (j1, j2) = mapKeysSV (subtract j1) $ extractSubRow m i (j1, j2)

  -- toSV $ extractSubRowSM_RK m i (j1, j2)



-- | extract column interval
extractSubCol :: SpMatrix a -> IxCol -> (IxRow, IxRow) -> SpVector a
extractSubCol m j (i1, i2)  = toSV $ extractSubColSM m j (i1, i2)

-- | extract column interval, rebalance keys by subtracting lowest one
extractSubCol_RK :: SpMatrix a -> IxCol -> (IxRow, IxRow) -> SpVector a
extractSubCol_RK m j (i1, i2)  = toSV $ extractSubColSM_RK m j (i1, i2)




-- ** Matrix action on a vector

{- 
FIXME : matVec is more general than SpVector's :

\m v -> fmap (`dot` v) m
  :: (Normed f1, Num b, Functor f) => f (f1 b) -> f1 b -> f b
-}





-- | LinearVectorSpace instances
instance LinearVectorSpace (SpVector Double) where
  type MatrixType (SpVector Double) = SpMatrix Double
  (#>) = matVecSD
  (<#) = vecMatSD

instance LinearVectorSpace (SpVector (Complex Double)) where
  type MatrixType (SpVector (Complex Double)) = SpMatrix (Complex Double)
  (#>) = matVecSD
  (<#) = vecMatSD


matVecSD :: (InnerSpace (IM.IntMap t), s ~ Scalar (IM.IntMap t)) =>
    SpMatrix t -> SpVector t -> SpVector s
matVecSD (SM (nr, nc) mdata) (SV n sv)
  | nc == n = SV nr $ fmap (`dot` sv) mdata
  | otherwise = error $ "matVec : mismatched dimensions " ++ show (nc, n)

-- |Vector-on-matrix (FIXME : transposes matrix: more costly than `matVec`, I think)
-- vecMat, (<#) :: SpVector Double -> SpMatrix Double -> SpVector Double
vecMatSD :: (InnerSpace (IM.IntMap t), s ~ Scalar(IM.IntMap t)) =>
    SpVector t -> SpMatrix t -> SpVector s
vecMatSD (SV n sv) (SM (nr, nc) mdata)
  | n == nr = SV nc $ fmap (`dot` sv) (transposeIM2 mdata)
  | otherwise = error $ "vecMat : mismatching dimensions " ++ show (n, nr)





-- -- generalized matVec : we require a function `rowsf` that produces a functor of elements of a Hilbert space (the rows of `m`)
-- matVecG :: (Hilbert v, Functor f, f (Scalar v) ~ v) => (m -> f v) -> m -> v -> v
-- matVecG rowsf m v = fmap (`dot` v) (rowsf m)

-- matVecGA
--   :: (Hilbert v, Traversable t, t (Scalar v) ~ v) =>
--      (m -> t v) -> m -> v -> v
-- matVecGA rowsf m v = traverse (<.> v) (rowsf m)

-- -- --  Really, a matrix is just notation for a linear map between two finite-dimensional Hilbert spaces, i.e.
-- matVec :: (Hilbert u, Hilbert v) => (u -> v) -> u -> v
-- which is a specialization of a function application operator like ($) :: (a -> b) -> a -> b




-- -- --  from `vector-space`

-- data a -* b where
--   Dot :: VectorSpace b => b -> (b -* Scalar b)
--   (:&&) :: (a -* c) -> (a -* d) -> (a -* (c, d)) -- a,c,d should be constrained

-- apply :: Hilbert a => (a -* b) -> (a -> b)
-- apply (Dot b) = dot b
-- apply (f :&& g) = apply f &&& apply g
--   where (u &&& v) a = (u a, v a) -- (&&&) from Control.Arrow

-- -- type a :~ b = Scalar a ~ Scalar b


















-- | Pack a V.Vector of SpVectors as columns of an SpMatrix
fromCols :: V.Vector (SpVector a) -> SpMatrix a
fromCols qv = V.ifoldl' ins (zeroSM m n) qv where
  n = V.length qv
  m = dim $ V.head qv
  ins mm i c = insertCol mm c i









-- * Pretty printing


showNz :: (Epsilon a, Show a) => a -> String
showNz x | nearZero x = " _ "
        | otherwise = show x

-- toDenseRow :: Num a => SpMatrix a -> IM.Key -> [a]
toDenseRow sm irow =
  fmap (\icol -> sm @@ (irow,icol)) [0..ncol-1] where (_, ncol) = dim sm



-- toDenseRowClip :: (Show a, Num a) => SpMatrix a -> IM.Key -> Int -> String
toDenseRowClip sm irow ncomax
  | nco > ncomax = unwords (map showNz h) ++  " ... " ++ showNz t
  | otherwise = unwords $ showNz <$> dr
     where dr = toDenseRow sm irow
           h = take (ncomax - 2) dr
           t = last dr
           (_, nco) = dim sm


-- printDenseSM :: (Show t, Num t) => SpMatrix t -> IO ()
printDenseSM :: (ScIx c ~ (Int, Int), FDSize c ~ (Int, Int), SpContainer c a, Show a, Epsilon a) => c a -> IO ()
printDenseSM sm = do
  newline
  putStrLn $ sizeStr sm
  newline
  printDenseSM' sm 5 5
  newline
  where    
    -- printDenseSM' :: (Show t, Num t) => SpMatrix t -> Int -> Int -> IO ()
    printDenseSM' sm' nromax ncomax = mapM_ putStrLn rr_' where
      (nr, _) = dim sm'
      rr_ = map (\i -> toDenseRowClip sm' i ncomax) [0..nr - 1]
      rr_' | nr > nromax = take (nromax - 2) rr_ ++ [" ... "] ++[last rr_]
           | otherwise = rr_



toDenseListClip :: (Show a, Num a, Epsilon a) => SpVector a -> Int -> String
toDenseListClip sv ncomax
  | dim sv > ncomax = unwords (map showNz h) ++  " ... " ++ showNz t
  | otherwise = unwords $ showNz <$> dr
     where dr = toDenseListSV sv
           h = take (ncomax - 2) dr
           t = last dr

printDenseSV :: (Show t, Num t, Epsilon t) => SpVector t -> IO ()
printDenseSV sv = do
  newline
  putStrLn $ sizeStrSV sv
  newline
  printDenseSV' sv 5
  newline where
    printDenseSV' v nco = putStrLn rr_' where
      rr_ = toDenseListClip v nco :: String
      rr_' | dim sv > nco = unwords [take (nco - 2) rr_ , " ... " , [last rr_]]
           | otherwise = rr_

-- ** Pretty printer typeclass


instance (Show a, Num a, Epsilon a) => PrintDense (SpVector a) where
  prd = printDenseSV

instance (Show a, Num a, Epsilon a) => PrintDense (SpMatrix a) where
  prd = printDenseSM



  



-- instance (Elt a, Show a) => PrintDense (CsrMatrix a) where
--   prd = printDenseSM


















  
