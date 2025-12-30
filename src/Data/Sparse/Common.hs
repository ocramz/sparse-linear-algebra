{-# language TypeFamilies #-}
{-# language TypeOperators, GADTs #-}
{-# language FlexibleInstances, FlexibleContexts #-}
{-# OPTIONS_GHC -Wno-orphans #-}
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
         fromRowsL, fromRowsV, fromColsV, fromColsL, toRowsL, toColsL) where

-- import Control.Exception
-- import Control.Exception.Common
-- import Control.Monad.Catch

import Data.Sparse.Utils as X
import Data.Sparse.PPrint as X
import Data.Sparse.Types as X
import Data.Sparse.Internal.IntMap2 -- as X
import qualified Data.Sparse.Internal.IntM as I
import Data.Sparse.SpMatrix as X
import Data.Sparse.SpVector as X
-- import Data.Sparse.Internal.CSR as X

import Numeric.Eps as X
import Numeric.LinearAlgebra.Class as X

import qualified Data.IntMap.Strict as IM

import GHC.Exts
import Data.Complex

-- import Control.Applicative
-- import Data.Traversable

import Data.Maybe (fromMaybe)
import qualified Data.Vector as V








-- withBoundsSM m ij e f
--   | isValidIxSM m ij = f m ij
--   | otherwise = error e

-- | Modify the size of a SpVector. Do not use directly
resizeSV :: Int -> SpVector a -> SpVector a
resizeSV d2 (SV _ sv) = SV d2 sv
-- | Remap the keys of a SpVector. Do not use directly
mapKeysSV :: (IM.Key -> IM.Key) -> SpVector a -> SpVector a
mapKeysSV fk (SV d sv) = SV d $ I.mapKeys fk sv

-- * Insert row/column vector in matrix

-- | Insert row , using the provided row index transformation function
insertRowWith :: (IxCol -> IxCol) -> SpMatrix a -> SpVector a -> IM.Key -> SpMatrix a
insertRowWith fj (SM (m, n) im) (SV d sv) i
  | not (inBounds0 m i) = error "insertRowWith : index out of bounds"
  | n >= d = SM (m,n) $ I.insert i (insertOrUnion i sv' im) im
  | otherwise = error $ "insertRowWith : incompatible dimensions " ++ show (n, d)
    where sv' = I.mapKeys fj sv
          insertOrUnion i' sv'' im' = maybe sv'' (I.union sv'') (I.lookup i' im')

-- | Insert row          
insertRow :: SpMatrix a -> SpVector a -> IM.Key -> SpMatrix a
insertRow = insertRowWith id

-- | Insert column, using the provided row index transformation function
insertColWith :: (IxRow -> IxRow) -> SpMatrix a -> SpVector a -> IxCol -> SpMatrix a
insertColWith fi smm sv j
  | not (inBounds0 n j) = error "insertColWith : index out of bounds"
  | m >= mv = insIM2 smm vl j
  | otherwise = error $ "insertColWith : incompatible dimensions " ++ show (m,mv) where
      (m, n) = dim smm
      mv = dim sv
      vl = toListSV sv
      insIM2 im2 ((i,x):xs) j' = insIM2 (insertSpMatrix (fi i) j' x im2) xs j'
      insIM2 im2 [] _ = im2

-- | Insert column
insertCol :: SpMatrix a -> SpVector a -> IxCol -> SpMatrix a
insertCol = insertColWith id



-- * Outer vector product

-- | Outer product
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
svToSM (SV n d) = SM (n, 1) $ I.singleton 0 d

-- |Demote (n x 1) or (1 x n) SpMatrix to SpVector
toSV :: SpMatrix a -> SpVector a
toSV (SM (m, n) im) = SV d im' where
  im' | m < n = case toList im of
                  ((_,sv):_) -> sv
                  [] -> I.empty
      | otherwise = fmap extractFirst im
  extractFirst imInner = case toList imInner of
                           ((_,v):_) -> v
                           [] -> error "toSV: empty inner map"
  d | m==1 && n==1 = 1
    | m==1 && n>1 = n 
    | n==1 && m>1 = m
    | otherwise = error $ "toSV : incompatible matrix dimension " ++ show (m,n)  





-- | Lookup a row in a SpMatrix; returns an SpVector with the row, if this is non-empty
lookupRowSM :: SpMatrix a -> IxRow -> Maybe (SpVector a)
lookupRowSM sm i = SV (ncols sm) <$> I.lookup i (dat sm) 


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
  (_, n) = dim mm
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

instance (InnerSpace t, Scalar t ~ t) => LinearVectorSpace (SpVector t) where
  type MatrixType (SpVector t) = SpMatrix t
  (#>) = matVecSD
  (<#) = vecMatSD

matVecSD :: (InnerSpace t, Scalar t ~ t) => SpMatrix t -> SpVector t -> SpVector t
matVecSD (SM (nr, nc) mdata) (SV n sv)
  | nc == n = SV nr $ fmap (`dotu` sv) mdata
  | otherwise = error $ "matVec : mismatched dimensions " ++ show (nc, n)

-- |Vector-on-matrix (FIXME : transposes matrix: more costly than `matVec`, I think)
vecMatSD :: (InnerSpace t, Scalar t ~ t) => SpVector t -> SpMatrix t -> SpVector t
vecMatSD (SV n sv) (SM (nr, nc) mdata)
  | n == nr = SV nc $ fmap (`dotu` sv) (transposeIM2 mdata)
  | otherwise = error $ "vecMat : mismatching dimensions " ++ show (n, nr)

-- | Un-conjugated inner product, to be used in matrix-vector product
dotu :: (Foldable t, Num a, Set t) => t a -> t a -> a
dotu u v = sum $ liftI2 (*) u v 




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






-- | Pack a list of SpVectors as rows of an SpMatrix
fromRowsL :: [SpVector a] -> SpMatrix a
fromRowsL = fromRowsV . V.fromList

-- | Pack a list of SpVectors as columns an SpMatrix
fromColsL :: [SpVector a] -> SpMatrix a
fromColsL = fromColsV . V.fromList

-- | Unpack the rows of an SpMatrix into a list of SpVectors
toRowsL :: SpMatrix a -> [SpVector a]
toRowsL aa = map (extractRow aa) [0 .. m-1] where
  (m,_) = dim aa

-- | Unpack the columns of an SpMatrix into a list of SpVectors
toColsL :: SpMatrix a -> [SpVector a]
toColsL aa = map (extractCol aa) [0 .. n-1] where
  (_,n) = dim aa




-- | Pack a V.Vector of SpVectors as columns of an SpMatrix
fromColsV :: V.Vector (SpVector a) -> SpMatrix a
fromColsV qv = V.ifoldl' ins (zeroSM m n) qv where
  n = V.length qv
  m = dim $ V.head qv
  ins mm i c = insertCol mm c i


-- | Pack a V.Vector of SpVectors as rows of an SpMatrix
fromRowsV :: V.Vector (SpVector a) -> SpMatrix a
fromRowsV qv = V.ifoldl' ins (zeroSM m n) qv where
  m = V.length qv 
  n = dim $ V.head qv
  ins mm i c = insertRow mm c i





-- * Pretty printing


-- showNz :: (Epsilon a, Show a) => a -> String
-- showNz x | nearZero x = " _ "
--          | otherwise = show x

toDenseRow :: Num a => SpMatrix a -> IM.Key -> [a]
toDenseRow sm irow =
  fmap (\icol -> sm @@ (irow,icol)) [0..ncols sm-1]



prdR, prdC :: PPrintOptions
prdR = PPOpts 4 2 7   -- reals
prdC = PPOpts 4 2 16  -- complex values



-- -- printDenseSM :: (Show t, Num t) => SpMatrix t -> IO ()
-- printDenseSM :: (ScIx c ~ (Int, Int), FDSize c ~ (Int, Int), SpContainer c a, Show a, Epsilon a) => c a -> IO ()
printDenseSM :: (PrintDense (SpMatrix a), Show a, Epsilon a) => SpMatrix a -> IO ()
printDenseSM sm = do
  newline
  putStrLn $ sizeStrSM sm
  newline
  prd0 sm



printDenseSM0
  :: (SpMatrix a -> IxRow -> Int -> String) -> SpMatrix a -> IO ()
printDenseSM0 f sm = do
  printDenseSM' sm 5 5
  newline
  where
    printDenseSM' sm' nromax ncomax = mapM_ putStrLn rr_' where
      (nr, _) = (nrows sm, ncols sm)
      rr_ = map (\i -> f sm' i ncomax) [0..nr-1]
      rr_' | nr > nromax = take (nromax - 2) rr_ ++ [" ... "] ++[last rr_]
           | otherwise = rr_ 

printDenseSM0r :: SpMatrix Double -> IO ()
printDenseSM0r sm = printDenseSM0 g sm
  where
    g sm' irow ncolmax = printDN (ncols sm') ncolmax prdR $ toDenseRow sm' irow
printDenseSM0c :: SpMatrix (Complex Double) -> IO ()
printDenseSM0c sm = printDenseSM0 g sm
  where
    g sm' irow ncolmax = printCN (ncols sm') ncolmax prdC $ toDenseRow sm' irow


-- printDenseSV :: (Show t, Epsilon t) => SpVector t -> IO ()
printDenseSV :: PrintDense (SpVector a) => SpVector a -> IO ()
printDenseSV sv = do
  newline
  putStrLn $ sizeStrSV sv
  newline
  prd0 sv

printDenseSV0r :: SpVector Double -> IO ()
printDenseSV0r = printDenseSV0 g where
  g l n = printDN l n prdR
printDenseSV0c :: SpVector (Complex Double) -> IO ()
printDenseSV0c = printDenseSV0 g where
  g l n = printCN l n prdC

printDenseSV0 :: Num a =>
  (Int -> Int -> [a] -> String) -> SpVector a -> IO ()
printDenseSV0 f sv = do
  printDenseSV' (svDim sv) 5
  newline where
    printDenseSV' l n = putStrLn (f l n vd)
    vd = toDenseListSV sv


    
-- ** Pretty printer typeclass

instance PrintDense (SpVector Double) where
  prd = printDenseSV
  prd0 = printDenseSV0r

instance PrintDense (SpVector (Complex Double)) where
  prd = printDenseSV
  prd0 = printDenseSV0c

instance PrintDense (SpMatrix Double) where
  prd = printDenseSM
  prd0 = printDenseSM0r

instance PrintDense (SpMatrix (Complex Double)) where
  prd = printDenseSM
  prd0 = printDenseSM0c





  



-- instance (Elt a, Show a) => PrintDense (CsrMatrix a) where
--   prd = printDenseSM


















  
