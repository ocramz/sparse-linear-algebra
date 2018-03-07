{-# language FlexibleInstances, FlexibleContexts, TypeFamilies #-}
{-# language DeriveFunctor, DeriveFoldable #-}
-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------
module Data.Sparse.SpMatrix where

-- import Control.Exception.Common
-- import Data.Sparse.SpVector
import Data.Sparse.Utils
import Data.Sparse.Types

import Numeric.Eps
import Numeric.LinearAlgebra.Class

import Data.Sparse.Internal.IntM (IntM (..))
import qualified Data.Sparse.Internal.IntM as I
import Data.Sparse.Internal.IntMap2

import GHC.Exts
import Text.Printf

import qualified Data.IntMap.Strict as IM

import Data.Complex
import Data.Foldable (foldl')
import Data.Maybe

-- import Data.VectorSpace hiding (magnitude)


-- *

-- instance IxContainer SpMatrix a where
--   type Ix SpMatrix = (Int, Int)
--   type IxSz SpMatrix = (Int, Int)
--   ixcLookup m (i,j) = lookupSM m i j
--   ixcIfilter g = ifilterSM g' where g' i j = g (i, j)
--   ixcInsert (i, j) = insertSpMatrix i j
--   -- ixcFromList = fromListSM

-- instance SMatrix SpMatrix a where



-- * Sparse Matrix

data SpMatrix a = SM {smDim :: {-# UNPACK #-} !(Rows, Cols),
                      smData :: !(IntM (IntM a))}
                deriving (Eq, Functor, Foldable)


-- sizeStr :: (FDSize f ~ (a1, a2), Sparse f a, Show a2, Show a1) => f a -> String
sizeStrSM :: SpMatrix a -> String
sizeStrSM sm@(SM (nr, nc) _) =
  unwords ["(",show nr,"rows,",show nc,"columns ) ,",show nz,"NZ ( density",sys,")"] where
  -- (nr, nc) = dim sm
  nz = nnz sm
  sy = spy sm :: Double
  sys = printf "%1.3f %%" (sy * 100) :: String


instance Show a => Show (SpMatrix a) where
  show sm@(SM _ x) = unwords ["SM",sizeStrSM sm,show (toList $ toList <$> x)]
  -- show sm@(SM _ x) = show x

instance Set SpMatrix where
  liftU2 f2 (SM n1 x1) (SM n2 x2) = SM (maxTup n1 n2) ((liftU2.liftU2) f2 x1 x2)
  liftI2 f2 (SM n1 x1) (SM n2 x2) = SM (minTup n1 n2) ((liftI2.liftI2) f2 x1 x2)

-- | 'SpMatrix'es form an additive group, in that they can have an invertible associtative operation (matrix sum)
instance AdditiveGroup a => AdditiveGroup (SpMatrix a) where
  zeroV = SM (0,0) I.empty
  (^+^) = liftU2 (^+^)
  negateV = fmap negateV

instance VectorSpace a => VectorSpace (SpMatrix a) where
  type Scalar (SpMatrix a) = Scalar a
  n .* v = fmap (n .*) v




-- | 'SpMatrix'es are maps between finite-dimensional spaces
instance FiniteDim (SpMatrix a) where
  type FDSize (SpMatrix a) = (Rows, Cols)
  dim = smDim

instance HasData (SpMatrix a) where
  type HDData (SpMatrix a) = IntM (IntM a)
  nnz = nzSM
  dat = smData

instance Sparse (SpMatrix a) where
  spy = spySM

-- | 'SpMatrix'es are sparse containers too, i.e. any specific component may be missing (so it is assumed to be 0)
instance Num a => SpContainer (SpMatrix a) where
  type ScIx (SpMatrix a) = (Rows, Cols)
  type ScElem (SpMatrix a) = a
  scInsert (i,j) = insertSpMatrix i j
  scLookup m (i, j) = lookupSM m i j
  m @@ d | isValidIxSM m d = m @@! d
         | otherwise = error $ "@@ : incompatible indices : matrix size is " ++ show (dim m) ++ ", but user looked up " ++ show d









-- ** Creation

-- | `zeroSM m n` : Empty SpMatrix of size (m, n)
zeroSM :: Rows -> Cols -> SpMatrix a
zeroSM m n = SM (m,n) I.empty


-- *** Diagonal matrix
-- | `mkDiagonal n ll` : create a diagonal matrix of size `n` from a list `ll` of elements
mkDiagonal :: Int -> [a] -> SpMatrix a
mkDiagonal n = mkSubDiagonal n 0

-- *** Identity matrix
-- | `eye n` : identity matrix of rank `n`
eye :: Num a => Int -> SpMatrix a
eye n = mkDiagonal n (replicate n 1)


-- *** Permutation matrix

-- | Permutation matrix from a (possibly incomplete) list of row swaps starting from row 0
-- e.g. `permutationSM 5 [1,3]` first swaps rows (0, 1) and then rows (1, 3) :
-- 
-- >>> prd (permutationSM 5 [1,3] :: SpMatrix Double)
-- 
-- @
-- ( 5 rows, 5 columns ) , 5 NZ ( density 20.000 % )
-- 
-- _      , 1.00   , _      , _      , _      
-- _      , _      , _      , 1.00   , _      
-- _      , _      , 1.00   , _      , _      
-- 1.00   , _      , _      , _      , _      
-- _      , _      , _      , _      , 1.00
-- @
permutationSM :: Num a => Int -> [IxRow] -> SpMatrix a
permutationSM n iis = permutPairsSM n (zip [0 .. n-1] iis)

-- | Permutation matrix from a (possibly incomplete) list of row pair swaps
-- e.g. `permutPairs 5 [(2,4)]` swaps rows 2 and 4 :
-- 
-- >>> prd (permutPairsSM 5 [(2,4)] :: SpMatrix Double)
-- 
-- @
-- ( 5 rows, 5 columns ) , 5 NZ ( density 20.000 % )
-- 
-- 1.00   , _      , _      , _      , _      
-- _      , 1.00   , _      , _      , _      
-- _      , _      , _      , _      , 1.00   
-- _      , _      , _      , 1.00   , _      
-- _      , _      , 1.00   , _      , _
-- @
permutPairsSM :: Num a => Int -> [(IxRow, IxRow)] -> SpMatrix a
permutPairsSM n iix = go iix (eye n) where
  go ((i1, i2):iis) m = go iis (swapRows i1 i2 m)
  go [] m = m







-- *** Super- or sub- diagonal matrix
-- | `mkSubDiagonal n o xx` creates a square SpMatrix of size `n` with `xx` on the `o`th subdiagonal
mkSubDiagonal :: Int -> Int -> [a] -> SpMatrix a
mkSubDiagonal n o xx | abs o < n = if o >= 0
                                   then fz ii jj xx
                                   else fz jj ii xx
                     | otherwise = error "mkSubDiagonal : offset > dimension" where
  ii = [0 .. n-1]
  jj = [abs o .. n - 1]
  fz a b x = fromListSM (n,n) (zip3 a b x)


-- fromList :: [(Key,a)] -> IntMap a
-- fromList xs
--   = foldlStrict ins empty xs
--   where
--     ins t (k,x)  = insert k x t



-- ** Element insertion

-- | Insert an element in a preexisting Spmatrix at the specified indices
insertSpMatrix :: IxRow -> IxCol -> a -> SpMatrix a -> SpMatrix a
insertSpMatrix i j x s
  | inBounds02 d (i,j) = SM d $ insertIM2 i j x smd 
  | otherwise = error "insertSpMatrix : index out of bounds" where
      smd = immSM s
      d = dim s




-- ** fromList

-- | Add to existing SpMatrix using data from list (row, col, value)
fromListSM' :: Foldable t => t (IxRow, IxCol, a) -> SpMatrix a -> SpMatrix a
fromListSM' iix sm = foldl' ins sm iix where
  ins t (i,j,x) = insertSpMatrix i j x t

-- | Create new SpMatrix using data from a Foldable (e.g. a list) in (row, col, value) form
fromListSM :: Foldable t => (Int, Int) -> t (IxRow, IxCol, a) -> SpMatrix a
fromListSM (m,n) iix = fromListSM' iix (zeroSM m n)


mkSpMR :: Foldable t =>
   (Int, Int) -> t (IxRow, IxCol, Double) -> SpMatrix Double
mkSpMR d ixv = fromListSM d ixv :: SpMatrix Double

mkSpMC :: Foldable t =>
   (Int, Int) -> t (IxRow, IxCol, Complex Double) -> SpMatrix (Complex Double)
mkSpMC d ixv = fromListSM d ixv :: SpMatrix (Complex Double)




-- | Create new SpMatrix assuming contiguous, 0-based indexing of elements
fromListDenseSM :: Int -> [a] -> SpMatrix a
fromListDenseSM m ll = fromListSM (m, n) $ indexed2 m ll where
  n = length ll `div` m


-- ** toList

-- | Populate list with SpMatrix contents
toListSM :: SpMatrix t -> [(IxRow, IxCol, t)]
toListSM = ifoldlSM buildf [] where
  buildf i j x y = (i, j, x) : y 


-- | Populate list with SpMatrix contents and populate missing entries with 0
toDenseListSM :: Num t => SpMatrix t -> [(IxRow, IxCol, t)]
toDenseListSM m =
  [(i, j, m @@ (i, j)) | i <- [0 .. nrows m - 1], j <- [0 .. ncols m- 1]]










-- ** Lookup




lookupSM :: SpMatrix a -> IxRow -> IxCol -> Maybe a
lookupSM (SM _ im) i j = I.lookup i im >>= I.lookup j

-- | Looks up an element in the matrix with a default (if the element is not found, zero is returned)

lookupWD_SM, (@@!):: Num a => SpMatrix a -> (IxRow, IxCol) -> a
lookupWD_SM sm (i,j) =
  fromMaybe 0 (lookupSM sm i j)



-- | Zero-default lookup, infix form (no bound checking)
(@@!) = lookupWD_SM







-- FIXME : to throw an exception or just ignore the out-of-bound access ?







-- ** Sub-matrices

-- | Indexed filtering function
filterSM :: (IM.Key -> IM.Key -> a -> Bool) -> SpMatrix a -> SpMatrix a
filterSM f sm = SM (dim sm) $ ifilterIM2 f (dat sm)

-- | Diagonal, subdiagonal, superdiagonal partitions of a SpMatrix (useful for writing preconditioners)
extractDiag, extractSuperDiag, extractSubDiag :: SpMatrix a -> SpMatrix a
extractSubDiag = filterSM (\i j _ -> i > j)

extractSuperDiag = filterSM (\i j _ -> i < j)

extractDiag = filterSM (\i j _ -> i == j)




-- | Extract a submatrix given the specified index bounds, rebalancing keys with the two supplied functions
extractSubmatrixSM ::
  (IM.Key -> IM.Key) ->   -- row index function
  (IM.Key -> IM.Key) ->   -- column "  "
  SpMatrix a ->
  (IxRow, IxRow) -> (IxCol, IxCol) ->
  SpMatrix a
extractSubmatrixSM fi gi (SM (r, c) im) (i1, i2) (j1, j2)
  | q = SM (m', n') imm'
  | otherwise = error $ "extractSubmatrixSM : invalid index " ++ show (i1, i2) ++ ", " ++ show (j1, j2) where
  imm' = mapKeysIM2 fi gi $                          -- rebalance keys
          I.filterI (not . null) $                -- remove all-null rows
          ifilterIM2 ff im                           -- keep `submatrix`
  ff i j _ = i1 <= i &&
             i <= i2 &&
             j1 <= j &&
             j <= j2
  (m', n') = (i2-i1 + 1, j2-j1 + 1)
  q = inBounds0 r i1  &&
      inBounds0 r i2 &&
      inBounds0 c j1  &&
      inBounds0 c j2 &&      
      i2 >= i1

-- | Extract a submatrix given the specified index bounds
-- NB : subtracts (i1, j1) from the indices
extractSubmatrixRebalanceKeys ::
  SpMatrix a -> (IxRow, IxRow) -> (IxCol, IxCol) -> SpMatrix a
extractSubmatrixRebalanceKeys mm (i1,i2) (j1,j2) =
  extractSubmatrixSM (\i -> i - i1) (\j -> j - j1) mm (i1,i2) (j1,j2)

-- | Extract a submatrix given the specified index bounds
-- NB : submatrix indices are _preserved_
extractSubmatrix :: SpMatrix a -> (IxRow, IxRow) -> (IxCol, IxCol) -> SpMatrix a
extractSubmatrix = extractSubmatrixSM id id


takeRows :: IxRow -> SpMatrix a -> SpMatrix a
takeRows n mm = extractSubmatrix mm (0, n-1) (0, ncols mm - 1)

takeCols :: IxCol -> SpMatrix a -> SpMatrix a
takeCols n mm = extractSubmatrix mm (0, nrows mm - 1) (0, n - 1)



-- *** Extract i'th row
-- | Extract whole row
-- -- moved to Data.Sparse.Common




-- *** Extract j'th column
-- | Extract whole column
extractColSM :: SpMatrix a -> IxCol -> SpMatrix a
extractColSM sm j = extractSubmatrix sm (0, nrows sm - 1) (j, j)

-- | Extract column within a row range
extractSubColSM :: SpMatrix a -> IxCol -> (IxRow, IxRow) -> SpMatrix a
extractSubColSM sm j (i1, i2) = extractSubmatrix sm (i1, i2) (j, j)


-- | Extract column within a row range, rebalance keys
extractSubColSM_RK :: SpMatrix a -> IxCol -> (IxRow, IxRow) -> SpMatrix a
extractSubColSM_RK sm j (i1, i2) =
  extractSubmatrixRebalanceKeys sm (i1, i2) (j, j) 







  







-- ** Predicates
-- |Are the supplied indices within matrix bounds?
isValidIxSM :: SpMatrix a -> (Int, Int) -> Bool
isValidIxSM mm = inBounds02 (dim mm)

-- |Is the matrix square?
isSquareSM :: SpMatrix a -> Bool
isSquareSM m = nrows m == ncols m

-- |Is the matrix diagonal?
isDiagonalSM :: SpMatrix a -> Bool
isDiagonalSM m = I.size d == nrows m where
  d = I.filterWithKey ff (immSM m)
  ff irow row = I.size row == 1 &&
                I.size (I.filterWithKey (\j _ -> j == irow) row) == 1

-- | Is the matrix lower/upper triangular?
isLowerTriSM, isUpperTriSM :: Eq a => SpMatrix a -> Bool
isLowerTriSM m = m == lm where
  lm = ifilterSM (\i j _ -> i >= j) m

isUpperTriSM m = m == lm where
  lm = ifilterSM (\i j _ -> i <= j) m

-- |Is the matrix orthogonal? i.e. Q^t ## Q == I
isOrthogonalSM :: (Eq a, Epsilon a, MatrixRing (SpMatrix a)) => SpMatrix a -> Bool
isOrthogonalSM sm@(SM (_,n) _) = rsm == eye n where
  rsm = roundZeroOneSM $ transpose sm ## sm











-- ** Matrix data and metadata

-- | Data in internal representation (do not export)
-- immSM :: SpMatrix t -> IM.IntMap (IM.IntMap t)
immSM (SM _ imm) = imm

-- | (Number of rows, Number of columns)
dimSM :: SpMatrix t -> (Rows, Cols)
dimSM (SM d _) = d

-- | Number of rows times number of columns
nelSM :: SpMatrix t -> Int
nelSM (SM (nr,nc) _) = nr*nc

-- | Number of rows
nrows :: SpMatrix a -> Rows
nrows = fst . dim

-- | Number of columns
ncols :: SpMatrix a -> Cols
ncols = snd . dim

data SMInfo = SMInfo { smNz :: Int,
                       smSpy :: Double} deriving (Eq, Show)

infoSM :: SpMatrix a -> SMInfo
infoSM s = SMInfo (nzSM s) (spySM s)

nzSM :: SpMatrix a -> Int
nzSM s = sum $ fmap I.size (immSM s)

spySM :: Fractional b => SpMatrix a -> b
spySM s = fromIntegral (nzSM s) / fromIntegral (nelSM s)



-- ** Non-zero elements in a row

nzRow :: SpMatrix a -> IM.Key -> Int
nzRow s i | inBounds0 (nrows s) i = nzRowU s i
          | otherwise = error "nzRow : index out of bounds" where
              nzRowU :: SpMatrix a -> IM.Key -> Int
              nzRowU s i = maybe 0 I.size (I.lookup i $ immSM s)




-- ** Bandwidth bounds (min, max)

bwMinSM :: SpMatrix a -> Int
bwMinSM = fst . bwBoundsSM

bwMaxSM :: SpMatrix a -> Int
bwMaxSM = snd . bwBoundsSM

bwBoundsSM :: SpMatrix a -> (Int, Int)
bwBoundsSM s = -- b
                (snd $ I.findMin b,
                snd $ I.findMax b)
  where
  ss = immSM s
  fmi = fst . I.findMin
  fma = fst . I.findMax
  b = fmap (\x -> fma x - fmi x + 1:: Int) ss





















  






























-- ** Matrix stacking

-- | Vertical stacking of matrix blocks
vertStackSM, (-=-) :: SpMatrix a -> SpMatrix a -> SpMatrix a
vertStackSM mm1 mm2 = SM (m, n) $ I.union u1 u2 where
  nro1 = nrows mm1
  m = nro1 + nrows mm2
  n = max (ncols mm1) (ncols mm2)
  u1 = immSM mm1
  u2 = I.mapKeys (+ nro1) (immSM mm2)

(-=-) = vertStackSM

-- | Horizontal stacking of matrix blocks
horizStackSM, (-||-) :: SpMatrix a -> SpMatrix a -> SpMatrix a
horizStackSM mm1 mm2 = t (t mm1 -=- t mm2) where
  t = transposeSM

(-||-) = horizStackSM



-- | Assemble a block-diagonal square matrix from a list of square matrices, arranging these along the main diagonal
fromBlocksDiag :: [SpMatrix a] -> SpMatrix a
fromBlocksDiag mml = fromListSM (n, n) lstot where
  dims = map nrows mml
  n = sum dims
  shifts = init $ scanl (+) 0 dims
  lstot = concat $ zipWith shiftDims shifts $ map toListSM mml --lsts
  shiftDims s = map (\(i,j,x) -> (i + s, j + s, x))








-- ** Misc. SpMatrix operations

-- | Indexed filter over SpMatrix
{-# INLINE ifilterSM #-}
ifilterSM :: (IM.Key -> IM.Key -> a -> Bool) -> SpMatrix a -> SpMatrix a
ifilterSM f (SM d im) = SM d $ ifilterIM2 f im

 

-- | Left fold over SpMatrix
{-# INLINE foldlSM #-}
foldlSM :: (a -> b -> b) -> b -> SpMatrix a -> b
foldlSM f n (SM _ m)= foldlIM2 f n m

-- | Indexed left fold over SpMatrix
{-# INLINE ifoldlSM #-}
ifoldlSM :: (IM.Key -> IM.Key -> a -> b -> b) -> b -> SpMatrix a -> b
ifoldlSM f n (SM _ m) = ifoldlIM2' f n m







-- |Count sub-diagonal nonzeros
countSubdiagonalNZSM :: SpMatrix a -> Int
countSubdiagonalNZSM (SM _ im) = countSubdiagonalNZ im


-- extractDiagonalSM :: (Num a, Eq a) => SpMatrix a -> SpVector a
-- extractDiagonalSM (SM (m,n) im) = mkSpVectorD m $ extractDiagonalIM2 im


  




  

-- |Filter the index subset that lies below the diagonal (used in the QR decomposition, for example)
subdiagIndicesSM :: SpMatrix a -> [(IxRow, IxCol)]
subdiagIndicesSM (SM _ im) = subdiagIndices im





-- ** Sparsify : remove almost-0 elements (|x| < eps)
sparsifyIM2 ::
  Epsilon a => I.IntM (I.IntM a) -> I.IntM (I.IntM a)
sparsifyIM2 = ifilterIM2 (\_ _ x -> isNz x)

-- | Sparsify an SpMatrix
sparsifySM :: Epsilon a => SpMatrix a -> SpMatrix a
sparsifySM (SM d im) = SM d $ sparsifyIM2 im




-- ** Value rounding
-- | Round almost-0 and almost-1 to 0 and 1 respectively
roundZeroOneSM :: Epsilon a => SpMatrix a -> SpMatrix a
roundZeroOneSM (SM d im) = sparsifySM $ SM d $ mapIM2 roundZeroOne im  








-- | Modify (row, column) keys, leaving data intact. Be careful when using this!
-- modifyKeysSM' :: (IxRow -> IxRow) -> (IxCol -> IxCol) -> SpMatrix a -> SpMatrix a
modifyKeysSM' :: (IxRow -> a) -> (IxCol -> b) -> SpMatrix c -> [(a, b, c)]
modifyKeysSM' fi fj mm =  zip3 (fi <$> ii) (fj <$> jj) xx where
  (ii, jj, xx) = unzip3 $ toListSM mm


modifyKeysSM :: (IxRow -> IxRow) -> (IxCol -> IxCol) -> SpMatrix a -> SpMatrix a
modifyKeysSM fi fj mm = fromListSM (dim mm) $ zip3 (fi <$> ii) (fj <$> jj) xx where
  (ii, jj, xx) = unzip3 $ toListSM mm










-- * Primitive algebra operations


-- ** Matrix row swap
-- | Swap two rows of a SpMatrix (bounds not checked)
swapRows :: IxRow -> IxRow -> SpMatrix a -> SpMatrix a
swapRows i1 i2 (SM d im) = SM d $ I.insert i1 ro2 im' where
  ro1 = im I.! i1
  ro2 = im I.! i2
  im' = I.insert i2 ro1 im

-- | Swap two rows of a SpMatrix (bounds checked)  
swapRowsSafe :: IxRow -> IxRow -> SpMatrix a -> SpMatrix a
swapRowsSafe i1 i2 m
  | inBounds02 (nro, nro) (i1, i2) = swapRows i1 i2 m
  | otherwise =
     error $ "swapRowsSafe : index out of bounds " ++ show (i1, i2)
      where nro = nrows m  






-- ** Matrix transpose
-- | transposeSM : Matrix transpose
transposeSM :: SpMatrix a -> SpMatrix a
transposeSM (SM (m, n) im) = SM (n, m) (transposeIM2 im)

-- | Hermitian conjugate
hermitianConj :: Num a => SpMatrix (Complex a) -> SpMatrix (Complex a)
hermitianConj m = conjugate <$> transposeSM m






-- ** Multiply matrix by a scalar
matScale :: Num a => a -> SpMatrix a -> SpMatrix a
matScale a = fmap (* a)








-- ** Trace

-- | Matrix trace
trace :: Num b => SpMatrix b -> b
trace m = foldlSM (+) 0 $ extractDiag m




-- ** Frobenius norm

normFrobeniusSM :: (MatrixRing (SpMatrix a), Floating a) => SpMatrix a -> a
normFrobeniusSM m = sqrt $ trace (m ##^ m)

normFrobeniusSMC ::
  (MatrixRing (SpMatrix (Complex a)), RealFloat a) => SpMatrix (Complex a) -> a
normFrobeniusSMC m = sqrt $ magnitude $ trace (m ##^ m)









-- ** Matrix-matrix product

instance MatrixRing (SpMatrix Double) where
  type MatrixNorm (SpMatrix Double) = Double
  (##) = matMat_ AB
  (##^) = matMat_ ABt
  transpose = transposeSM
  normFrobenius = normFrobeniusSM

instance MatrixRing (SpMatrix (Complex Double)) where
  type MatrixNorm (SpMatrix (Complex Double)) = Double
  (##) = matMat_ AB
  (##^) = matMat_ ABt
  transpose = hermitianConj
  normFrobenius = normFrobeniusSMC


-- | Internal implementation
data MatProd_ = AB | ABt deriving (Eq, Show)

{-# INLINE matMat_ #-}
matMat_ pt mm1 mm2 =
  case pt of AB -> matMatCheck (matMatUnsafeWith transposeIM2) mm1 mm2
             ABt -> matMatCheck (matMatUnsafeWith id) mm1 (trDim mm2)
   where
     trDim (SM (a, b) x) = SM (b, a) x
     matMatCheck mmf m1 m2
       | c1 == r2 = mmf m1 m2
       | otherwise = error $ "matMat : incompatible matrix sizes" ++ show (d1, d2)
         where
           d1@(_, c1) = dim m1
           d2@(r2, _) = dim m2


-- | Matrix product without dimension checks
{-# INLINE matMatUnsafeWith #-}
-- matMatUnsafeWith :: Num a =>
--    (IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)) ->
--    SpMatrix a ->
--    SpMatrix a ->
--    SpMatrix a
matMatUnsafeWith ff2 m1 m2 = SM (nrows m1, ncols m2) (overRows2 <$> immSM m1) where
          overRows2 vm1 = (`dott` vm1) <$> ff2 (immSM m2)
          dott x y = sum $ liftI2 (*) x y    -- NB !! no complex conjugation






-- ** Matrix-matrix product, sparsified
-- | After multiplying the two matrices, all elements `x` for which `| x | <= eps` are removed.
matMatSparsified, (#~#) :: (MatrixRing (SpMatrix a), Epsilon a) =>
    SpMatrix a -> SpMatrix a -> SpMatrix a
matMatSparsified m1 m2 = sparsifySM $ m1 ## m2

(#~#) = matMatSparsified




-- *** Sparsified matrix products of two matrices

-- | Sparsifying A^T B
(#~#^) :: (MatrixRing (SpMatrix a), Epsilon a) =>
     SpMatrix a -> SpMatrix a -> SpMatrix a
a #~^# b = transpose a #~# b


-- | Sparsifying A B^T
(#~^#) :: (MatrixRing (SpMatrix a), Epsilon a) =>
     SpMatrix a -> SpMatrix a -> SpMatrix a
a #~#^ b = a #~# transpose b













-- ** Partial inner product
-- | Contract row `i` of A with column `j` of B up to an index `n`, i.e. summing over repeated indices: 
-- Aij Bjk , for j in [0 .. n] 
contractSub :: Elt a => SpMatrix a -> SpMatrix a -> IxRow -> IxCol -> Int -> a
contractSub a b i j n
  | ncols a == nrows b &&
    isValidIxSM a (i,j) &&
    n <= ncols a = sum $ map (\i' -> (a@@!(i,i'))*b@@!(i',j)) [0 .. n]
  | otherwise = error "contractSub : n must be <= i"







-- -- * Misc.utilities

-- encode :: (Int, Int) -> (Rows, Cols) -> Int
-- encode (nr,_) (i,j) = i + (j * nr)

-- decode :: (Int, Int) -> Int -> (Rows, Cols)
-- decode (nr, _) ci = (r, c) where (c,r ) = quotRem ci nr
