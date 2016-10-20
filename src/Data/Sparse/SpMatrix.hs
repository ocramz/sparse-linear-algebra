{-# language FlexibleInstances, MultiParamTypeClasses, TypeFamilies #-}
module Data.Sparse.SpMatrix where

import Data.Sparse.Utils

import Data.Sparse.Utils

import Numeric.Eps
import Numeric.LinearAlgebra.Class
import Numeric.LinearAlgebra.Data

import Numeric.LinearAlgebra.Sparse.IntMap

import qualified Data.IntMap as IM

import Data.Maybe



-- * Sparse Matrix

data SpMatrix a = SM {smDim :: (Rows, Cols),
                      smData :: IM.IntMap (IM.IntMap a)} deriving Eq


sizeStr :: SpMatrix a -> String
sizeStr sm =
  unwords ["(",show (nrows sm),"rows,",show (ncols sm),"columns ) ,",show nz,"NZ ( sparsity",show sy,")"] where
  (SMInfo nz sy) = infoSM sm

instance Show a => Show (SpMatrix a) where
  show sm@(SM _ x) = "SM " ++ sizeStr sm ++ " "++ show (IM.toList x)

instance Functor SpMatrix where
  fmap f (SM d md) = SM d ((fmap . fmap) f md)

instance Set SpMatrix where
  liftU2 f2 (SM n1 x1) (SM n2 x2) = SM (maxTup n1 n2) ((liftU2.liftU2) f2 x1 x2)
  liftI2 f2 (SM n1 x1) (SM n2 x2) = SM (minTup n1 n2) ((liftI2.liftI2) f2 x1 x2)
  
instance Additive SpMatrix where
  zero = SM (0,0) IM.empty
  (^+^) = liftU2 (+)


instance FiniteDim SpMatrix where
  type FDSize SpMatrix = (Rows, Cols)
  dim = smDim

instance HasData SpMatrix a where
  type HDData SpMatrix a = IM.IntMap (IM.IntMap a)
  dat = smData

instance Sparse SpMatrix a where
  spy = spySM






-- ** Creation

-- | `zeroSM m n` : Empty SpMatrix of size (m, n)
zeroSM :: Rows -> Cols -> SpMatrix a
zeroSM m n = SM (m,n) IM.empty


-- *** Diagonal matrix
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
-- [0,1,0,0,0]
-- [0,0,0,1,0]
-- [0,0,1,0,0]
-- [1,0,0,0,0]
-- [0,0,0,0,1]
permutationSM :: Num a => Int -> [IxRow] -> SpMatrix a
permutationSM n iis = permutPairsSM n (zip [0 .. n-1] iis)

-- | Permutation matrix from a (possibly incomplete) list of row pair swaps
-- e.g. `permutPairs 5 [(2,4)]` swaps rows (2, 4) :
--
-- [1,0,0,0,0]
-- [0,1,0,0,0]
-- [0,0,0,0,1]
-- [0,0,0,1,0]
-- [0,0,1,0,0]
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
fromListSM' iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertSpMatrix i j x t

-- | Create new SpMatrix using data from list (row, col, value)
fromListSM :: Foldable t => (Int, Int) -> t (IxRow, IxCol, a) -> SpMatrix a
fromListSM (m,n) iix = fromListSM' iix (zeroSM m n)


-- | Create new SpMatrix assuming contiguous, 0-based indexing of elements
fromListDenseSM :: Int -> [a] -> SpMatrix a
fromListDenseSM m ll = fromListSM (m, n) $ denseIxArray2 m ll where
  n = length ll `div` m


-- ** toList

-- |Populate list with SpMatrix contents and populate missing entries with 0
toDenseListSM :: Num t => SpMatrix t -> [(IxRow, IxCol, t)]
toDenseListSM m =
  [(i, j, m @@ (i, j)) | i <- [0 .. nrows m - 1], j <- [0 .. ncols m- 1]]












-- ** Lookup

lookupSM :: SpMatrix a -> IxRow -> IxCol -> Maybe a
lookupSM (SM _ im) i j = IM.lookup i im >>= IM.lookup j

-- | Looks up an element in the matrix with a default (if the element is not found, zero is returned)

lookupWD_SM, (@@) :: Num a => SpMatrix a -> (IxRow, IxCol) -> a
lookupWD_SM sm (i,j) =
  fromMaybe 0 (lookupSM sm i j)

lookupWD_IM :: Num a => IM.IntMap (IM.IntMap a) -> (IxRow, IxCol) -> a
lookupWD_IM im (i,j) = fromMaybe 0 (IM.lookup i im >>= IM.lookup j)

-- | Zero-default lookup, infix form
(@@) = lookupWD_SM




-- FIXME : to throw an exception or just ignore the out-of-bound access ?







-- ** Sub-matrices

-- | Extract a submatrix given the specified index bounds
extractSubmatrixSM :: SpMatrix a -> (IxRow, IxCol) -> (IxRow, IxCol) -> SpMatrix a
extractSubmatrixSM (SM (r, c) im) (i1, i2) (j1, j2)
  | q = SM (m', n') imm'
  | otherwise = error $ "extractSubmatrixSM : invalid indexing " ++ show (i1, i2) ++ ", " ++ show (j1, j2) where
  imm' = mapKeysIM2 (\i -> i - i1) (\j -> j - j1) $  -- rebalance keys
          IM.filter (not . IM.null) $                -- remove all-null rows
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




-- *** Extract j'th column
extractColSM :: SpMatrix a -> IxCol -> SpMatrix a
extractColSM sm j = extractSubmatrixSM sm (0, nrows sm - 1) (j, j)




-- *** Extract i'th row
extractRowSM :: SpMatrix a -> IxRow -> SpMatrix a
extractRowSM sm i = extractSubmatrixSM sm (i, i) (0, ncols sm - 1)






  







-- ** Predicates
-- |Are the supplied indices within matrix bounds?
validIxSM :: SpMatrix a -> (Int, Int) -> Bool
validIxSM mm = inBounds02 (dim mm)

-- |Is the matrix square?
isSquareSM :: SpMatrix a -> Bool
isSquareSM m = nrows m == ncols m

-- |Is the matrix diagonal?
isDiagonalSM :: SpMatrix a -> Bool
isDiagonalSM m = IM.size d == nrows m where
  d = IM.filterWithKey ff (immSM m)
  ff irow row = IM.size row == 1 &&
                IM.size (IM.filterWithKey (\j _ -> j == irow) row) == 1

-- |is the matrix orthogonal? i.e. Q^t ## Q == I
isOrthogonalSM :: SpMatrix Double -> Bool
isOrthogonalSM sm@(SM (_,n) _) = rsm == eye n where
  rsm = roundZeroOneSM $ transposeSM sm ## sm












-- ** Matrix data and metadata

-- | Data in internal representation (do not export)
immSM :: SpMatrix t -> IM.IntMap (IM.IntMap t)
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
nzSM s = sum $ fmap IM.size (immSM s)

spySM :: Fractional b => SpMatrix a -> b
spySM s = fromIntegral (nzSM s) / fromIntegral (nelSM s)



-- ** Non-zero elements in a row

nzRow :: SpMatrix a -> IM.Key -> Int
nzRow s i | inBounds0 (nrows s) i = nzRowU s i
          | otherwise = error "nzRow : index out of bounds" where
              nzRowU :: SpMatrix a -> IM.Key -> Int
              nzRowU s i = maybe 0 IM.size (IM.lookup i $ immSM s)




-- ** Bandwidth bounds (min, max)

bwMinSM :: SpMatrix a -> Int
bwMinSM = fst . bwBoundsSM

bwMaxSM :: SpMatrix a -> Int
bwMaxSM = snd . bwBoundsSM

bwBoundsSM :: SpMatrix a -> (Int, Int)
bwBoundsSM s = -- b
                (snd $ IM.findMin b,
                snd $ IM.findMax b)
  where
  ss = immSM s
  fmi = fst . IM.findMin
  fma = fst . IM.findMax
  b = fmap (\x -> fma x - fmi x + 1:: Int) ss





















  













-- encode :: (Int, Int) -> (Rows, Cols) -> Int
-- encode (nr,_) (i,j) = i + (j * nr)

-- decode :: (Int, Int) -> Int -> (Rows, Cols)
-- decode (nr, _) ci = (r, c) where (c,r ) = quotRem ci nr
















-- ** Matrix stacking

-- | Vertical stacking
vertStackSM, (-=-) :: SpMatrix a -> SpMatrix a -> SpMatrix a
vertStackSM mm1 mm2 = SM (m, n) $ IM.union u1 u2 where
  nro1 = nrows mm1
  m = nro1 + nrows mm2
  n = max (ncols mm1) (ncols mm2)
  u1 = immSM mm1
  u2 = IM.mapKeys (+ nro1) (immSM mm2)

(-=-) = vertStackSM

-- | Horizontal stacking
horizStackSM, (-||-) :: SpMatrix a -> SpMatrix a -> SpMatrix a
horizStackSM mm1 mm2 = t (t mm1 -=- t mm2) where
  t = transposeSM

(-||-) = horizStackSM

















-- ** Misc. SpMatrix operations

-- | Left fold over SpMatrix
foldlSM :: (a -> b -> b) -> b -> SpMatrix a -> b
foldlSM f n (SM _ m)= foldlIM2 f n m

-- | Indexed left fold over SpMatrix
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
sparsifyIM2 :: IM.IntMap (IM.IntMap Double) -> IM.IntMap (IM.IntMap Double)
sparsifyIM2 = ifilterIM2 (\_ _ x -> abs x >= eps)

-- | Sparsify an SpMatrix
sparsifySM :: SpMatrix Double -> SpMatrix Double
sparsifySM (SM d im) = SM d $ sparsifyIM2 im




-- ** Value rounding
-- | Round almost-0 and almost-1 to 0 and 1 respectively
roundZeroOneSM :: SpMatrix Double -> SpMatrix Double
roundZeroOneSM (SM d im) = sparsifySM $ SM d $ mapIM2 roundZeroOne im  






-- * Primitive algebra operations


-- ** Matrix row swap
-- | swap two rows of a SpMatrix (bounds not checked)
swapRows :: IxRow -> IxRow -> SpMatrix a -> SpMatrix a
swapRows i1 i2 (SM d im) = SM d $ IM.insert i1 ro2 im' where
  ro1 = im IM.! i1
  ro2 = im IM.! i2
  im' = IM.insert i2 ro1 im

-- | swap two rows of a SpMatrix (bounds checked)  
swapRowsSafe :: IxRow -> IxRow -> SpMatrix a -> SpMatrix a
swapRowsSafe i1 i2 m
  | inBounds02 (nro, nro) (i1, i2) = swapRows i1 i2 m
  | otherwise =
     error $ "swapRowsSafe : index out of bounds " ++ show (i1, i2)
      where nro = nrows m  






-- ** Matrix transpose
-- | transposeSM, (#^) : Matrix transpose
transposeSM, (#^) :: SpMatrix a -> SpMatrix a
transposeSM (SM (m, n) im) = SM (n, m) (transposeIM2 im)

(#^) = transposeSM






-- ** Multiply matrix by a scalar
matScale :: Num a => a -> SpMatrix a -> SpMatrix a
matScale a = fmap (*a)

-- ** Frobenius norm
normFrobenius :: SpMatrix Double -> Double
normFrobenius m = sqrt $ foldlSM (+) 0 m' where
  m' | nrows m > ncols m = transposeSM m ## m
     | otherwise = m ## transposeSM m 



















-- ** Matrix-matrix product

matMat, (##) :: Num a => SpMatrix a -> SpMatrix a -> SpMatrix a
matMat m1 m2
  | c1 == r2 = matMatU m1 m2
  | otherwise = error $ "matMat : incompatible matrix sizes" ++ show (d1, d2) where
      d1@(r1, c1) = dim m1
      d2@(r2, c2) = dim m2
      matMatU :: Num a => SpMatrix a -> SpMatrix a -> SpMatrix a
      matMatU m1 m2 =
        SM (nrows m1, ncols m2) im where
          im = fmap (\vm1 -> (`dot` vm1) <$> transposeIM2 (immSM m2)) (immSM m1)
    

(##) = matMat

-- matMat m1 m2 =
--   withDim2 m1 m2
--     (\(r1,c1) (r2,c2) _ _ -> c1 == r2)
--     matMatU
--     "matMat : incompatible matrix sizes"
--     (\m1 m2 -> unwords [show (dim m1), show (dim m2)])





-- ** Matrix-matrix product, sparsified
-- | Removes all elements `x` for which `| x | <= eps`)
matMatSparsified, (#~#)  :: SpMatrix Double -> SpMatrix Double -> SpMatrix Double
matMatSparsified m1 m2 = sparsifySM $ matMat m1 m2

(#~#) = matMatSparsified




-- *** Sparsified matrix products of two matrices

-- | A^T B
(#^#) :: SpMatrix Double -> SpMatrix Double -> SpMatrix Double
a #^# b = transposeSM a #~# b


-- | A B^T
(##^) :: SpMatrix Double -> SpMatrix Double -> SpMatrix Double
a ##^ b = a #~# transposeSM b














