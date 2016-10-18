{-# LANGUAGE FlexibleContexts, TypeFamilies, MultiParamTypeClasses, FlexibleInstances #-}
-- {-# OPTIONS_GHC -O2 -rtsopts -with-rtsopts=-K32m -prof#-}

module Numeric.LinearAlgebra.Sparse where

import Numeric.LinearAlgebra.Sparse.IntMap 


import Control.Monad.Primitive

import Control.Monad (mapM_, forM_, replicateM)

import Control.Monad.State.Strict
-- import Control.Monad.Writer
-- import Control.Monad.Trans

-- import Control.Monad.Trans.State (runStateT)
-- import Control.Monad.Trans.Writer (runWriterT)

import qualified Data.IntMap.Strict as IM
-- import Data.Utils.StrictFold (foldlStrict) -- hidden in `containers`

import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWC

import Data.Monoid
import qualified Data.Foldable as F
import qualified Data.Traversable as T


-- import qualified Data.List as L
import Data.Maybe



{-|  CLASSES and common operations -}

-- * Additive ring 
class Functor f => Additive f where
  -- | Ring zero element
  zero :: Num a => f a
  
  -- | Ring +
  (^+^) :: Num a => f a -> f a -> f a




-- | negate the values in a functor
negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate

-- | subtract two Additive objects
(^-^) :: (Additive f, Num a) => f a -> f a -> f a
x ^-^ y = x ^+^ negated y






-- * Vector space
class Additive f => VectorSpace f where
  -- | multiplication by a scalar
  (.*) :: Num a => a -> f a -> f a
  

-- |linear interpolation
lerp :: (VectorSpace f, Num a) => a -> f a -> f a -> f a
lerp a u v = a .* u ^+^ ((1-a) .* v)


-- * Hilbert space (inner product)
class VectorSpace f => Hilbert f where
  -- | inner product
  dot :: Num a => f a -> f a -> a


-- ** Hilbert-space distance function
-- |`hilbertDistSq x y = || x - y ||^2`
hilbertDistSq :: (Hilbert f, Num a) => f a -> f a -> a
hilbertDistSq x y = dot t t where
  t = x ^-^ y

  


-- * Normed vector space
class Hilbert f => Normed f where
  norm :: (Floating a, Eq a) => a -> f a -> a




-- ** Norms and related results

-- | Squared 2-norm
normSq :: (Hilbert f, Num a) => f a -> a
normSq v = v `dot` v


-- |L1 norm
norm1 :: (Foldable t, Num a, Functor t) => t a -> a
norm1 v = sum (fmap abs v)

-- |Euclidean norm
norm2 :: (Hilbert f, Floating a) => f a -> a
norm2 v = sqrt (normSq v)

-- |Lp norm (p > 0)
normP :: (Foldable t, Functor t, Floating a) => a -> t a -> a
normP p v = sum u**(1/p) where
  u = fmap (**p) v

-- |Infinity-norm
normInfty :: (Foldable t, Ord a) => t a -> a
normInfty = maximum



-- |Normalize w.r.t. p-norm (p finite)
normalize :: (Normed f, Floating a, Eq a) => a -> f a -> f a
normalize p v = (1 / norm p v) .* v







-- |Lp inner product (p > 0)
dotLp :: (Set t, Foldable t, Floating a) => a -> t a -> t a ->  a
dotLp p v1 v2 = sum u**(1/p) where
  f a b = (a*b)**p
  u = liftI2 f v1 v2


-- |Reciprocal
reciprocal :: (Functor f, Fractional b) => f b -> f b
reciprocal = fmap recip


-- |Scale
scale :: (Num b, Functor f) => b -> f b -> f b
scale n = fmap (* n)







-- * FiniteDim : finite-dimensional objects

class Additive f => FiniteDim f where
  type FDSize f :: *
  dim :: f a -> FDSize f


-- | unary dimension-checking bracket
withDim :: (FiniteDim f, Show e) =>
     f a
     -> (FDSize f -> f a -> Bool)
     -> (f a -> c)
     -> String
     -> (f a -> e)
     -> c
withDim x p f e ef | p (dim x) x = f x
                   | otherwise = error e' where e' = e ++ show (ef x)

-- | binary dimension-checking bracket
withDim2 :: (FiniteDim f, FiniteDim g, Show e) =>
     f a
     -> g b
     -> (FDSize f -> FDSize g -> f a -> g b -> Bool)
     -> (f a -> g b -> c)
     -> String
     -> (f a -> g b -> e)
     -> c
withDim2 x y p f e ef | p (dim x) (dim y) x y = f x y
                      | otherwise = error e' where e' = e ++ show (ef x y)






-- * HasData : accessing inner data (do not export)

class Additive f => HasData f a where
  type HDData f a :: * 
  dat :: f a -> HDData f a


-- * Sparse : sparse datastructures

class (FiniteDim f, HasData f a) => Sparse f a where
  spy :: Fractional b => f a -> b




-- * Set : things that behave as sets

class Functor f => Set f where
  -- |union binary lift : apply function on _union_ of two Sets
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- |intersection binary lift : apply function on _intersection_ of two Sets
  liftI2 :: (a -> b -> c) -> f a -> f b -> f c  



-- class (Set f, Sparse f a) => SparseSet f a

-- instance SparseSet SpVector a where




instance Set IM.IntMap where
  liftU2 = IM.unionWith
  {-# INLINE liftU2 #-}
  liftI2 = IM.intersectionWith
  {-# INLINE liftI2 #-}

instance Additive IM.IntMap where
  zero = IM.empty
  {-# INLINE zero #-}
  (^+^) = liftU2 (+)
  {-# INLINE (^+^) #-}


instance VectorSpace IM.IntMap where
  n .* im = IM.map (* n) im
  
instance Hilbert IM.IntMap where
   a `dot` b = sum $ liftI2 (*) a b
              

instance Normed IM.IntMap where
  norm p v | p==1 = norm1 v
           | p==2 = norm2 v
           | otherwise = normP p v





-- * Sparse Vector

data SpVector a = SV { svDim :: Int ,
                       svData :: IM.IntMap a} deriving Eq

-- | SpVector sparsity
spySV :: Fractional b => SpVector a -> b
spySV s = fromIntegral (IM.size (dat s)) / fromIntegral (dim s)






instance Functor SpVector where
  fmap f (SV n x) = SV n (fmap f x)

instance Set SpVector where  
  liftU2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftU2 f2 x1 x2)
  liftI2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftI2 f2 x1 x2)
  
instance Foldable SpVector where
    foldr f d v = F.foldr f d (svData v)

instance Additive SpVector where
  zero = SV 0 IM.empty
  (^+^) = liftU2 (+)


                      
instance VectorSpace SpVector where
  n .* v = scale n v


instance FiniteDim SpVector where
  type FDSize SpVector = Int
  dim = svDim  

instance HasData SpVector a where
  type HDData SpVector a = IM.IntMap a
  dat = svData

instance Sparse SpVector a where
  spy = spySV


instance Hilbert SpVector where
  a `dot` b | dim a == dim b = dot (dat a) (dat b)
            | otherwise =
                     error $ "dot : sizes must coincide, instead we got " ++
                           show (dim a, dim b)


instance Normed SpVector where
  norm p (SV _ v) = norm p v








-- ** Creation

-- | empty sparse vector (length n, no entries)
zeroSV :: Int -> SpVector a
zeroSV n = SV n IM.empty


-- | singleton sparse vector (length 1)
singletonSV :: a -> SpVector a
singletonSV x = SV 1 (IM.singleton 0 x)



-- | create a sparse vector from an association list while discarding all zero entries
mkSpVector :: (Num a, Eq a) => Int -> IM.IntMap a -> SpVector a
mkSpVector d im = SV d $ IM.filterWithKey (\k v -> v /= 0 && inBounds0 d k) im

-- | ", from logically dense array (consecutive indices)
mkSpVectorD :: (Num a, Eq a) => Int -> [a] -> SpVector a
mkSpVectorD d ll = mkSpVector d (IM.fromList $ denseIxArray (take d ll))

-- ", don't filter zero elements
mkSpVector1 :: Int -> IM.IntMap a -> SpVector a
mkSpVector1 d ll = SV d $ IM.filterWithKey (\ k _ -> inBounds0 d k) ll

-- | Create new sparse vector, assumin 0-based, contiguous indexing
fromListDenseSV :: Int -> [a] -> SpVector a
fromListDenseSV d ll = SV d (IM.fromList $ denseIxArray (take d ll))


-- | one-hot encoding : `oneHotSV n k` produces a SpVector of length n having 1 at the k-th position
oneHotSVU :: Num a => Int -> IxRow -> SpVector a
oneHotSVU n k = SV n (IM.singleton k 1)

oneHotSV :: Num a => Int -> IxRow -> SpVector a
oneHotSV n k |inBounds0 n k = oneHotSVU n k
             |otherwise = error "`oneHotSV n k` must satisfy 0 <= k <= n"


-- | DENSE vector of `1`s
onesSV :: Num a => Int -> SpVector a
onesSV d = SV d $ IM.fromList $ denseIxArray $ replicate d 1

-- | DENSE vector of `0`s
zerosSV :: Num a => Int -> SpVector a
zerosSV d = SV d $ IM.fromList $ denseIxArray $ replicate d 0



-- ** Element insertion

-- |insert element `x` at index `i` in a preexisting SpVector
insertSpVector :: Int -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim)
  | inBounds0 d i = SV d (IM.insert i x xim)
  | otherwise = error "insertSpVector : index out of bounds"



-- ** fromList
fromListSV :: Int -> [(Int, a)] -> SpVector a
fromListSV d iix = SV d (IM.fromList (filter (inBounds0 d . fst) iix ))

-- ** toList
toListSV :: SpVector a -> [(IM.Key, a)]
toListSV sv = IM.toList (dat sv)

-- |To dense list (default = 0)
toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d im) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]









  
instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (IM.toList x)


-- ** Lookup

-- |lookup an index in a SpVector (returns 0 if lookup fails)
lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV i (SV _ im) = IM.findWithDefault 0 i im 






-- ** Sub-vectors
-- | Tail elements
tailSV :: SpVector a -> SpVector a
tailSV (SV n sv) = SV (n-1) ta where
  ta = IM.mapKeys (\i -> i - 1) $ IM.delete 0 sv
  
-- | Head element
headSV :: Num a => SpVector a -> a
headSV sv = fromMaybe 0 (IM.lookup 0 (dat sv))



-- | concatenate two sparse vectors
concatSV :: SpVector a -> SpVector a -> SpVector a
concatSV (SV n1 s1) (SV n2 s2) = SV (n1+n2) (IM.union s1 s2') where
  s2' = IM.mapKeys (+ n1) s2










-- | promote a SV to SM
svToSM :: SpVector a -> SpMatrix a
svToSM (SV n d) = SM (n, 1) $ IM.singleton 0 d



    

-- ** Outer vector product

outerProdSV, (><) :: Num a => SpVector a -> SpVector a -> SpMatrix a
outerProdSV v1 v2 = fromListSM (m, n) ixy where
  m = dim v1
  n = dim v2
  ixy = [(i,j, x * y) | (i,x) <- toListSV v1 , (j, y) <- toListSV v2]

(><) = outerProdSV



-- * Sparsify : remove almost-0 elements (|x| < eps)
-- | Sparsify an SpVector
sparsifySV :: SpVector Double -> SpVector Double
sparsifySV (SV d im) = SV d $ IM.filter (\x -> abs x >= eps) im








-- * Sparse Matrix



data SpMatrix a = SM {smDim :: (Rows, Cols),
                      smData :: IM.IntMap (IM.IntMap a)} deriving Eq


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
  

-- | Componentwise tuple operations
-- TODO : use semilattice properties instead
maxTup, minTup :: Ord t => (t, t) -> (t, t) -> (t, t)
maxTup (x1,y1) (x2,y2) = (max x1 x2, max y1 y2)
minTup (x1,y1) (x2,y2) = (min x1 x2, min y1 y2)






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

-- | permutation matrix from a (possibly incomplete) list of row swaps
permutationSM :: Num a => Int -> [IxRow] -> SpMatrix a
permutationSM n iis = permut (zip [0 .. n-1] iis) (eye n) where
  permut ((i1,i2):iis) m = permut iis (swapRows i1 i2 m)
  permut [] m = m




              

-- permutationIx


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

-- |Demote (n x 1) or (1 x n) SpMatrix to SpVector
toSV :: SpMatrix a -> SpVector a
toSV (SM (m,n) im) = SV d $ snd . head $ IM.toList im where
  d | m==1 && n==1 = 1
    | m==1 && n>1 = n 
    | n==1 && m>1 = m
    | otherwise = error $ "toSV : incompatible dimensions " ++ show (m,n)


-- *** Extract j'th column
extractColSM :: SpMatrix a -> IxCol -> SpMatrix a
extractColSM sm j = extractSubmatrixSM sm (0, nrows sm - 1) (j, j)

-- |", and place into SpVector
extractCol :: SpMatrix a -> IxCol -> SpVector a
extractCol m j = toSV $ extractColSM m j


-- *** Extract i'th row
extractRowSM :: SpMatrix a -> IxRow -> SpMatrix a
extractRowSM sm i = extractSubmatrixSM sm (i, i) (0, ncols sm - 1)

-- |", and place into SpVector
extractRow :: SpMatrix a -> IxRow -> SpVector a
extractRow m i = toSV $ extractRowSM m i




  







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

-- | Extract the diagonal as a SpVector (with default 0)
extractDiagonalDSM :: Num a => SpMatrix a -> SpVector a
extractDiagonalDSM mm = fromListDenseSV n $ foldr ins [] ll  where
  ll = [0 .. n - 1]
  n = nrows mm
  ins i acc = mm@@(i,i) : acc
  




  

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







-- * Numerical tolerance for "near-0" tests
-- | eps = 1e-8 
eps :: Double
eps = 1e-8

  



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







-- ** Matrix action on a vector

{- 
FIXME : matVec is more general than SpVector's :

\m v -> fmap (`dot` v) m
  :: (Normed f1, Num b, Functor f) => f (f1 b) -> f1 b -> f b
-}



-- |Matrix-on-vector
matVec, (#>) :: Num a => SpMatrix a -> SpVector a -> SpVector a
matVec (SM (nr, nc) mdata) (SV n sv)
  | nc == n = SV nr $ fmap (`dot` sv) mdata
  | otherwise = error $ "matVec : mismatching dimensions " ++ show (nc, n)

(#>) = matVec

-- |Vector-on-matrix (FIXME : transposes matrix: more costly than `matVec`, I think)
vecMat, (<#) :: Num a => SpVector a -> SpMatrix a -> SpVector a  
vecMat (SV n sv) (SM (nr, nc) mdata)
  | n == nr = SV nc $ fmap (`dot` sv) (transposeIM2 mdata)
  | otherwise = error $ "vecMat : mismatching dimensions " ++ show (n, nr)

(<#) = vecMat











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
















-- * Matrix condition number

-- |uses the R matrix from the QR factorization
conditionNumberSM :: SpMatrix Double -> Double
conditionNumberSM m | isInfinite kappa = error "Infinite condition number : rank-deficient system"
                    | otherwise = kappa where
  kappa = lmax / lmin
  (_, r) = qr m
  u = extractDiagonalDSM r  -- FIXME : need to extract with default element 0 
  lmax = abs (maximum u)
  lmin = abs (minimum u)







-- * Householder transformation

hhMat :: Num a => a -> SpVector a -> SpMatrix a
hhMat beta x = eye n ^-^ scale beta (x >< x) where
  n = dim x


{-| a vector `x` uniquely defines an orthogonal plane; the Householder operator reflects any point `v` with respect to this plane:
 v' = (I - 2 x >< x) v
-}
hhRefl :: SpVector Double -> SpMatrix Double
hhRefl = hhMat 2.0











-- * Givens rotation matrix


hypot :: Floating a => a -> a -> a
hypot x y = abs x * (sqrt (1 + y/x)**2)

sign :: (Ord a, Num a) => a -> a
sign x
  | x > 0 = 1
  | x == 0 = 0
  | otherwise = -1 

-- | Givens coefficients (using stable algorithm shown in  Anderson, Edward (4 December 2000). "Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem". LAPACK Working Note)
givensCoef :: (Ord a, Floating a) => a -> a -> (a, a, a)
givensCoef a b  -- returns (c, s, r) where r = norm (a, b)
  | b==0 = (sign a, 0, abs a)
  | a==0 = (0, sign b, abs b)
  | abs a > abs b = let t = b/a
                        u = sign a * abs ( sqrt (1 + t**2))
                      in (1/u, - t/u, a*u)
  | otherwise = let t = a/b
                    u = sign b * abs ( sqrt (1 + t**2))
                in (t/u, - 1/u, b*u)


{- |
Givens method, row version: choose other row index i' s.t. i' is :
* below the diagonal
* corresponding element is nonzero

QR.C1 ) To zero out entry A(i, j) we must find row k such that A(k, j) is
non-zero but A has zeros in row k for all columns less than j.
-}

givens :: SpMatrix Double -> IxRow -> IxCol -> SpMatrix Double
givens mm i j 
  | validIxSM mm (i,j) && isSquareSM mm =
       sparsifySM $ fromListSM' [(i,i,c),(j,j,c),(j,i,-s),(i,j,s)] (eye (nrows mm))
  | otherwise = error "givens : indices out of bounds"      
  where
    (c, s, _) = givensCoef a b
    i' = head $ fromMaybe (error $ "givens: no compatible rows for entry " ++ show (i,j)) (candidateRows (immSM mm) i j)
    a = mm @@ (i', j)
    b = mm @@ (i, j)   -- element to zero out

-- |Is the `k`th the first nonzero column in the row?
firstNonZeroColumn :: IM.IntMap a -> IxRow -> Bool
firstNonZeroColumn mm k = isJust (IM.lookup k mm) &&
                          isNothing (IM.lookupLT k mm)

-- |Returns a set of rows {k} that satisfy QR.C1
candidateRows :: IM.IntMap (IM.IntMap a) -> IxRow -> IxCol -> Maybe [IM.Key]
candidateRows mm i j | IM.null u = Nothing
                     | otherwise = Just (IM.keys u) where
  u = IM.filterWithKey (\irow row -> irow /= i &&
                                     firstNonZeroColumn row j) mm





-- * QR decomposition


-- | Applies Givens rotation iteratively to zero out sub-diagonal elements
qr :: SpMatrix Double -> (SpMatrix Double, SpMatrix Double)
qr mm = (transposeSM qmatt, rmat)  where
  qmatt = F.foldl' (#~#) ee $ gmats mm -- Q^T = (G_n * G_n-1 ... * G_1)
  rmat = qmatt #~# mm                  -- R = Q^T A
  ee = eye (nrows mm)
      
-- | Givens matrices in order [G1, G2, .. , G_N ]
gmats :: SpMatrix Double -> [SpMatrix Double]
gmats mm = gm mm (subdiagIndicesSM mm) where
 gm m ((i,j):is) = let g = givens m i j
                   in g : gm (g #~# m) is
 gm _ [] = []





-- -- | QR algorithm, state transformer version
-- gmatST0 (m, (i,j):is) = (m', is) where    -- WRONG, possible access to []
--   g = givens m i j                        
--   m' = g #~# m
-- gmatST0 (m, []) = (eye (nrows m), [])

-- gmatST m = gmatST0 (m, subdiagIndicesSM m)






-- * Eigenvalue algorithms

-- ** All eigenvalues (QR algorithm)

-- | `eigsQR n mm` performs `n` iterations of the QR algorithm on matrix `mm` 
eigsQR :: Int -> SpMatrix Double -> SpVector Double
eigsQR nitermax m = extractDiagonalDSM $ execState (convergtest eigsStep) m where
  eigsStep m = r #~# q where (q, r) = qr m
  convergtest g = modifyInspectN nitermax f g where
    f [m1, m2] = let dm1 = extractDiagonalDSM m1
                     dm2 = extractDiagonalDSM m2
                 in norm2 (dm1 ^-^ dm2) <= eps






-- ** One eigenvalue and eigenvector (Rayleigh iteration)


-- | `eigsRayleigh n mm` performs `n` iterations of the Rayleigh algorithm on matrix `mm`. Cubic-order convergence, but it requires a mildly educated guess on the initial eigenpair
eigRayleigh :: Int                -- max # iterations
     -> SpMatrix Double           -- matrix
     -> (SpVector Double, Double) -- initial guess of (eigenvector, eigenvalue)
     -> (SpVector Double, Double) -- final estimate of (eigenvector, eigenvalue)
eigRayleigh nitermax m = execState (convergtest (rayleighStep m)) where
  convergtest g = modifyInspectN nitermax f g where
    f [(b1, _), (b2, _)] = norm2 (b2 ^-^ b1) <= eps 
  rayleighStep aa (b, mu) = (b', mu') where
      ii = eye (nrows aa)
      nom = (aa ^-^ (mu `matScale` ii)) <\> b
      b' = normalize 2 nom
      mu' = b' `dot` (aa #> b') / (b' `dot` b')




-- * Householder vector 

-- (Golub & Van Loan, Alg. 5.1.1, function `house`)
hhV :: SpVector Double -> (SpVector Double, Double)
hhV x = (v, beta) where
  n = dim x
  tx = tailSV x
  sigma = tx `dot` tx
  vtemp = singletonSV 1 `concatSV` tx
  (v, beta) | sigma <= eps = (vtemp, 0)
            | otherwise = let mu = sqrt (headSV x**2 + sigma)
                              xh = headSV x
                              vh | xh <= 1 = xh - mu
                                 | otherwise = - sigma / (xh + mu)
                              vnew = (1 / vh) .* insertSpVector 0 vh vtemp     
                          in (vnew, 2 * xh**2 / (sigma + vh**2))

                         



-- * Householder bidiagonalization

{- G & VL Alg. 5.4.2 -}






-- * SVD

{- Golub & Van Loan, sec 8.6.2 (p 452 segg.)

SVD of A, Golub-Kahan method

* reduce A to upper bidiagonal form B (Alg. 5.4.2, Householder bidiagonalization)
* compute SVD of B (implicit-shift QR step applied to B^T B, Alg. 8.3.2)

-}















-- * LU factorization















-- * Iterative linear solvers



-- -- | residual of candidate solution x0 of a linear system
-- residual :: Num a => SpMatrix a -> SpVector a -> SpVector a -> SpVector a
-- residual aa b x0 = b ^-^ (aa #> x0)

-- converged :: SpMatrix Double -> SpVector Double -> SpVector Double -> Bool
-- converged aa b x0 = normSq (residual aa b x0) <= eps



-- ** CGS

-- | one step of CGS
cgsStep :: SpMatrix Double -> SpVector Double -> CGS -> CGS
cgsStep aa rhat (CGS x r p u) = CGS xj1 rj1 pj1 uj1
  where
  aap = aa #> p
  alphaj = (r `dot` rhat) / (aap `dot` rhat)
  q = u ^-^ (alphaj .* aap)
  xj1 = x ^+^ (alphaj .* (u ^+^ q))  -- updated solution
  rj1 = r ^-^ (alphaj .* (aa #> (u ^+^ q)))-- updated residual
  betaj = (rj1 `dot` rhat) / (r `dot` rhat)
  uj1 = rj1 ^+^ (betaj .* q)
  pj1 = uj1 ^+^ (betaj .* (q ^+^ (betaj .* p)))

data CGS = CGS { _x :: SpVector Double,
                 _r :: SpVector Double,
                 _p :: SpVector Double,
                 _u :: SpVector Double } deriving Eq

-- | iterate solver until convergence or until max # of iterations is reached
cgs ::
  SpMatrix Double ->
  SpVector Double ->
  SpVector Double ->
  SpVector Double ->
  CGS
cgs aa b x0 rhat =
  execState (untilConverged _x (cgsStep aa rhat)) cgsInit where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  p0 = r0
  u0 = r0
  cgsInit = CGS x0 r0 p0 u0


instance Show CGS where
  show (CGS x r p u) = "x = " ++ show x ++ "\n" ++
                                "r = " ++ show r ++ "\n" ++
                                "p = " ++ show p ++ "\n" ++
                                "u = " ++ show u ++ "\n"





-- ** BiCGSTAB

-- _aa :: SpMatrix Double,    -- matrix
-- _b :: SpVector Double,     -- rhs
-- _r0 :: SpVector Double,    -- initial residual
-- _r0hat :: SpVector Double, -- candidate solution: r0hat `dot` r0 >= 0

-- | one step of BiCGSTAB
bicgstabStep :: SpMatrix Double -> SpVector Double -> BICGSTAB -> BICGSTAB
bicgstabStep aa r0hat (BICGSTAB x r p) = BICGSTAB xj1 rj1 pj1 where
  aap = aa #> p
  alphaj = (r `dot` r0hat) / (aap `dot` r0hat)
  sj = r ^-^ (alphaj .* aap)
  aasj = aa #> sj
  omegaj = (aasj `dot` sj) / (aasj `dot` aasj)
  xj1 = x ^+^ (alphaj .* p) ^+^ (omegaj .* sj)
  rj1 = sj ^-^ (omegaj .* aasj)
  betaj = (rj1 `dot` r0hat)/(r `dot` r0hat) * alphaj / omegaj
  pj1 = rj1 ^+^ (betaj .* (p ^-^ (omegaj .* aap)))

data BICGSTAB = BICGSTAB { _xBicgstab :: SpVector Double,
                           _rBicgstab :: SpVector Double,
                           _pBicgstab :: SpVector Double} deriving Eq

-- | iterate solver until convergence or until max # of iterations is reached
bicgstab
  :: SpMatrix Double
     -> SpVector Double
     -> SpVector Double
     -> SpVector Double
     -> BICGSTAB
bicgstab aa b x0 r0hat =
  execState (untilConverged _xBicgstab (bicgstabStep aa r0hat)) bicgsInit where
   r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
   p0 = r0
   bicgsInit = BICGSTAB x0 r0 p0

instance Show BICGSTAB where
  show (BICGSTAB x r p) = "x = " ++ show x ++ "\n" ++
                                "r = " ++ show r ++ "\n" ++
                                "p = " ++ show p ++ "\n"





-- * Moore-Penrose pseudoinverse
-- | Least-squares approximation of a rectangular system of equaitons. Uses <\\> for the linear solve
pinv :: SpMatrix Double -> SpVector Double -> SpVector Double
pinv aa b = aa #^# aa <\> atb where
  atb = transposeSM aa #> b







-- * Linear solver interface

data LinSolveMethod = CGS_ | BICGSTAB_ deriving (Eq, Show) 

-- | Linear solve with _random_ starting vector
linSolveM ::
  PrimMonad m =>
    LinSolveMethod -> SpMatrix Double -> SpVector Double -> m (SpVector Double)
linSolveM method aa b = do
  let (m, n) = dim aa
      nb     = dim b
  if n /= nb then error "linSolve : operand dimensions mismatch" else do
    x0 <- randVec nb
    case method of CGS_ -> return $ _xBicgstab (bicgstab aa b x0 x0)
                   BICGSTAB_ -> return $ _x (cgs aa b x0 x0)

-- | Linear solve with _deterministic_ starting vector (every component at 0.1) 
linSolve ::
  LinSolveMethod -> SpMatrix Double -> SpVector Double -> SpVector Double
linSolve method aa b
  | n /= nb = error "linSolve : operand dimensions mismatch"
  | otherwise = solve aa b where
      solve aa' b' | isDiagonalSM aa = (reciprocal aa') #> b'
                   | otherwise = solveWith aa' b' 
      solveWith aa' b' = case method of
                                CGS_ ->  _xBicgstab (bicgstab aa' b' x0 x0)
                                BICGSTAB_ -> _x (cgs aa' b' x0 x0)
      x0 = mkSpVectorD n $ replicate n 0.1 
      (m, n) = dim aa
      nb     = dim b

-- | <\\> : linSolve using the BiCGSTAB method as default
(<\>) :: SpMatrix Double -> SpVector Double -> SpVector Double      
(<\>) = linSolve BICGSTAB_ 
  







-- | TODO : if system is poorly conditioned, is it better to warn the user or just switch solvers (e.g. via the pseudoinverse) ?

-- linSolveQR aa b init f1 stepf
--   | isInfinite k = do
--        tell "linSolveQR : rank-deficient system"
--   | otherwise = do
--        solv aa b init
--     where
--      (q, r) = qr aa
--      k = conditionNumberSM r
--      solv aa b init = execState (untilConverged f1 stepf) init




















-- * Control primitives for bounded iteration with convergence check

-- | transform state until a condition is met
modifyUntil :: MonadState s m => (s -> Bool) -> (s -> s) -> m s
modifyUntil q f = do
  x <- get
  let y = f x
  put y
  if q y then return y
         else modifyUntil q f     

-- | Keep a moving window buffer (length 2) of state `x` to assess convergence, stop when either a condition on that list is satisfied or when max # of iterations is reached  
loopUntilAcc :: Int -> ([t] -> Bool) -> (t -> t)  -> t -> t
loopUntilAcc nitermax q f x = go 0 [] x where
  go i ll xx | length ll < 2 = go (i + 1) (y : ll) y 
             | otherwise = if q ll || i == nitermax
                           then xx
                           else go (i + 1) (take 2 $ y:ll) y
                where y = f xx

-- | Keep a moving window buffer (length 2) of state `x` to assess convergence, stop when either a condition on that list is satisfied or when max # of iterations is reached (runs in State monad)
modifyInspectN ::
  MonadState s m =>
    Int ->           -- iteration budget
    ([s] -> Bool) -> -- convergence criterion
    (s -> s) ->      -- state stepping function
    m s
modifyInspectN nitermax q f 
  | nitermax > 0 = go 0 []
  | otherwise = error "modifyInspectN : n must be > 0" where
      go i ll = do
        x <- get
        let y = f x
        if length ll < 2
          then do put y
                  go (i + 1) (y : ll)
          else if q ll || i == nitermax
               then do put y
                       return y
               else do put y
                       go (i + 1) (take 2 $ y : ll)


-- helper functions for estimating convergence
meanl :: (Foldable t, Fractional a) => t a -> a
meanl xx = 1/fromIntegral (length xx) * sum xx

norm2l :: (Foldable t, Functor t, Floating a) => t a -> a
norm2l xx = sqrt $ sum (fmap (**2) xx)

diffSqL :: Floating a => [a] -> a
diffSqL xx = (x1 - x2)**2 where [x1, x2] = [head xx, xx!!1]







-- | iterate until convergence is verified or we run out of a fixed iteration budget
untilConverged :: MonadState a m => (a -> SpVector Double) -> (a -> a) -> m a
untilConverged fproj = modifyInspectN 100 (normDiffConverged fproj)

-- | convergence check (FIXME)
normDiffConverged :: (Foldable t, Functor t) =>
     (a -> SpVector Double) -> t a -> Bool
normDiffConverged fp xx = normSq (foldrMap fp (^-^) (zeroSV 0) xx) <= eps


  



-- | run `niter` iterations and append the state `x` to a list `xs`, stop when either the `xs` satisfies a predicate `q` or when the counter reaches 0
runAppendN :: ([t] -> Bool) -> (t -> t) -> Int -> t -> [t]
runAppendN qq ff niter x0 | niter<0 = error "runAppendN : niter must be > 0"
                          | otherwise = go qq ff niter x0 [] where
  go q f n z xs = 
    let x = f z in
    if n <= 0 || q xs then xs
                      else go q f (n-1) x (x : xs)

-- | ", NO convergence check 
runAppendN' :: (t -> t) -> Int -> t -> [t]
runAppendN' ff niter x0 | niter<0 = error "runAppendN : niter must be > 0"
                        | otherwise = go ff niter x0 [] where
  go f n z xs = 
    let x = f z in
    if n <= 0 then xs
              else go f (n-1) x (x : xs)

  






-- * Rounding operations

-- | Rounding rule
almostZero, almostOne :: Double -> Bool
almostZero x = abs x <= eps
almostOne x = x >= (1-eps) && x < (1+eps)

withDefault :: (t -> Bool) -> t -> t -> t
withDefault q d x | q x = d
                  | otherwise = x

roundZero, roundOne :: Double -> Double
roundZero = withDefault almostZero 0
roundOne = withDefault almostOne 1

with2Defaults :: (t -> Bool) -> (t -> Bool) -> t -> t -> t -> t
with2Defaults q1 q2 d1 d2 x | q1 x = d1
                            | q2 x = d2
                            | otherwise = x

-- | Round to respectively 0 or 1 within some predefined numerical precision eps
roundZeroOne :: Double -> Double
roundZeroOne = with2Defaults almostZero almostOne 0 1







-- * Random matrices and vectors

-- |Dense SpMatrix
randMat :: PrimMonad m => Int -> m (SpMatrix Double)
randMat n = do
  g <- MWC.create
  aav <- replicateM (n^2) (MWC.normal 0 1 g)
  let ii_ = [0 .. n-1]
      (ix_,iy_) = unzip $ concatMap (zip ii_ . replicate n) ii_
  return $ fromListSM (n,n) $ zip3 ix_ iy_ aav

-- | Dense SpVector  
randVec :: PrimMonad m => Int -> m (SpVector Double)
randVec n = do
  g <- MWC.create
  bv <- replicateM n (MWC.normal 0 1 g)
  let ii_ = [0..n-1]
  return $ fromListSV n $ zip ii_ bv



-- | Sparse SpMatrix
randSpMat :: Int -> Int -> IO (SpMatrix Double)
randSpMat n nsp | nsp > n^2 = error "randSpMat : nsp must be < n^2 "
                | otherwise = do
  g <- MWC.create
  aav <- replicateM nsp (MWC.normal 0 1 g)
  ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  jj <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  return $ fromListSM (n,n) $ zip3 ii jj aav

-- | Sparse SpVector
randSpVec :: Int -> Int -> IO (SpVector Double)
randSpVec n nsp | nsp > n = error "randSpVec : nsp must be < n"
                | otherwise = do
  g <- MWC.create
  aav <- replicateM nsp (MWC.normal 0 1 g)
  ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  return $ fromListSV n $ zip ii aav





-- * Pretty printing


sizeStr :: SpMatrix a -> String
sizeStr sm =
  unwords ["(",show (nrows sm),"rows,",show (ncols sm),"columns ) ,",show nz,"NZ ( sparsity",show sy,")"] where
  (SMInfo nz sy) = infoSM sm 


showNonZero :: (Show a, Num a, Eq a) => a -> String
showNonZero x  = if x == 0 then " " else show x


toDenseRow :: Num a => SpMatrix a -> IM.Key -> [a]
toDenseRow (SM (_,ncol) im) irow =
  fmap (\icol -> im `lookupWD_IM` (irow,icol)) [0..ncol-1]

toDenseRowClip :: (Show a, Num a) => SpMatrix a -> IM.Key -> Int -> String
toDenseRowClip sm irow ncomax
  | ncols sm > ncomax = unwords (map show h) ++  " ... " ++ show t
  | otherwise = show dr
     where dr = toDenseRow sm irow
           h = take (ncomax - 2) dr
           t = last dr

newline :: IO ()
newline = putStrLn ""

printDenseSM :: (Show t, Num t) => SpMatrix t -> IO ()
printDenseSM sm = do
  newline
  putStrLn $ sizeStr sm
  newline
  printDenseSM' sm 5 5
  newline
  where    
    printDenseSM' :: (Show t, Num t) => SpMatrix t -> Int -> Int -> IO ()
    printDenseSM' sm'@(SM (nr,_) _) nromax ncomax = mapM_ putStrLn rr_' where
      rr_ = map (\i -> toDenseRowClip sm' i ncomax) [0..nr - 1]
      rr_' | nrows sm > nromax = take (nromax - 2) rr_ ++ [" ... "] ++[last rr_]
           | otherwise = rr_


toDenseListClip :: (Show a, Num a) => SpVector a -> Int -> String
toDenseListClip sv ncomax
  | dim sv > ncomax = unwords (map show h) ++  " ... " ++ show t
  | otherwise = show dr
     where dr = toDenseListSV sv
           h = take (ncomax - 2) dr
           t = last dr

printDenseSV :: (Show t, Num t) => SpVector t -> IO ()
printDenseSV sv = do
  newline
  printDenseSV' sv 5
  newline where
    printDenseSV' v nco = putStrLn rr_' where
      rr_ = toDenseListClip v nco :: String
      rr_' | dim sv > nco = unwords [take (nco - 2) rr_ , " ... " , [last rr_]]
           | otherwise = rr_

-- ** Pretty printer typeclass
class PrintDense a where
  prd :: a -> IO ()

instance (Show a, Num a) => PrintDense (SpVector a) where
  prd = printDenseSV

instance (Show a, Num a) => PrintDense (SpMatrix a) where
  prd = printDenseSM



  





-- * Misc. utilities

-- | integer-indexed ziplist
denseIxArray :: [b] -> [(Int, b)]
denseIxArray xs = zip [0..length xs-1] xs 

-- | ", 2d arrays
denseIxArray2 :: Int -> [c] -> [(Int, Int, c)]
denseIxArray2 m xs = zip3 (concat $ replicate n ii_) jj_ xs where
  ii_ = [0 .. m-1]
  jj_ = concatMap (replicate m) [0 .. n-1]
  ln = length xs
  n = ln `div` m


-- folds

-- | foldr over the results of a fmap
foldrMap :: (Foldable t, Functor t) => (a -> b) -> (b -> c -> c) -> c -> t a -> c
foldrMap ff gg x0 = foldr gg x0 . fmap ff

-- | strict left fold
foldlStrict :: (a -> b -> a) -> a -> [b] -> a
foldlStrict f = go
  where
    go z []     = z
    go z (x:xs) = let z' = f z x in z' `seq` go z' xs

-- | indexed right fold
ifoldr :: Num i =>
     (a -> b -> b) -> b -> (i -> c -> d -> a) -> c -> [d] -> b  
ifoldr mjoin mneutral f  = go 0 where
  go i z (x:xs) = mjoin (f i z x) (go (i+1) z xs)
  go _ _ [] = mneutral


-- *** Bounds checking
type LB = Int
type UB = Int

inBounds :: LB -> UB -> Int -> Bool
inBounds ibl ibu i = i>= ibl && i<ibu

inBounds2 :: (LB, UB) -> (Int, Int) -> Bool
inBounds2 (ibl,ibu) (ix,iy) = inBounds ibl ibu ix && inBounds ibl ibu iy


-- ", lower bound = 0
inBounds0 :: UB -> Int -> Bool
inBounds0 = inBounds 0

inBounds02 :: (UB, UB) -> (Int, Int) -> Bool
inBounds02 (bx,by) (i,j) = inBounds0 bx i && inBounds0 by j




-- * Type synonyms for SpMatrix
type Rows = Int
type Cols = Int

type IxRow = Int
type IxCol = Int












  
