{-# LANGUAGE FlexibleContexts #-}

module Math.Linear.Sparse where

import Math.Linear.Sparse.IntMap 


import Control.Monad.Primitive

import Control.Monad (mapM_, forM_, replicateM)
import Control.Monad.Loops

import Control.Monad.State
import Control.Monad.Writer
import Control.Monad.Trans

import Control.Monad.Trans.State (runStateT)
import Control.Monad.Trans.Writer (runWriterT)

import qualified Data.IntMap.Strict as IM
-- import Data.Utils.StrictFold (foldlStrict) -- hidden in `containers`

import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWC

import Data.Monoid
import qualified Data.Foldable as F
import qualified Data.Traversable as T


import qualified Data.List as L
import Data.Maybe



-- | Additive ring 
class Functor f => Additive f where
  -- | zero element
  zero :: Num a => f a
  
  -- | componentwise operations
  (^+^) :: Num a => f a -> f a -> f a
  (^-^) :: Num a => f a -> f a -> f a

  -- |union binary lift
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- |intersection binary lift
  liftI2 :: (a -> b -> c) -> f a -> f b -> f c

-- | negate the values in a functor
negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate



-- | Vector space
class Additive f => VectorSpace f where
  -- | multiplication by a scalar
  (.*) :: Num a => a -> f a -> f a
  

-- |linear interpolation
lerp :: (VectorSpace f, Num a) => a -> f a -> f a -> f a
lerp a u v = a .* u ^+^ ((1-a) .* v)


-- | Hilbert space (inner product)
class VectorSpace f => Hilbert f where
  -- | inner product
  dot :: Num a => f a -> f a -> a

class Hilbert f => Normed f where
  norm :: (Floating a, Eq a) => a -> f a -> a


-- some norms and related results

-- squared norm 
normSq :: (Hilbert f, Num a) => f a -> a
normSq v = v `dot` v


-- L1 norm
norm1 :: (Foldable t, Num a, Functor t) => t a -> a
norm1 v = sum (fmap abs v)

-- Euclidean norm
norm2 :: (Hilbert f, Floating a) => f a -> a
norm2 v = sqrt (normSq v)

-- Lp norm (p > 0)
normP :: (Foldable t, Functor t, Floating a) => a -> t a -> a
normP p v = sum u**(1/p) where
  u = fmap (**p) v

-- infinity-norm
normInfty :: (Foldable t, Ord a) => t a -> a
normInfty = maximum



-- normalize
normalize :: (Normed f, Floating a, Eq a) => a -> f a -> f a
normalize n v = (1 / norm n v) .* v




-- -- Lp inner product (p > 0)
dotLp :: (Additive t, Foldable t, Floating a) => a -> t a -> t a ->  a
dotLp p v1 v2 = sum u**(1/p) where
  f a b = (a*b)**p
  u = liftI2 f v1 v2


-- reciprocal
reciprocal :: (Functor f, Fractional b) => f b -> f b
reciprocal = fmap (\x -> 1 / x)


-- scale
scale :: (Num b, Functor f) => b -> f b -> f b
scale n = fmap (* n)



-- | =======================================================

-- | IntMap implementation
instance Additive IM.IntMap where
  zero = IM.empty
  {-# INLINE zero #-}
  liftU2 = IM.unionWith
  {-# INLINE liftU2 #-}
  liftI2 = IM.intersectionWith
  {-# INLINE liftI2 #-}
  (^+^) = liftU2 (+)
  {-# INLINE (^+^) #-}
  x ^-^ y = x ^+^ negated y
  {-# INLINE (^-^) #-}

instance VectorSpace IM.IntMap where
  n .* im = IM.map (* n) im
  
instance Hilbert IM.IntMap where
   a `dot` b = sum $ liftI2 (*) a b 

instance Normed IM.IntMap where
  norm p v | p==1 = norm1 v
           | p==2 = norm2 v
           | otherwise = normP p v



-- | =======================================================

-- | Sparse Vector
data SpVector a = SV { svDim :: Int ,
                       svData :: IM.IntMap a} deriving Eq

dimSV :: SpVector a -> Int
dimSV = svDim

-- internal : projection functions, do not export
imSV :: SpVector a -> IM.IntMap a
imSV = svData


-- | instances for SpVector
instance Functor SpVector where
  fmap f (SV n x) = SV n (fmap f x)

instance Foldable SpVector where
    foldr f d v = F.foldr f d (svData v)

instance Additive SpVector where
  zero = SV 0 IM.empty
  (^+^) = liftU2 (+)
  (^-^) = liftU2 (-)
  liftU2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftU2 f2 x1 x2)
  liftI2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftI2 f2 x1 x2)
                      
instance VectorSpace SpVector where
  n .* v = scale n v

instance Hilbert SpVector where
  sv1 `dot` sv2
    | d1 == d2 = dot (svData sv1) (svData sv2)
    | otherwise = error $  "dot : vector sizes must coincide (instead : "++ show (d1, d2) ++ ")" where
        (d1, d2) = (svDim sv1, svDim sv2)

instance Normed SpVector where
  norm p (SV _ v) = norm p v
  



-- | empty sparse vector (size n, no entries)

zeroSV :: Int -> SpVector a
zeroSV n = SV n IM.empty


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

mkSpVector1D :: Int -> [a] -> SpVector a
mkSpVector1D d ll = mkSpVector1 d (IM.fromList $ denseIxArray (take d ll))



-- vector of `1`s
onesSV :: Num a => Int -> SpVector a
onesSV d = SV d $ IM.fromList $ denseIxArray $ replicate d 1

zerosSV :: Num a => Int -> SpVector a
zerosSV d = SV d $ IM.fromList $ denseIxArray $ replicate d 0




-- insert
insertSpVector :: Int -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim)
  | inBounds0 d i = SV d (IM.insert i x xim)
  | otherwise = error "insertSpVector : index out of bounds"


fromListSV :: Int -> [(Int, a)] -> SpVector a
fromListSV d iix = SV d (IM.fromList (filter (inBounds0 d . fst) iix ))

-- toList
toListSV :: SpVector a -> [(IM.Key, a)]
toListSV sv = IM.toList (imSV sv)

-- to dense list (default = 0)
toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d im) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]









  
instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (IM.toList x)


-- | lookup

lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV i (SV _ im) = IM.findWithDefault 0 i im 

findWithDefault0IM :: Num a => IM.Key -> IM.IntMap a -> a
findWithDefault0IM = IM.findWithDefault 0




-- | SV manipulation

tailSV :: SpVector a -> SpVector a
tailSV (SV n sv) = SV (n-1) ta where
  ta = IM.mapKeys (\i -> i - 1) $ IM.delete 0 sv
  

headSV :: Num a => SpVector a -> a
headSV sv = fromMaybe 0 (IM.lookup 0 (imSV sv))



-- | concatenate SpVector


concatSV :: SpVector a -> SpVector a -> SpVector a
concatSV (SV n1 s1) (SV n2 s2) = SV (n1+n2) (IM.union s1 s2') where
  s2' = IM.mapKeys (+ n1) s2










-- | promote a SV to SM

svToSM :: SpVector a -> SpMatrix a
svToSM (SV n d) = SM (n, 1) $ IM.singleton 0 d



    

-- | outer vector product

outerProdSV, (><) :: Num a => SpVector a -> SpVector a -> SpMatrix a
outerProdSV v1 v2 = fromListSM (m, n) ixy where
  m = svDim v1
  n = svDim v2
  ixy = [(i,j, x * y) | (i,x) <- toListSV v1 , (j, y) <- toListSV v2]

(><) = outerProdSV









-- | =======================================================


data SpMatrix a = SM {smDim :: (Int, Int),
                     smData :: IM.IntMap (IM.IntMap a)} deriving Eq



-- | instances for SpMatrix
instance Show a => Show (SpMatrix a) where
  show sm@(SM _ x) = "SM " ++ sizeStr sm ++ " "++ show (IM.toList x)

instance Functor SpMatrix where
  fmap f (SM d md) = SM d ((fmap . fmap) f md)

instance Additive SpMatrix where
  zero = SM (0,0) IM.empty
  (^+^) = liftU2 (+)
  (^-^) = liftU2 (-)
  liftU2 f2 (SM n1 x1) (SM n2 x2) = SM (maxTup n1 n2) ((liftU2.liftU2) f2 x1 x2)
  liftI2 f2 (SM n1 x1) (SM n2 x2) = SM (minTup n1 n2) ((liftI2.liftI2) f2 x1 x2)

-- | TODO : use semilattice properties instead
maxTup, minTup :: Ord t => (t, t) -> (t, t) -> (t, t)
maxTup (x1,y1) (x2,y2) = (max x1 x2, max y1 y2)
minTup (x1,y1) (x2,y2) = (min x1 x2, min y1 y2)

-- | empty matrix of size d
emptySpMatrix :: (Int, Int) -> SpMatrix a
emptySpMatrix d = SM d IM.empty



-- multiply matrix by a scalar
matScale :: Num a => a -> SpMatrix a -> SpMatrix a
matScale a = fmap (*a)

-- Frobenius norm (sqrt of trace of M^T M)
normFrobenius :: SpMatrix Double -> Double
normFrobenius m = sqrt $ foldlSM (+) 0 m' where
  m' | nrows m > ncols m = transposeSM m ## m
     | otherwise = m ## transposeSM m 
  





-- | ========= MATRIX METADATA

-- -- predicates
-- are the supplied indices within matrix bounds?
validIxSM :: SpMatrix a -> (Int, Int) -> Bool
validIxSM mm = inBounds02 (dimSM mm)

-- is the matrix square?
isSquareSM :: SpMatrix a -> Bool
isSquareSM m = nrows m == ncols m

-- is the matrix diagonal?
isDiagonalSM :: SpMatrix a -> Bool
isDiagonalSM m = IM.size d == nrows m where
  d = IM.filterWithKey ff (immSM m)
  ff irow row = IM.size row == 1 &&
                IM.size (IM.filterWithKey (\j e -> j == irow) row) == 1








-- -- internal projection functions, do not export:
immSM :: SpMatrix t -> IM.IntMap (IM.IntMap t)
immSM (SM _ imm) = imm

dimSM :: SpMatrix t -> (Int, Int)
dimSM (SM d _) = d

nelSM :: SpMatrix t -> Int
nelSM (SM (nr,nc) _) = nr*nc

-- | nrows, ncols : size accessors
nrows, ncols :: SpMatrix a -> Int
nrows = fst . smDim
ncols = snd . smDim




data SMInfo = SMInfo { smNz :: Int,
                       smSpy :: Double} deriving (Eq, Show)

infoSM :: SpMatrix a -> SMInfo
infoSM s = SMInfo nz spy where
  nz = IM.foldr (+) 0 $ IM.map IM.size (immSM s)
  spy = fromIntegral nz / fromIntegral (nelSM s)


-- # NZ in row i

nzRowU :: SpMatrix a -> IM.Key -> Int
nzRowU s i = maybe 0 IM.size (IM.lookup i $ immSM s)

nzRow :: SpMatrix a -> IM.Key -> Int
nzRow s i | inBounds0 (nrows s) i = nzRowU s i
          | otherwise = error "nzRow : index out of bounds"










-- | ========= SPARSE MATRIX BUILDERS

zeroSM :: Int -> Int -> SpMatrix a
zeroSM m n = SM (m,n) IM.empty 


insertSpMatrix :: Int -> Int -> a -> SpMatrix a -> SpMatrix a
insertSpMatrix i j x s
  | inBounds02 d (i,j) = SM d $ insertIM2 i j x smd 
  | otherwise = error "insertSpMatrix : index out of bounds" where
      smd = immSM s
      d = dimSM s


-- | from list (row, col, value)
fromListSM' :: Foldable t => t (Int, Int, a) -> SpMatrix a -> SpMatrix a
fromListSM' iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertSpMatrix i j x t

fromListSM :: Foldable t => (Int, Int) -> t (Int, Int, a) -> SpMatrix a  
fromListSM (m,n) iix = fromListSM' iix (zeroSM m n)


fromListDenseSM :: Int -> [a] -> SpMatrix a
fromListDenseSM m ll = fromListSM (m, n) $ denseIxArray2 m ll where
  n = length ll `div` m
  


-- | to List

-- toDenseListSM : populate missing entries with 0
toDenseListSM :: Num t => SpMatrix t -> [(Int, Int, t)]
toDenseListSM m =
  [(i, j, m @@ (i, j)) | i <- [0 .. nrows m - 1], j <- [0 .. ncols m- 1]]





-- -- create diagonal and identity matrix
mkDiagonal :: Int -> [a] -> SpMatrix a
mkDiagonal n xx = fromListSM (n,n) $ zip3 ii ii xx where
  ii = [0..n-1]


eye :: Num a => Int -> SpMatrix a
eye n = mkDiagonal n (ones n)

ones :: Num a => Int -> [a]
ones n = replicate n 1
  

-- fromList :: [(Key,a)] -> IntMap a
-- fromList xs
--   = foldlStrict ins empty xs
--   where
--     ins t (k,x)  = insert k x t




encode :: (Int, Int) -> (Rows, Cols) -> Int
encode (nr,_) (i,j) = i + (j * nr)

decode :: (Int, Int) -> Int -> (Rows, Cols)
decode (nr, _) ci = (r, c) where (c,r ) = quotRem ci nr

type Rows = Int
type Cols = Int
-- newtype Ix = Ix {unIx :: (Rows, Cols)} deriving Eq
-- instance Show Ix where show (Ix ii) = show ii

type IxRow = Int
type IxCol = Int







-- | ========= SUB-MATRICES


extractSubmatrixSM :: SpMatrix a -> (Int, Int) -> (Int, Int) -> SpMatrix a
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

-- extract row / column
extractRowSM :: SpMatrix a -> Int -> SpMatrix a
extractRowSM sm i = extractSubmatrixSM sm (i, i) (0, ncols sm - 1)

extractColSM :: SpMatrix a -> Int -> SpMatrix a
extractColSM sm j = extractSubmatrixSM sm (0, nrows sm - 1) (j, j)


toSV :: SpMatrix a -> SpVector a
toSV (SM (m,n) im) = SV d $ snd . head $ IM.toList im where
  d | m==1 && n>1 = n 
    | n==1 && m>1 = m
    | otherwise = error $ "toSV : incompatible dimensions " ++ show (m,n)


-- extract row or column and place into SpVector
extractCol :: SpMatrix a -> Int -> SpVector a
extractCol m i = toSV $ extractColSM m i

extractRow :: SpMatrix a -> Int -> SpVector a
extractRow m j = toSV $ extractRowSM m j







-- | ========= MATRIX STACKING

vertStackSM, (-=-) :: SpMatrix a -> SpMatrix a -> SpMatrix a
vertStackSM mm1 mm2 = SM (m, n) $ IM.union u1 u2 where
  nro1 = nrows mm1
  m = nro1 + nrows mm2
  n = max (ncols mm1) (ncols mm2)
  u1 = immSM mm1
  u2 = IM.mapKeys (+ nro1) (immSM mm2)

(-=-) = vertStackSM


horizStackSM, (-||-) :: SpMatrix a -> SpMatrix a -> SpMatrix a
horizStackSM mm1 mm2 = t (t mm1 -=- t mm2) where
  t = transposeSM

(-||-) = horizStackSM







-- | ========= LOOKUP

lookupSM :: SpMatrix a -> IM.Key -> IM.Key -> Maybe a
lookupSM (SM d im) i j = IM.lookup i im >>= IM.lookup j

-- | Looks up an element in the matrix (if not found, zero is returned)

lookupWD_SM, (@@) :: Num a => SpMatrix a -> (IM.Key, IM.Key) -> a
lookupWD_SM (SM d m) (i,j) =
  fromMaybe 0 (IM.lookup i m >>= IM.lookup j)

lookupWD_IM :: Num a => IM.IntMap (IM.IntMap a) -> (IM.Key, IM.Key) -> a
lookupWD_IM im (i,j) = fromMaybe 0 (IM.lookup i im >>= IM.lookup j)

(@@) = lookupWD_SM




-- FIXME : to throw an exception or just ignore the out-of-bound access ?

-- inB1 i d s f
--   | inBounds0 i d = f
--   | otherwise = error s

-- inB2 i d s f
--   | inBounds02 i d = f
--   | otherwise = error s  






-- | ========= MISC SpMatrix OPERATIONS

foldlSM :: (a -> b -> b) -> b -> SpMatrix a -> b
foldlSM f n (SM _ m)= foldlIM2 f n m

ifoldlSM :: (IM.Key -> IM.Key -> a -> b -> b) -> b -> SpMatrix a -> b
ifoldlSM f n (SM _ m) = ifoldlIM2' f n m






-- | mapping 





-- | folding

-- count sub-diagonal nonzeros
countSubdiagonalNZSM :: SpMatrix a -> Int
countSubdiagonalNZSM (SM _ im) = countSubdiagonalNZ im


-- | filtering

-- extractDiagonalSM :: (Num a, Eq a) => SpMatrix a -> SpVector a
-- extractDiagonalSM (SM (m,n) im) = mkSpVectorD m $ extractDiagonalIM2 im

-- extract with default 0
extractDiagonalDSM :: (Num a, Eq a) => SpMatrix a -> SpVector a
extractDiagonalDSM mm = mkSpVector1D n $ foldr ins [] ll  where
  ll = [0 .. n - 1]
  n = nrows mm
  ins i acc = mm@@(i,i) : acc
  




  

--  filtering the index subset that lies below the diagonal

subdiagIndicesSM :: SpMatrix a -> [(IM.Key, IM.Key)]
subdiagIndicesSM (SM _ im) = subdiagIndices im





-- | sparsify : remove 0s (!!!)

sparsifyIM2 :: IM.IntMap (IM.IntMap Double) -> IM.IntMap (IM.IntMap Double)
sparsifyIM2 = ifilterIM2 (\_ _ x -> abs x >= eps)

sparsifySM :: SpMatrix Double -> SpMatrix Double
sparsifySM (SM d im) = SM d $ sparsifyIM2 im



-- | ROUNDING operations (!!!)
                              
roundZeroOneSM :: SpMatrix Double -> SpMatrix Double
roundZeroOneSM (SM d im) = sparsifySM $ SM d $ mapIM2 roundZeroOne im







  



-- | ========= ALGEBRAIC PRIMITIVE OPERATIONS


-- | transpose


transposeSM :: SpMatrix a -> SpMatrix a
transposeSM (SM (m, n) im) = SM (n, m) (transposeIM2 im)









-- | matrix action on a vector

{- 
FIXME : matVec is more generic than SpVector's :

\m v -> fmap (`dot` v) m
  :: (Normed f1, Num b, Functor f) => f (f1 b) -> f1 b -> f b
-}

-- matrix on vector
matVec, (#>) :: Num a => SpMatrix a -> SpVector a -> SpVector a
matVec (SM (nrows,_) mdata) (SV n sv) = SV nrows $ fmap (`dot` sv) mdata

(#>) = matVec

-- vector on matrix (FIXME : transposes matrix: more costly than `matVec`)
vecMat, (<#) :: Num a => SpVector a -> SpMatrix a -> SpVector a  
vecMat (SV n sv) (SM (_, ncols) im) = SV ncols $ fmap (`dot` sv) (transposeIM2 im)

(<#) = vecMat










-- | matrix-matrix product

matMat, (##) :: Num a => SpMatrix a -> SpMatrix a -> SpMatrix a
matMat (SM (nr1,nc1) m1) (SM (nr2,nc2) m2)
  | nc1 == nr2 = SM (nr1, nc2) $
      fmap (\vm1 -> fmap (`dot` vm1) (transposeIM2 m2)) m1
  | otherwise = error "matMat : incompatible matrix sizes"

(##) = matMat


-- | sparsified matrix-matrix product (prunes all elements `x` for which `abs x <= eps`)
matMatSparsified, (#~#)  :: SpMatrix Double -> SpMatrix Double -> SpMatrix Double
matMatSparsified m1 m2 = sparsifySM $ matMat m1 m2

(#~#) = matMatSparsified






-- | ========= predicates

-- is the matrix orthogonal? i.e. Q^t ## Q == I
isOrthogonalSM :: SpMatrix Double -> Bool
isOrthogonalSM sm@(SM (_,n) _) = rsm == eye n where
  rsm = roundZeroOneSM $ transposeSM sm ## sm









-- | ========= condition number

conditionNumberSM :: SpMatrix Double -> Double
conditionNumberSM m | isInfinite kappa = error "Infinite condition number : rank-deficient system"
                    | otherwise = kappa where
  kappa = lmax / lmin
  (_, r) = qr m
  u = extractDiagonalDSM r  -- FIXME : need to extract with default element 0 
  lmax = abs (maximum u)
  lmin = abs (minimum u)







-- | ========= Householder transformation

householderMatrix :: Num a => a -> SpVector a -> SpMatrix a
householderMatrix beta x = eye n ^-^ scale beta (x >< x) where
  n = svDim x


-- a vector `x` uniquely defines an orthogonal plane; the Householder operator reflects any point `v` with respect to this plane:
-- v' = (I - 2 x >< x) v 
householderRefl :: SpVector Double -> SpMatrix Double
householderRefl = householderMatrix 2.0











-- | ========= Givens rotation matrix


hypot :: Floating a => a -> a -> a
hypot x y = abs x * (sqrt (1 + y/x)**2)

sign :: (Ord a, Num a) => a -> a
sign x
  | x > 0 = 1
  | x == 0 = 0
  | otherwise = -1 

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


{-
Givens method, row version: choose other row index i' s.t. i' is :
* below the diagonal
* corresponding element is nonzero

QR.C1 ) To zero out entry A(i, j) we must find row k such that A(k, j) is
non-zero but A has zeros in row k for all columns less than j.
-}

givens :: SpMatrix Double -> Int -> Int -> SpMatrix Double
givens mm i j 
  | validIxSM mm (i,j) && isSquareSM mm =
       sparsifySM $ fromListSM' [(i,i,c),(j,j,c),(j,i,-s),(i,j,s)] (eye (nrows mm))
  | otherwise = error "givens : indices out of bounds"      
  where
    (c, s, _) = givensCoef a b
    i' = head $ fromMaybe (error $ "givens: no compatible rows for entry " ++ show (i,j)) (candidateRows (immSM mm) i j)
    a = mm @@ (i', j)
    b = mm @@ (i, j)   -- element to zero out

-- is the `k`th the first nonzero column in the row?
firstNonZeroColumn :: IM.IntMap a -> IM.Key -> Bool
firstNonZeroColumn mm k = isJust (IM.lookup k mm) &&
                          isNothing (IM.lookupLT k mm)

-- returns a set of rows {k} that satisfy QR.C1
candidateRows :: IM.IntMap (IM.IntMap a) -> IM.Key -> IM.Key -> Maybe [IM.Key]
candidateRows mm i j | IM.null u = Nothing
                     | otherwise = Just (IM.keys u) where
  u = IM.filterWithKey (\irow row -> irow /= i &&
                                     firstNonZeroColumn row j) mm





-- | ========= QR algorithm

{-
applies Givens rotation iteratively to zero out sub-diagonal elements
-}

qr :: SpMatrix Double -> (SpMatrix Double, SpMatrix Double)
qr mm = (qmat, rmat)  where
  qmat = transposeSM $ foldr (#~#) ee $ gmats mm -- Q = (G_n * G_n-1 ... * G_1)^T
  rmat = (transposeSM qmat) #~# mm               -- R = Q^T A
  ee = eye (nrows mm)
      
-- Givens matrices in order [GN, G(N-1), .. ]
gmats :: SpMatrix Double -> [SpMatrix Double]
gmats mm = reverse $ gm mm (subdiagIndicesSM mm) where
 gm m ((i,j):is) = let g = givens m i j
                   in g : gm (g #~# m) is
 gm _ [] = []






-- | ========= Householder vector (G & VL Alg.5.1.1)

house :: (Ord a, Floating a) => SpVector a -> (SpVector a, a)
house x = (v, beta) where
  n = svDim x
  tx = tailSV x
  sigma = tx `dot` tx
  vtemp = singletonSV 1 `concatSV` tx
  (v, beta) | sigma == 0 = (vtemp, 0)
            | otherwise = let mu = sqrt (headSV x**2 + sigma)
                              xh = headSV x
                              vh | xh <= 1 = xh - mu
                                 | otherwise = - sigma / (xh + mu)
                              vnew = (1 / vh) .* insertSpVector 0 vh vtemp     
                          in (vnew, 2 * xh**2 / (sigma + vh**2))

                         






-- | ========= SVD

{- Golub & Van Loan, sec 8.6.2 (p 452 segg.)

SVD of A :

* reduce A to upper bidiagonal form B (Alg. 5.4.2)
* compute SVD of B (implicit-shift QR step, Alg. 8.3.2)

-}


















-- | =======================================================

-- | LINEAR SOLVERS : solve A x = b

-- | numerical tolerance for e.g. solution convergence
eps :: Double
eps = 1e-8

-- | residual of candidate solution x0
residual :: Num a => SpMatrix a -> SpVector a -> SpVector a -> SpVector a
residual aa b x0 = b ^-^ (aa #> x0)

converged :: SpMatrix Double -> SpVector Double -> SpVector Double -> Bool
converged aa b x0 = normSq (residual aa b x0) <= eps



-- | CGS

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

cgs ::
  SpMatrix Double ->
  SpVector Double ->
  SpVector Double ->
  SpVector Double ->
  -- Int ->
  CGS
cgs aa b x0 rhat =
  -- execState (replicateM n (modify (cgsStep aa rhat))) cgsInit where
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



  

-- | BiCSSTAB

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


bicgstab
  :: SpMatrix Double
     -> SpVector Double
     -> SpVector Double
     -> SpVector Double
     -- -> Int
     -> BICGSTAB
bicgstab aa b x0 r0hat =
  -- execState (replicateM n (modify (bicgstabStep aa r0hat))) bicgsInit where
  execState (untilConverged _xBicgstab (bicgstabStep aa r0hat)) bicgsInit where
   r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
   p0 = r0
   bicgsInit = BICGSTAB x0 r0 p0
   -- q (BICGSTAB xi _ _) = nor

instance Show BICGSTAB where
  show (BICGSTAB x r p) = "x = " ++ show x ++ "\n" ++
                                "r = " ++ show r ++ "\n" ++
                                "p = " ++ show p ++ "\n"











-- | ========= LINEAR SOLVERS INTERFACE

data LinSolveMethod = CGS_ | BICGSTAB_ deriving (Eq, Show) 

-- random starting vector
linSolveM ::
  PrimMonad m =>
    LinSolveMethod -> SpMatrix Double -> SpVector Double -> m (SpVector Double)
linSolveM method aa b = do
  let (m,n) = dimSM aa
      mb = dimSV b
  if m/=mb then error "linSolve : operand dimensions mismatch" else do
    x0 <- randVec mb
    case method of CGS_ -> return $ _xBicgstab (bicgstab aa b x0 x0)
                   BICGSTAB_ -> return $ _x (cgs aa b x0 x0)

-- deterministic starting vector (every component at 0.1) 
linSolve ::
  LinSolveMethod -> SpMatrix Double -> SpVector Double -> SpVector Double
linSolve method aa b
  | m/=mb = error "linSolve : operand dimensions mismatch"
  | otherwise = solve aa b where
      solve aa' b' | isDiagonalSM aa = (reciprocal aa') #> b'
                   | otherwise = solveWith aa' b' 
      solveWith aa' b' = case method of
                                CGS_ ->  _xBicgstab (bicgstab aa' b' x0 x0)
                                BICGSTAB_ -> _x (cgs aa' b' x0 x0)
      x0 = mkSpVectorD n $ replicate n 0.1 
      (m,n) = dimSM aa
      mb = dimSV b

-- <\> : sets default solver method 

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










-- | ========= PRETTY PRINTING

-- | Show details and contents of sparse matrix

sizeStr :: SpMatrix a -> String
sizeStr sm =
  unwords ["(",show (nrows sm),"rows,",show (ncols sm),"columns ) ,",show nz,"NZ (sparsity",show spy,")"] where
  (SMInfo nz spy) = infoSM sm 


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
  | svDim sv > ncomax = unwords (map show h) ++  " ... " ++ show t
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
      rr_' | svDim sv > nco = unwords [take (nco - 2) rr_ , " ... " , [last rr_]]
           | otherwise = rr_

class PrintDense a where
  prd :: a -> IO ()

instance (Show a, Num a) => PrintDense (SpVector a) where
  prd = printDenseSV

instance (Show a, Num a) => PrintDense (SpMatrix a) where
  prd = printDenseSM










-- | rounding to 0 or 1 within some predefined numerical precision


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

roundZeroOne :: Double -> Double
roundZeroOne = with2Defaults almostZero almostOne 0 1





-- | transform state until a condition is met

modifyUntil :: MonadState s m => (s -> Bool) -> (s -> s) -> m s
modifyUntil q f = do
  x <- get
  let y = f x
  put y
  if q y then return y
         else modifyUntil q f     
  

-- modify state and append, until max # of iterations is reached
modifyInspectN :: MonadState s m => Int -> ([s] -> Bool) -> (s -> s) -> m s
modifyInspectN nn q ff = go nn ff [] where
   go n f xs = do
    x <- get
    if n <= 0
      then do
       put x
       return x
      else do
       let y = f x
           ys = y : xs
       put y
       if (length ys == n) && q ys
         then return y
         else go (n-1) f ys


untilConverged :: MonadState a m => (a -> SpVector Double) -> (a -> a) -> m a
untilConverged fproj = modifyInspectN 2 (normDiffConverged fproj)

-- convergence check (FIXME)
normDiffConverged :: (Foldable t, Functor t) =>
     (a -> SpVector Double) -> t a -> Bool
normDiffConverged fp xx = normSq (foldrMap (^-^) (zeroSV 0) fp xx) <= eps
              


-- run `niter` iterations and append the state `x` to a list `xs`, stop when either the `xs` satisfies a predicate `q` or when the counter reaches 0

runAppendN :: ([t] -> Bool) -> (t -> t) -> Int -> t -> [t]
runAppendN qq ff niter x0 | niter<0 = error "runAppendN : niter must be > 0"
                          | otherwise = go qq ff niter x0 [] where
  go q f n z xs = 
    let x = f z in
    if n <= 0 || q xs then xs
                      else go q f (n-1) x (x : xs)

-- ", NO convergence check 
runAppendN' :: (t -> t) -> Int -> t -> [t]
runAppendN' ff niter x0 | niter<0 = error "runAppendN : niter must be > 0"
                        | otherwise = go ff niter x0 [] where
  go f n z xs = 
    let x = f z in
    if n <= 0 then xs
              else go f (n-1) x (x : xs)

-- runN :: Int -> (a -> a) -> a -> a
-- runN n stepf x0 = runAppendN' stepf n x0
  









-- | random matrices and vectors

-- dense

randMat :: PrimMonad m => Int -> m (SpMatrix Double)
randMat n = do
  g <- MWC.create
  aav <- replicateM (n^2) (MWC.normal 0 1 g)
  let ii_ = [0 .. n-1]
      (ix_,iy_) = unzip $ concatMap (zip ii_ . replicate n) ii_
  return $ fromListSM (n,n) $ zip3 ix_ iy_ aav
  
randVec :: PrimMonad m => Int -> m (SpVector Double)
randVec n = do
  g <- MWC.create
  bv <- replicateM n (MWC.normal 0 1 g)
  let ii_ = [0..n-1]
  return $ fromListSV n $ zip ii_ bv



-- sparse

randSpMat :: Int -> Int -> IO (SpMatrix Double)
randSpMat n nsp | nsp > n^2 = error "randSpMat : nsp must be < n^2 "
                | otherwise = do
  g <- MWC.create
  aav <- replicateM nsp (MWC.normal 0 1 g)
  ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  jj <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  return $ fromListSM (n,n) $ zip3 ii jj aav


randSpVec :: Int -> Int -> IO (SpVector Double)
randSpVec n nsp | nsp > n = error "randSpVec : nsp must be < n"
                | otherwise = do
  g <- MWC.create
  aav <- replicateM nsp (MWC.normal 0 1 g)
  ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  return $ fromListSV n $ zip ii aav






-- | misc utils



-- | integer-indexed ziplist
denseIxArray :: [b] -> [(Int, b)]
denseIxArray xs = zip [0..length xs-1] xs 

-- ", 2d arrays
denseIxArray2 :: Int -> [c] -> [(Int, Int, c)]
denseIxArray2 m xs = zip3 (concat $ replicate n ii_) jj_ xs where
  ii_ = [0 .. m-1]
  jj_ = concatMap (replicate m) [0 .. n-1]
  ln = length xs
  n = ln `div` m


-- folds

foldrMap :: (Foldable t, Functor t) => (a -> c -> c) -> c -> (a1 -> a) -> t a1 -> c
foldrMap ff x0 pp = foldr ff x0 . fmap pp

foldlStrict :: (a -> b -> a) -> a -> [b] -> a
foldlStrict f = go
  where
    go z []     = z
    go z (x:xs) = let z' = f z x in z' `seq` go z' xs

-- indexed fold?
ifoldr :: Num i =>
     (a -> b -> b) -> b -> (i -> c -> d -> a) -> c -> [d] -> b  
ifoldr mjoin mneutral f  = go 0 where
  go i z (x:xs) = mjoin (f i z x) (go (i+1) z xs)
  go _ _ [] = mneutral


-- bounds checking

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









--

tm0, tm1, tm2, tm3, tm4 :: SpMatrix Double
tm0 = fromListSM (2,2) [(0,0,pi), (1,0,sqrt 2), (0,1, exp 1), (1,1,sqrt 5)]

tv0 :: SpVector Double
tv0 = mkSpVectorD 2 [5, 6]

-- wikipedia test matrix for Givens rotation

tm1 = sparsifySM $ fromListDenseSM 3 [6,5,0,5,1,4,0,4,3]

tm1g1 = givens tm1 1 0
tm1a2 = tm1g1 ## tm1

tm1g2 = givens tm1a2 2 1
tm1a3 = tm1g2 ## tm1a2

tm1q = transposeSM (tm1g2 ## tm1g1)


-- wp test matrix for QR decomposition via Givens rotation

tm2 = fromListDenseSM 3 [12, 6, -4, -51, 167, 24, 4, -68, -41]




tm3 = transposeSM $ fromListDenseSM 3 [1 .. 9]

tm3g1 = fromListDenseSM 3 [1, 0,0, 0,c,-s, 0, s, c]
  where c= 0.4961
        s = 0.8682


--

tm4 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]


-- playground

-- | terminate after n iterations or when q becomes true, whichever comes first
untilC :: (a -> Bool) -> Int ->  (a -> a) -> a -> a
untilC p n f = go n
  where
    go m x | p x || m <= 0 = x
           | otherwise     = go (m-1) (f x)





-- testing State


-- data T0 = T0 {unT :: Int} deriving Eq
-- instance Show T0 where
--   show (T0 x) = show x

-- -- modifyT :: MonadState T0 m => (Int -> Int) -> m String
-- modifyT f = state (\(T0 i) -> (i, T0 (f i)))
  

-- t00 = T0 0

-- testT n = execState $ replicateM n (modifyT (+1)) 


-- testT2 = execState $ when 
  

-- replicateSwitch p m f = loop m where
--       loop n | n <= 0 || p = pure (#)
--              | otherwise = f *> loop (n-1)



-- testing Writer
               
-- asdfw n = runWriter $ do
--   tell $ "potato " ++ show n
--   tell "jam"
--   return (n+1)


-- --


-- testing State and Writer



-- runMyApp runA k maxDepth =
--     let config = maxDepth
--         state =  0
--     in runStateT (runWriterT (runA k) config) state


  
