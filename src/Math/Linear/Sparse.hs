{-# LANGUAGE FlexibleContexts #-}

module Math.Linear.Sparse where


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

-- | Normed vector space
class VectorSpace f => Normed f where
  -- | inner product
  dot :: Num a => f a -> f a -> a
  -- norm :: Num a => a -> f a -> a



-- | =======================================================

-- | IntMap implementation (can be swapped out with different backends in case)
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
  
instance Normed IM.IntMap where
   a `dot` b = sum $ liftI2 (*) a b 

normSq :: (Normed f, Num a) => f a -> a
normSq v = v `dot` v

   

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
  n .* v = fmap (*n) v

instance Normed SpVector where
  sv1 `dot` sv2
    | d1 == d2 = dot (svData sv1) (svData sv2)
    | otherwise = error $  "dot : vector sizes must coincide (instead : "++ show (d1, d2) ++ ")" where
        (d1, d2) = (svDim sv1, svDim sv2)

  

-- | empty sparse vector (size n, no entries)

zeroSV :: Int -> SpVector a
zeroSV n = SV n IM.empty

-- | create a sparse vector from an association list while discarding all zero entries
mkSpVector :: (Num a, Eq a) => Int -> IM.IntMap a -> SpVector a
mkSpVector d im = SV d $ IM.filterWithKey (\k v -> v /= 0 && inBounds0 d k) im

-- | ", from logically dense array (consecutive indices)
mkSpVectorD :: (Num a, Eq a) => Int -> [a] -> SpVector a
mkSpVectorD d ll = mkSpVector d (IM.fromList $ denseIxArray (take d ll))




insertSpVector :: Int -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim)
  | inBounds0 d i = SV d (IM.insert i x xim)
  | otherwise = error "insertSpVector : index out of bounds"


fromListSV :: Int -> [(Int, a)] -> SpVector a
fromListSV d iix = SV d (IM.fromList (filter (inBounds0 d . fst) iix ))

toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d im) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]


  
instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (IM.toList x)

lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV i (SV _ im) = IM.findWithDefault 0 i im 

findWithDefault0IM :: Num a => IM.Key -> IM.IntMap a -> a
findWithDefault0IM = IM.findWithDefault 0



    
                      



-- | =======================================================

-- | Sparse Matrices
-- data SpMatrix a = SM {smDim :: (Int, Int),
--                       smData :: IM.IntMap (SpVector a)} deriving Eq

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



-- | ========= MATRIX METADATA

-- -- predicates
validIxSM :: SpMatrix a -> (Int, Int) -> Bool
validIxSM mm = inBounds02 (dimSM mm)

isSquareSM :: SpMatrix a -> Bool
isSquareSM m = nrows m == ncols m


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






-- | ========= IntMap-of-IntMap (IM2) stuff


-- insert an element
insertIM2 ::
  IM.Key -> IM.Key -> a -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
insertIM2 i j x imm = IM.insert i (IM.insert j x ro) imm where
  ro = maybe (IM.singleton j x) (IM.insert j x) (IM.lookup i imm)

-- lookup a key
lookupIM2 ::
  IM.Key -> IM.Key -> IM.IntMap (IM.IntMap a) -> Maybe a
lookupIM2 i j imm = IM.lookup i imm >>= IM.lookup j

-- populate an IM2 from a list of (row index, column index, value)  
fromListIM2 ::
  Foldable t =>
     t (IM.Key, IM.Key, a) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
fromListIM2 iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertIM2 i j x t

-- indexed fold over an IM2
ifoldlIM2 ::
  (IM.Key -> IM.Key -> t -> IM.IntMap a -> IM.IntMap a) ->
  IM.IntMap (IM.IntMap t) ->  
  IM.IntMap a
ifoldlIM2 f m         = IM.foldlWithKey' accRow IM.empty m where
  accRow    acc i row = IM.foldlWithKey' (accElem i) acc row
  accElem i acc j x   = f i j x acc

-- transposeIM2 : inner indices become outer ones and vice versa. No loss of information because both inner and outer IntMaps are nubbed.
transposeIM2 :: IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
transposeIM2 = ifoldlIM2 (flip insertIM2)


-- map over outer IM and filter all inner IM's
ifilterIM2 ::
  (IM.Key -> IM.Key -> a -> Bool) ->
  IM.IntMap (IM.IntMap a) ->
  IM.IntMap (IM.IntMap a)
ifilterIM2 f  =
  IM.mapWithKey (\irow row -> IM.filterWithKey (f irow) row) 


-- map over IM2

mapIM2 :: (a -> b) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap b)
mapIM2 = IM.map . IM.map   -- imapIM2 (\_ _ x -> f x)


-- indexed map over IM2
imapIM2 ::
  (IM.Key -> IM.Key -> a -> b) ->
  IM.IntMap (IM.IntMap a) ->
  IM.IntMap (IM.IntMap b)
imapIM2 f im = IM.mapWithKey ff im where
  ff j x = IM.mapWithKey (`f` j) x




-- map over a single `column`

mapColumnIM2 :: (b -> b) -> IM.IntMap (IM.IntMap b) -> Int -> IM.IntMap (IM.IntMap b)
mapColumnIM2 f im jj = imapIM2 (\i j x -> if j == jj then f x else x) im


  

build g = g (:) []

toList t = build (\c n -> foldr c n t)












-- | ========= SPARSE MATRIX BUILDERS

zeroSM :: Int -> Int -> SpMatrix a
zeroSM m n = SM (m,n) IM.empty 


insertSpMatrix :: Int -> Int -> a -> SpMatrix a -> SpMatrix a
insertSpMatrix i j x s
  | inBounds02 d (i,j) = SM d $ insertIM2 i j x smd 
  | otherwise = error "insertSpMatrix : index out of bounds" where
      smd = immSM s
      d = dimSM s


fromListSM' :: Foldable t => t (Int, Int, a) -> SpMatrix a -> SpMatrix a
fromListSM' iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertSpMatrix i j x t

fromListSM :: Foldable t => (Int, Int) -> t (Int, Int, a) -> SpMatrix a  
fromListSM (m,n) iix = fromListSM' iix (zeroSM m n)


fromListDenseSM :: Int -> [a] -> SpMatrix a
fromListDenseSM m ll = fromListSM (m, n) $ denseIxArray2 m ll where
  n = length ll `div` m
  


-- -- diagonal and identity matrix
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

-- unsafe : no bounds checking
extractRowSMU :: SpMatrix a -> IxRow -> SpMatrix a
extractRowSMU s irow =
  SM (1, ncols s) $ IM.filterWithKey (\i _ -> irow == i) (immSM s)

extractRowsSMU :: SpMatrix a -> LB -> UB -> SpMatrix a
extractRowsSMU s i1 i2 =
  SM (1, ncols s) $ IM.filterWithKey (\i _ -> inBounds i1 i2 i) (immSM s)

extractColSMU :: SpMatrix a -> IxCol -> SpMatrix a
extractColSMU s jcol = SM d' imm' where
  imm' = IM.map ff (immSM s)
  d' = (nrows s, 1)
  ff = IM.filterWithKey (\j _ -> j==jcol)

  
-- safe : out-of-bounds indexing is thrown as an exception
-- extractRowsSM :: SpMatrix a -> IxRow -> Int -> SpMatrix a
extractRowsSM (SM (nro,nco) im) i1 i2
  | inBounds0 nro i1  && inBounds0 nro i2 && i2 >= i1 = SM (i2-i1,nco) imf
  | otherwise = error $ "rowsSM : invalid indexing " ++ show (i1, i2) where
      imf = IM.filterWithKey (\i _ -> inBounds i1 i2 i) im










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










  



-- | ========= ALGEBRAIC PRIMITIVE OPERATIONS


-- | transpose


transposeSM :: SpMatrix a -> SpMatrix a
transposeSM (SM (m, n) im) = SM (n, m) (transposeIM2 im)





-- | mapping 




-- | folding

-- count sub-diagonal nonzeros



countSubdiagonalNZ :: IM.IntMap (IM.IntMap a) -> Int
countSubdiagonalNZ im =
  IM.size $ IM.filter (not . IM.null) (ifilterIM2 (\i j _ -> i>j) im)

countSubdiagonalNZSM :: SpMatrix a -> Int
countSubdiagonalNZSM (SM _ im) = countSubdiagonalNZ im
  



-- | sparsify : remove 0s (!!!)

sparsifyIM2 :: IM.IntMap (IM.IntMap Double) -> IM.IntMap (IM.IntMap Double)
sparsifyIM2 = ifilterIM2 (\_ _ x -> x /= 0.0)

sparsifySM :: SpMatrix Double -> SpMatrix Double
sparsifySM (SM d im) = SM d $ sparsifyIM2 im



-- | ROUNDING operations (!!!)
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

roundZeroOneSM :: SpMatrix Double -> SpMatrix Double
roundZeroOneSM (SM d im) = sparsifySM $ SM d $ mapIM2 roundZeroOne im




-- | matrix action on a vector

{- 
FIXME : matVec is more generic than SpVector's :

\m v -> fmap (`dot` v) m
  :: (Normed f1, Num b, Functor f) => f (f1 b) -> f1 b -> f b
-}

matVec, (#>) :: Num a => SpMatrix a -> SpVector a -> SpVector a
matVec (SM (nrows,_) mdata) (SV n sv) = SV nrows $ fmap (`dot` sv) mdata

(#>) = matVec


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








-- | ========= predicates

isOrthogonalSM :: SpMatrix Double -> Bool
isOrthogonalSM sm@(SM (_,n) _) = rsm == eye n where
  rsm = roundZeroOneSM $ transposeSM sm ## sm












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

givens :: SpMatrix Double -> Int -> Int -> SpMatrix Double
givens mm i j 
  | validIxSM mm (i,j) && isSquareSM mm =
       fromListSM' [(i,i,c),(j,j,c),(j,i,-s),(i,j,s)] (eye (nrows mm))
  | otherwise = error "givens : indices out of bounds"      
  where
    (c, s, _) = givensCoef a b
    a = mm @@ (i-1,j) -- FIXME to be chosen from column of b, below diagonal
    b = mm @@ (i,j)   -- element to zero out






{-
Givens method, row version: choose other row index i' s.t. i' is :
* below the diagonal
* corresponding element is nonzero
-} 

-- chooseNZix m (i,j) nr
--   | nr < nrows m - j = undefined
--    where
--      i_ = [j .. j + nr - 1] -- indices below the diagonal
--      elems = IM.filter

-- filterMaybe q ll
--   | null ll' = Nothing
--   | otherwise = Just ll' where
--   ll' = L.filter q ll

       
{-
%%%%Van Loan's Function, Chapter 7%%%%%%%%
  function [c,s] = GivensRotation(x1,x2);
% Pre:
%   x1,x2   scalars
% Post:
%   c,s     c^2+s^2=1 so -s*x1 + c*x2 = 0.
%
   if x2==0
      c = 1;
          s = 0;
   else
      if abs(x2)>=abs(x1)
             cotangent = x1/x2;
                 s = 1/sqrt(1+cotangent^2);
                 c = s*cotangent;
          else
             tangent = x2/x1;
                 c = 1/sqrt(1+tangent^2);
                 s = c*tangent;
          end
   end

-}



-- | ========= QR algorithm





{-
applies Givens rotation iteratively to zero out sub-diagonal elements
-}







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
  | otherwise = case method of
      CGS_ ->  _xBicgstab (bicgstab aa b x0 x0)
      BICGSTAB_ -> _x (cgs aa b x0 x0)
     where
      x0 = mkSpVectorD n $ replicate n 0.1 
      (m,n) = dimSM aa
      mb = dimSV b

-- <\> : sets default solver method 

(<\>) :: SpMatrix Double -> SpVector Double -> SpVector Double      
(<\>) = linSolve BICGSTAB_
  











--







-- | Show details and contents of sparse matrix

sizeStr :: SpMatrix a -> String
sizeStr sm =
  unwords ["(",show (nrows sm),"rows,",show (ncols sm),"columns ) ,",show nz,"NZ (sparsity",show spy,")"] where
  (SMInfo nz spy) = infoSM sm 





-- -- showSparseMatrix :: (Show α, Eq α, Num α) => [[α]] -> String
-- showSparseMatrix [] = "(0,0):\n[]\n"
-- showSparseMatrix m = show (length m, length (head m))++": \n"++
--     (unlines $ L.map (("["++) . (++"]") . L.intercalate "|")
--              $ L.transpose $ L.map column $ L.transpose m)

-- column :: (Show a, Num a, Eq a) => [a] -> [[Char]]
-- column c = let c'       = L.map showNonZero c
--                width    = L.maximum $ L.map length c'
--                offset x = replicate (width - (length x)) ' ' ++ x
--            in L.map offset c'

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


toDenseListClip :: (Show a, Num a) => SpVector a -> Int -> [Char]
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
    printDenseSV' v@(SV nv vd) nco = putStrLn rr_' where
      rr_ = toDenseListClip v nco :: String
      rr_' | svDim sv > nco = unwords [take (nco - 2) rr_ , " ... " , [last rr_]]
           | otherwise = rr_

class PrintDense a where
  prd :: a -> IO ()

instance (Show a, Num a) => PrintDense (SpVector a) where
  prd = printDenseSV

instance (Show a, Num a) => PrintDense (SpMatrix a) where
  prd = printDenseSM










-- | utilities

-- transform state until a condition is met

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

tm0, tm1, tm2 :: SpMatrix Double
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


  
