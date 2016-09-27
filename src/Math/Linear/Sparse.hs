{-# LANGUAGE FlexibleContexts #-}

module Math.Linear.Sparse where

import Control.Monad (mapM_, forM_, replicateM)
import Control.Monad.Loops

import Control.Monad.State

import qualified Data.IntMap as IM
-- import Data.Utils.StrictFold (foldlStrict) -- hidden in `containers`

import qualified Data.Foldable as F
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
emptySVector :: Int -> SpVector a
emptySVector n = SV n IM.empty

zeroSV :: Int -> SpVector a
zeroSV n = SV n IM.empty

-- | create a sparse vector from an association list while discarding all zero entries
mkSpVector :: (Num a, Eq a) => Int -> IM.IntMap a -> SpVector a
mkSpVector d im = SV d $ IM.filterWithKey (\k v -> v /= 0 && inBounds0 d k) im

-- | ", from logically dense array (consecutive indices)
mkSpVectorD :: (Num a, Eq a) => Int -> [a] -> SpVector a
mkSpVectorD d ll = mkSpVector d (IM.fromList $ ixArray (take d ll))

-- | integer-indexed ziplist
ixArray :: [b] -> [(Int, b)]
ixArray xs = zip [0..length xs-1] xs 


insertSpVector :: Int -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim)
  | inBounds0 d i = SV d (IM.insert i x xim)
  | otherwise = error "insertSpVector : index out of bounds"

-- fromListSpUnsafe d iix = SV d (IM.fromList iix)

fromListSV :: Int -> [(Int, a)] -> SpVector a
fromListSV d iix = SV d (IM.fromList (filter (inBounds0 d . fst) iix )) 
 
instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (IM.toList x)

lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV i (SV _ im) = IM.findWithDefault 0 i im 

findWithDefault0IM :: Num a => IM.Key -> IM.IntMap a -> a
findWithDefault0IM i = IM.findWithDefault 0 i

toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d im) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]

    
                      





-- | =======================================================

-- | Sparse Matrices
-- data SpMatrix a = SM {smDim :: (Int, Int),
--                       smData :: IM.IntMap (SpVector a)} deriving Eq

data SpMatrix a = SM {smDim :: (Int, Int),
                     smData :: IM.IntMap (IM.IntMap a)} deriving Eq

-- | instances for SpMatrix
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


-- | nrows, ncols : size accessors
nrows, ncols :: SpMatrix a -> Int
nrows = fst . smDim
ncols = snd . smDim


data SMInfo = SMInfo { smNz :: Int,
                       smSpy :: Double} deriving (Eq, Show)

infoSM :: SpMatrix a -> SMInfo
infoSM (SM (nr,nc) im) = SMInfo nz spy where
  nz = IM.foldr (+) 0 $ IM.map IM.size im
  spy = fromIntegral nz / fromIntegral (nr*nc)



-- | ========= DISPLAY

-- | Show details and contents of sparse matrix

sizeStr :: SpMatrix a -> String
sizeStr sm =
  unwords ["SM:",show (nrows sm),"rows,",show (ncols sm),"columns,",show nz,"NZ (sparsity",show spy,")"] where
  (SMInfo nz spy) = infoSM sm 



  

instance Show a => Show (SpMatrix a) where
  show sm@(SM d x) = "SM " ++ sizeStr sm ++ " "++ show (IM.toList x)


-- showSparseMatrix :: (Show α, Eq α, Num α) => [[α]] -> String
showSparseMatrix [] = "(0,0):\n[]\n"
showSparseMatrix m = show (length m, length (head m))++": \n"++
    (unlines $ L.map (("["++) . (++"]") . L.intercalate "|")
             $ L.transpose $ L.map column $ L.transpose m)

column :: (Show a, Num a, Eq a) => [a] -> [[Char]]
column c = let c'       = L.map showNonZero c
               width    = L.maximum $ L.map length c'
               offset x = replicate (width - (length x)) ' ' ++ x
           in L.map offset c'

showNonZero :: (Show a, Num a, Eq a) => a -> [Char]
showNonZero x  = if x == 0 then " " else show x

-- | Converts sparse matrix to plain list-matrix with all zeroes restored
-- fillMx :: (Num α) => SparseMatrix α -> [[α]]
-- fillMx m = [ [ m # (i,j) | j <- [1 .. nrows  m] ]
--                          | i <- [1 .. ncols m] ]
                                         

toDenseRow :: Num a => SpMatrix a -> IM.Key -> [a]
toDenseRow (SM (_,ncol) im) irow =
  fmap (\icol -> im `lookupWithDIM` (irow,icol)) [0..ncol-1]

toDenseRowClip :: (Show a, Num a) => SpMatrix a -> IM.Key -> Int -> String
toDenseRowClip sm irow ncomax
  | ncols sm > ncomax = unwords (map show h) ++  " ... " ++ show t
  | otherwise = show dr
     where dr = toDenseRow sm irow
           h = take (ncomax - 2) dr
           t = last dr

printDenseSM' sm@(SM (nr,nc) im) nromax ncomax = mapM_ putStrLn rr_' where
  rr_ = map (\i -> toDenseRowClip sm i ncomax) [0..nr - 1]
  rr_' | nrows sm > nromax = take (nromax - 2) rr_ ++ [" ... "] ++[last rr_]
       | otherwise = rr_

printDenseSM :: (Show t, Num t) => SpMatrix t -> IO ()
printDenseSM sm = do
  putStrLn ""
  putStrLn $ sizeStr sm
  putStrLn ""
  printDenseSM' sm 5 5
  putStrLn ""

  



-- | ========= BUILDERS

zeroSM :: Int -> Int -> SpMatrix a
zeroSM m n = SM (m,n) IM.empty 



insertSpMatrix :: Int -> Int -> a -> SpMatrix a -> SpMatrix a
insertSpMatrix i j x (SM dims smd)
  | inBounds02 dims (i,j) = SM dims (IM.insert i (IM.insert j x ri) smd) 
  | otherwise = error "insertSpMatrix : index out of bounds" where
      -- (dx, dy) = dims
      ri = fromMaybe IM.empty (IM.lookup i smd)


fromListSM' :: Foldable t => t (Int, Int, a) -> SpMatrix a -> SpMatrix a
fromListSM' iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertSpMatrix i j x t

fromListSM :: Foldable t => (Int, Int) -> t (Int, Int, a) -> SpMatrix a  
fromListSM (m,n) iix = fromListSM' iix (zeroSM m n)
  


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
encode (nr,nc) (i,j) = i + (j * nr)

decode :: (Int, Int) -> Int -> (Rows, Cols)
decode (nr, nc) ci = (r, c) where (c,r ) = quotRem ci nr

type Rows = Int
type Cols = Int
-- newtype Ix = Ix {unIx :: (Rows, Cols)} deriving Eq
-- instance Show Ix where show (Ix ii) = show ii






-- | ========= SUB-MATRICES



-- rowsSM (SM d im) i1 i2
--   | inBounds2 d (i1, i2) = IM.filterWithKey (\i x -> inBounds)




-- | ========= LOOKUP

lookupSM :: SpMatrix a -> IM.Key -> IM.Key -> Maybe a
lookupSM (SM d im) i j = IM.lookup i im >>= IM.lookup j

-- | Looks up an element in the matrix (if not found, zero is returned)

lookupWithDefaultSpm, (#) :: Num a => SpMatrix a -> (IM.Key, IM.Key) -> a
lookupWithDefaultSpm (SM d m) (i,j) =
  fromMaybe 0 (IM.lookup i m >>= IM.lookup j)

lookupWithDIM :: Num a => IM.IntMap (IM.IntMap a) -> (IM.Key, IM.Key) -> a
lookupWithDIM im (i,j) = fromMaybe 0 (IM.lookup i im >>= IM.lookup j)

(#) = lookupWithDefaultSpm




-- FIXME : to throw an exception or just ignore the out-of-bound access ?

-- inB1 i d s f
--   | inBounds0 i d = f
--   | otherwise = error s

-- inB2 i d s f
--   | inBounds02 i d = f
--   | otherwise = error s  








  
  

  



-- | ========= ALGEBRAIC PRIMITIVE OPERATIONS

-- | matrix action on a vector

{- 
FIXME : matVec is more generic than SpVector's :

\m v -> fmap (`dot` v) m
  :: (Normed f1, Num b, Functor f) => f (f1 b) -> f1 b -> f b
-}

matVec, (#>) :: Num a => SpMatrix a -> SpVector a -> SpVector a
matVec (SM (nrows,_) mdata) (SV n sv) = SV nrows $ fmap (`dot` sv) mdata

(#>) = matVec

  


-- | matrix-matrix product

matMat :: Num a => SpMatrix a -> SpMatrix a -> SpMatrix a
matMat (SM (nr1,nc1) m1) (SM (nr2,nc2) m2)
  | nc1 == nr2 = SM (nr1, nc2) (fmap (\vm1 -> fmap (`dot` vm1) m2) m1)
  | otherwise = error "matMat : incompatible matrix sizes"




-- | diagonal and identity matrices






-- | Givens rotation matrix

data Givens = Givens !Double !Double deriving (Eq, Show)

-- -- givens :: Floating a => Int -> Int -> Int -> a -> SpMatrix a
-- givens n i j theta
--   | inBounds2 (i,j) (n,n) =
--        fromListSM' [(i,i,c),(j,j,c),(j,i,-s),(i,j,s)] (eye n)
--   | otherwise = error "givens : indices out of bounds"      
--   where
--    -- c = cos theta
--    -- s = sin theta
--     -- x1 = mm 
--     Givens c s 
--       | x2 == 0 = Givens 1 0 
--       | otherwise = let cot = x1/x2
--                         ta = x2/x1
--                         s1 = (1/sqrt (1 + cot**2))
--                         c1 = s1*cot
--                         c2 = 1/sqrt (1 + ta**2)
--                         s2 = c2 * ta
--                     in
--                     if abs x2 >= abs x1 then Givens c1 s1 
--                                         else Givens c2 s2
                  


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





-- | =======================================================

-- | LINEAR SOLVERS : solve A x = b

-- | numerical tolerance for e.g. solution convergence
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
  where
    normDiffConverged fp xx = normSq (foldrMap (^-^) (zeroSV 0) fp xx) <= eps
              

foldrMap :: (Foldable t, Functor t) => (a -> c -> c) -> c -> (a1 -> a) -> t a1 -> c
foldrMap ff x0 pp = foldr ff x0 . fmap pp





-- inspectBicgstabStates :: [BICGSTAB] -> Bool
-- inspectBicgstabStates = withProjs (\[x0,x1] -> normSq (x0 ^-^ x1) <= eps) _xBicgstab
 


-- gox n f xs = 


-- misc

foldlStrict :: (a -> b -> a) -> a -> [b] -> a
foldlStrict f = go
  where
    go z []     = z
    go z (x:xs) = let z' = f z x in z' `seq` go z' xs



-- bounds checking

type LB = Int
type UB = Int

inBounds :: LB -> UB -> Int -> Bool
inBounds ibl ibu i = i>= ibl && i<=ibu

inBounds2 :: (LB, UB) -> (Int, Int) -> Bool
inBounds2 (ibl,ibu) (ix,iy) = inBounds ibl ibu ix && inBounds ibl ibu iy

inBounds0 :: UB -> Int -> Bool
inBounds0 = inBounds 0


-- inBounds0 :: (Ord a, Num a) => a -> a -> Bool
-- inBounds0 i b = i>=0 && i<=b

-- inBounds02 :: (Ord a, Num a) => (a, a) -> (a, a) -> Bool
inBounds02 (i,j) (bx,by) = inBounds0 i bx && inBounds0 j by




-- playground

-- | terminate after n iterations or when q becomes true, whichever comes first
untilC :: (a -> Bool) -> Int ->  (a -> a) -> a -> a
untilC p n f = go n
  where
    go m x | p x || m <= 0 = x
           | otherwise     = go (m-1) (f x)





-- testing State


data T0 = T0 {unT :: Int} deriving Eq
instance Show T0 where
  show (T0 x) = show x

-- modifyT :: MonadState T0 m => (Int -> Int) -> m String
modifyT f = state (\(T0 i) -> (i, T0 (f i)))
  

t00 = T0 0

testT n = execState $ replicateM n (modifyT (+1)) 


-- testT2 = execState $ when 
  

-- replicateSwitch p m f = loop m where
--       loop n | n <= 0 || p = pure (#)
--              | otherwise = f *> loop (n-1)


               
