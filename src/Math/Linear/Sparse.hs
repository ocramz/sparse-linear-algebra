{-# LANGUAGE FlexibleContexts #-}

module Math.Linear.Sparse where

import Control.Monad (mapM_, forM_, replicateM)
import Control.Monad.Loops

import Control.Monad.State

import qualified Data.IntMap as IM
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

  -- union binary lift
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- intersection binary lift
  liftI2 :: (a -> b -> c) -> f a -> f b -> f c


negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate



-- | Vector space
class Additive f => VectorSpace f where
  (.*) :: Num a => a -> f a -> f a

-- linear interpolation
lerp :: (VectorSpace f, Num a) => a -> f a -> f a -> f a
lerp a u v = a .* u ^+^ ((1-a) .* v)

-- | Normed vector space
class VectorSpace f => Normed f where
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

data SpVector a = SV { svDim :: Int ,
                       svData :: IM.IntMap a} deriving Eq

emptySVector :: Int -> SpVector a
emptySVector n = SV n IM.empty

-- | create a sparse vector from an association list while discarding all zero entries
mkSpVector :: (Num a, Eq a) => Int -> IM.IntMap a -> SpVector a
mkSpVector d im = SV d $ IM.filterWithKey (\k v -> v /= 0 && inBounds k d) im

-- | ", from logically dense array (consecutive indices)
mkSpVectorD :: (Num a, Eq a) => Int -> [a] -> SpVector a
mkSpVectorD d ll = mkSpVector d (IM.fromList $ ixArray (take d ll))


ixArray :: [b] -> [(Int, b)]
ixArray xs = zip [0..length xs-1] xs 


insertSpVector :: Int -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim)
  | inBounds i d = SV d (IM.insert i x xim)
  | otherwise = error "insertSpVector : index out of bounds"

instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (snd.unzip.IM.toList $ x)

lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV i (SV _ im) = IM.findWithDefault 0 i im 

findWithDefault0IM :: Num a => IM.Key -> IM.IntMap a -> a
findWithDefault0IM i = IM.findWithDefault 0 i

toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d im) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]

    
                      
-- instances for SparseVector
instance Functor SpVector where
  fmap f (SV n x) = SV n (fmap f x)

-- | fold functions are applied to non-zero values
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
  sv1 `dot` sv2 = dot (svData sv1) (svData sv2)




-- | =======================================================

-- Sparse Matrices
data SpMatrix a = SM {smDim :: (Int, Int),
                      smData :: IM.IntMap (SpVector a)} deriving Eq
nrows = fst . smDim
ncols = snd . smDim

                  
emptySpMatrix :: (Int, Int) -> SpMatrix a
emptySpMatrix d = SM d IM.empty



-- instance Show a => Show (SpMatrix a) where
--   show (SM d x) = "SM " ++ show d ++ " "++ fmap show (IM.toList x)


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
                                         

-- | Looks up an element in the matrix (if not found, zero is returned)
-- (#) :: (Num a) => SpMatrix a -> (Int,Int) -> a
-- m # (i,j) = maybe 0 (IM.findWithDefault 0 j) (IM.lookup i (smData m))


inBounds :: (Ord a, Num a) => a -> a -> Bool
inBounds i b = i>=0 && i<=b
inBounds2 :: (Ord a, Num a) => (a, a) -> (a, a) -> Bool
inBounds2 (i,j) (bx,by) = inBounds i bx && inBounds j by


insertSpMatrix :: Int -> Int -> a -> SpMatrix a -> SpMatrix a
insertSpMatrix i j x (SM dims smd)
  | inBounds2 (i,j) dims = SM dims (IM.insert i (insertSpVector j x ri) smd) 
  | otherwise = error "insertSpMatrix : index out of bounds" where
      (dx, dy) = dims
      ri = fromMaybe (emptySVector dy) (IM.lookup i smd)


-- indexSpM i j (SM _ im) a = maybe 0 (IM.findWithDefault 0 j) (IM.lookup i (svData im))





inB1 i d s f
  | inBounds i d = f
  | otherwise = error s

inB2 i d s f
  | inBounds2 i d = f
  | otherwise = error s  


lookupSM i j (SM d im) =
  IM.lookup i im >>= \(SV nv ve) -> case IM.lookup j ve of Just c -> return c
                                                           Nothing -> return 0


lookupColSV j (SV d im) = IM.lookup j im

lookupSM' i j (SM d im) = IM.lookup i im >>= \(SV nv ve) ->
  IM.findWithDefault 0 j ve


-- spMatrixFromList d xs = fmap (\(ii, x) -> insertSpMatrix ii x (mkSpMatrix d)) xs 

spMatrixFromList :: Monad m => (Int, Int) -> [(Int, Int, t)] -> m (SpMatrix t)
spMatrixFromList d xx = go xx (emptySpMatrix d) where
  go ((i,j,x):xs) mat = do
    let mat' = insertSpMatrix i j x mat
    go xs mat'
  go [] m = return m
  
  

  

instance Functor SpMatrix where
  fmap f (SM d md) = SM d ((fmap . fmap) f md)



instance Additive SpMatrix where
  zero = SM (0,0) IM.empty
  (^+^) = liftU2 (+)
  (^-^) = liftU2 (-)
  liftU2 f2 (SM n1 x1) (SM n2 x2) = SM (maxTup n1 n2) ((liftU2.liftU2) f2 x1 x2)
  liftI2 f2 (SM n1 x1) (SM n2 x2) = SM (minTup n1 n2) ((liftI2.liftI2) f2 x1 x2)

maxTup, minTup :: Ord t => (t, t) -> (t, t) -> (t, t)
maxTup (x1,y1) (x2,y2) = (max x1 x2, max y1 y2)
minTup (x1,y1) (x2,y2) = (min x1 x2, min y1 y2)


matVec, (#>) :: Num a => SpMatrix a -> SpVector a -> SpVector a
matVec (SM (nrows,_) mdata) sv = SV nrows $ fmap (`dot` sv) mdata

(#>) = matVec

  




-- | LINEAR SOLVERS : solve A x = b

eps = 1e-8

-- initial residua
residual :: Num a => SpMatrix a -> SpVector a -> SpVector a -> SpVector a
residual aa b x0 = b ^-^ (aa #> x0)

converged :: SpMatrix Double -> SpVector Double -> SpVector Double -> Bool
converged aa b x0 = normSq (residual aa b x0) <= eps



-- | CGS

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

cgsN ::
  SpMatrix Double ->
  SpVector Double ->
  SpVector Double ->
  SpVector Double ->
  Int ->
  CGS
cgsN aa b x0 rhat n =
  execState (replicateM n (modify (cgsStep aa rhat))) cgsInit where
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


bicgstabN aa b x0 r0hat n =
  execState (replicateM n (modify (bicgstabStep aa r0hat))) bicgsInit where
   r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
   p0 = r0
   bicgsInit = BICGSTAB x0 r0 p0

instance Show BICGSTAB where
  show (BICGSTAB x r p) = "x = " ++ show x ++ "\n" ++
                                "r = " ++ show r ++ "\n" ++
                                "p = " ++ show p ++ "\n"







-- | utilities
zeroSV (SV n _) = SV n IM.empty










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
  

replicateSwitch p m f = loop m where
      loop n | n <= 0 || p = pure ()
             | otherwise = f *> loop (n-1)

-- untilM p 
               
