{-# LANGUAGE FlexibleContexts #-}
{-# language TemplateHaskell #-}
module Math.Linear.Sparse where

import Control.Monad (mapM_, forM_, replicateM)

-- import Control.Monad.State.Class
import Control.Monad.State

import Control.Lens hiding ((#))

-- import Data.Functor ((<$>))
import qualified Data.IntMap as IM
import qualified Data.Foldable as F
import qualified Data.List as L
import Data.Maybe

-- additive ring 
class Functor f => Additive f where
  -- | zero element
  zero :: Num a => f a
  
  -- | componentwise operations
  (^+^) :: Num a => f a -> f a -> f a
  (^-^) :: Num a => f a -> f a -> f a

  -- | linear interpolation. Doesn't really belong here though
  
  -- lerp :: Num a => a -> f a -> f a -> f a

  -- union binary lift
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- intersection binary lift
  liftI2 :: (a -> b -> c) -> f a -> f b -> f c


negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate


class Additive f => VectorSpace f where
  (.*) :: Num a => a -> f a -> f a

class VectorSpace f => Normed f where
  dot :: Num a => f a -> f a -> a
  -- norm :: Num a => a -> f a -> a


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






-- Sparse Matrices
data SpMatrix a = SM {smDim :: (Int, Int),
                      smData :: IM.IntMap (SpVector a)} deriving Eq
nrows = fst . smDim
ncols = snd . smDim

                  
mkSpMatrix :: (Int, Int) -> SpMatrix a
mkSpMatrix d = SM d IM.empty



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
spMatrixFromList d xx = go xx (mkSpMatrix d) where
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

  
-- -- testing testing

v1 :: SpVector Int
v1 = SV 5 (IM.fromList [(0, 4), (1, 3)])
v2 :: SpVector Double
v2 = SV 4 (IM.fromList [(0, 2), (1, 3.2), (3, 15)])

m1 = mkSpMatrix (3,4)
m2 = insertSpMatrix 1 3  pi m1
m3 = insertSpMatrix 3 4 1 m2



-- smInsert (ii,jj) x (SM (nrows, ncols))


normSq v = v `dot` v

-- | BiCSSTAB
-- solve A x = b

-- initial residual
residual aa b x0 = b ^-^ (aa #> x0)

converged :: SpMatrix Double -> SpVector Double -> SpVector Double -> Bool
converged aa b x0 = normSq (residual aa b x0) <= eps

eps = 1e-8

bicgIter aa b r0 r0hat rim rhoim alphaim omegaim pim vim xim
  | normSq hres <= eps = h
  | normSq xires <= eps = xi
  | otherwise = ri
  where
  rhoi = r0 `dot` rim
  beta = rhoi * alphaim /(rhoim * omegaim)
  pii = rim ^+^ (beta .* (pim ^-^ (omegaim .* vim)))
  vi = aa #> pii
  alphai = rhoi / (r0hat `dot` vi)
  h = xim ^+^ (alphai .* vi)
  hres = (aa #> h) ^-^ b  -- residual of candidate solution
  s = rim ^-^ (alphai .* vi)
  t = aa #> s
  omegai = (t `dot` s) / (t `dot` t)
  xi = h ^+^ (omegai .* s)
  xires = (aa #> xi) ^-^ b -- residual of second candidate soln
  ri = s ^-^ (omegai .* t)

bicgStep :: SpMatrix Double
                  -> SpVector Double
                  -> SpVector Double
                  -> SpVector Double
                  -> BICG
                  -> BICG
bicgStep aa b r0 r0hat (BICG rim rhoim alphaim omegaim pim vim xim) =
  BICG ri rhoi alphai omegai pii vi xi
  where
  rhoi = r0 `dot` rim
  beta = rhoi * alphaim /(rhoim * omegaim)
  pii = rim ^+^ (beta .* (pim ^-^ (omegaim .* vim)))
  vi = aa #> pii
  alphai = rhoi / (r0hat `dot` vi)
  h = xim ^+^ (alphai .* vi)
  hres = (aa #> h) ^-^ b  -- residual of candidate solution
  s = rim ^-^ (alphai .* vi)
  t = aa #> s
  omegai = (t `dot` s) / (t `dot` t)
  xi = h ^+^ (omegai .* s)
  xires = (aa #> xi) ^-^ b -- residual of second candidate soln
  ri = s ^-^ (omegai .* t)

bicgStepState
  :: MonadState BICG m =>
     SpMatrix Double
     -> SpVector Double -> SpVector Double -> SpVector Double -> m ()
bicgStepState aa b r0 r0hat = modify (bicgStep aa b r0 r0hat)


bicgNSteps
  :: SpMatrix Double
     -> SpVector Double
     -> SpVector Double
     -> SpVector Double
     -> Int
     -> BICG
     -> BICG
bicgNSteps aa b r0 r0hat n =
  execState $ replicateM n (bicgStepState aa b r0 r0hat)



bicgIter'' :: MonadState BICG m =>
           SpMatrix Double ->   -- matrix
           SpVector Double ->   -- rhs
           SpVector Double ->   -- initial 
           SpVector Double ->
           m (SpVector Double)
bicgIter'' aa b r0 r0hat = do
  (BICG rim rhoim alphaim omegaim pim vim xim) <- get
  let
    rhoi = r0 `dot` rim
    beta = rhoi * alphaim /(rhoim * omegaim)
    pii = rim ^+^ (beta .* (pim ^-^ (omegaim .* vim)))
    vi = aa #> pii
    alphai = rhoi / (r0hat `dot` vi)
    h = xim ^+^ (alphai .* vi)
    hres = (aa #> h) ^-^ b  -- residual of candidate solution
    s = rim ^-^ (alphai .* vi)
    t = aa #> s
    omegai = (t `dot` s) / (t `dot` t)
    xi = h ^+^ (omegai .* s)
    xires = (aa #> xi) ^-^ b -- residual of second candidate soln
    ri = s ^-^ (omegai .* t)
  put $ BICG ri rhoi alphai omegai pii vi xi  
  if normSq hres <= eps
     then return h
     else if normSq xires <= eps
          then return xi
          else return ri


-- _aa :: SpMatrix Double,    -- matrix
-- _b :: SpVector Double,     -- rhs
-- _r0 :: SpVector Double,    -- initial residual
-- _r0hat :: SpVector Double, -- candidate solution: r0hat `dot` r0 >= 0

data BICG =
  BICG { _rim :: SpVector Double,   -- r_(i-1)
         _rhoim :: Double,          -- rho_(i-1)
         _alphaim :: Double,        -- alpha_(i-1)
         _omegaim :: Double,        -- omega_(i-1)
         _pim :: SpVector Double,   -- p_(i-1)
         _vim :: SpVector Double,   -- v_(i-1)
         _xim :: SpVector Double    -- x_(i-1)
       } deriving (Eq, Show)

makeLenses ''BICG





zeroSV (SV n _) = SV n IM.empty


bicgsState0 aa b x0 =
  BICG
    (b ^-^ (aa #> x0)) 1 1 1 v0 p0 x0 where
      v0 = zeroSV x0
      p0 = zeroSV x0

-- n iterations of BiCGSTAB
bicgsSolveNsteps
  :: SpMatrix Double
     -> SpVector Double
     -> SpVector Double
     -> SpVector Double
     -> Int
     -> SpVector Double
bicgsSolveNsteps aa b x0 x0hat n =
  _xim $ bicgNSteps aa b x0 x0hat n (bicgsState0 aa b x0)



-- -- testing testing

