{-# language TypeFamilies, ScopedTypeVariables #-}
{-# language UndecidableInstances #-}
module Data.Array.Accelerate.Sparse.COOElem where

import Foreign.Storable (Storable(..), peek, poke, alignment, sizeOf)
import Foreign.Ptr
import qualified Data.Vector.Storable as VS

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Array.Data

import Data.Typeable

-- | Matrix element in COO form
newtype COOElem i e = CooE (i, i, e) deriving (Eq, Show)

getRow, getCol :: COOElem i e -> i
getRow (CooE (i, _, _)) = i
getCol (CooE (_, j, _)) = j

getElem :: COOElem i e -> e
getElem (CooE (_, _, e)) = e


-- | Lexicographic ordering of matrix elements, rows-first
instance (Ord i, Eq e) => Ord (COOElem i e) where
  (CooE (i, j, _)) <= (CooE (i', j', _))
    | i < i' = True
    | i == i' && j <= j' = True
    | otherwise = False

instance (Integral i, Storable i, Storable e) => VS.Storable (COOElem i e) where
  sizeOf _ = sizeOf (undefined :: e) + 2 * sizeOf (undefined :: i)
  alignment _ = max (alignment (undefined :: e)) (alignment (undefined :: i))


type instance EltRepr (COOElem i a) = EltRepr (i, i, a)

instance (ArrayElt (EltRepr e), Typeable (EltRepr e), ArrayElt (EltRepr i), Typeable (EltRepr i), Show i, Show e, Typeable i, Typeable e, Elt i, Elt e) => A.Elt (COOElem i e) where
  fromElt (CooE (i, j, x)) = fromElt (i, j, x)
  toElt c = let (i, j, x) = toElt c in CooE (i, j, x)
  eltType (_ :: COOElem i a) = eltType (undefined :: (i, i, a))

-- -- 

-- instance (Elt a, Elt b, Elt c) => Elt (a, b, c) where
--   eltType _             = PairTuple (eltType (undefined :: (a, b))) (eltType (undefined :: c))
--   fromElt (a, b, c)     = (fromElt (a, b), fromElt c)
--   toElt (ab, c)         = let (a, b) = toElt ab in (a, b, toElt c)

-- instance (Show e, Show i, Typeable i, Typeable e, Typeable (EltRepr i), Typeable (EltRepr e), ArrayElt (EltRepr i), ArrayElt (EltRepr e)) => A.Elt (COOElem i e) where
  -- eltType _ = 
