{-# language TypeFamilies, ScopedTypeVariables, MultiParamTypeClasses #-}
{-# language FlexibleInstances, UndecidableInstances #-}
module Data.Array.Accelerate.Sparse.COOElem where

import Foreign.Storable (Storable(..))
import Foreign.Ptr
-- import qualified Data.Vector.Storable as VS

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Smart
import Data.Array.Accelerate.Product  ( IsProduct(..), TupleIdx(..) )

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

-- | COOElem must be Storable if we want to sort vectors of these via vector-algorithms
instance (Integral i, Storable i, Storable e) => Storable (COOElem i e) where
  sizeOf _ = sizeOf (undefined :: e) + 2 * sizeOf (undefined :: i)
  alignment _ = max (alignment (undefined :: e)) (alignment (undefined :: i))
  peek p = do
    let szint = sizeOf (undefined :: i)
    i <- peekByteOff p 0
    j <- peekByteOff p szint
    let p' = castPtr p
    e <- peekByteOff p' szint
    return $ CooE (i, j, e)
  poke p (CooE (i,j,e)) = do
    let szint = sizeOf (undefined :: i)
    pokeByteOff p 0 i
    pokeByteOff p szint j
    let p' = castPtr p
    pokeByteOff p' szint e
    
    
  


-- * `accelerate`-related instances  

type instance EltRepr (COOElem i a) = EltRepr (i, i, a)

-- | Elt
instance (ArrayElt (EltRepr e), Typeable (EltRepr e), ArrayElt (EltRepr i), Typeable (EltRepr i), Show i, Show e, Typeable i, Typeable e, Elt i, Elt e) => A.Elt (COOElem i e) where
  fromElt (CooE (i, j, x)) = fromElt (i, j, x)
  toElt c = let (i, j, x) = toElt c in CooE (i, j, x)
  eltType (_ :: COOElem i a) = eltType (undefined :: (i, i, a))

-- | IsProduct
instance (Elt ix, Elt e) => IsProduct Elt (COOElem ix e) where
  type ProdRepr (COOElem ix e) = ((((), ix), ix), e)
  fromProd _ (CooE (i, j, x)) = ((((), i), j), x)
  toProd _ ((((),i), j), x) = CooE (i,j,x)
  prod cst _ = prod cst (undefined :: (i,i,e))

-- | Lift
instance (A.Lift A.Exp e, Elt (A.Plain e), A.Lift A.Exp i, Elt (A.Plain i)) => A.Lift A.Exp (COOElem i e) where
  type Plain (COOElem i e) = COOElem (A.Plain i) (A.Plain e)
  lift (CooE (i,j,x)) = Exp . Tuple $ NilTup `SnocTup` A.lift i `SnocTup` A.lift j `SnocTup` A.lift x

-- | Unlift
instance (A.Lift A.Exp e, Elt (A.Plain e), A.Lift A.Exp i, Elt (A.Plain i), Elt i, Elt e) => A.Unlift A.Exp (COOElem (A.Exp i) (A.Exp e)) where
  unlift c = let ii = Exp $ SuccTupIdx (SuccTupIdx ZeroTupIdx) `Prj` c
                 jj = Exp $ SuccTupIdx ZeroTupIdx `Prj` c
                 xx = Exp $ ZeroTupIdx `Prj` c
             in CooE (ii, jj, xx)



